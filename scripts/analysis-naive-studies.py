#! /usr/bin/env python

import os
import re
import csv
import sys
import base64
import hashlib
import requests
import xlsxwriter as xw
from decimal import Decimal
from itertools import groupby
import xml.etree.ElementTree as ET
from collections import OrderedDict, Counter
from analysis_functions import aggregate_naiveseqs_adindex

from lanl_reader import lanl_reader
from codonutils import get_codons, compare_codon, translate_codon
from data_writer import csv_writer, data_writer
from data_reader import (get_prevalence, data_reader,
                         possible_apobecs_reader,
                         ROOT, CONSENSUS)

EFETCH_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
NAIVES = ('Naive', 'PINaive', 'ProbablyNaive')
GENE_RANGE = {
    'gag': (790, 2289),
    'gp41': (7758, 8792)
}
PREC0 = Decimal('1')
PREC1 = Decimal('1.0')

REVIEW_TABLE_HEADERS = [
    'PubID', 'PubMedID', 'PubYear', 'NumPts', 'NumIsolates',
    'NumLANLIsolates', 'NumLANLIsolatesQCPassed', 'Title',
    'Authors', 'Subtypes', 'RxStatus', 'Notes'
]

CLEAN_TABLE_HEADERS = [
    'PubMedID', 'PubYear', 'Title', 'Authors', 'NumPatients'
]


def naive_papers(papers):
    return [p for p in papers if p['RxStatus'] in NAIVES]


def get_fact_table(gene):
    fact_table = data_reader(
        os.path.join(ROOT, 'internalFiles', 'papersReview',
                     '{}GenbankFact.csv'.format(gene.lower())))
    fact_table_merged = {}
    for f in fact_table:
        pubid = f['PubID']
        if pubid not in fact_table_merged:
            fact_table_merged[pubid] = f
        else:
            # TODO: sanity check of the title and authors
            ff = fact_table_merged[pubid]
            if not ff.get('PubYr'):
                ff['PubYr'] = f.get('PubYr')
            if not ff.get('PMID'):
                ff['PMID'] = f.get('PMID')
            if not ff.get('Notes'):
                ff['Notes'] = f.get('Notes')
            if 'NumPts' in ff:
                num_pts = int(ff['NumPts'] or 0)
                num_pts += int(f.get('NumPts') or 0)
                if num_pts:
                    ff['NumPts'] = str(num_pts)
            if 'NumIsolates' in ff:
                num_isos = int(ff['NumIsolates'] or 0)
                num_isos += int(f.get('NumIsolates') or 0)
                if num_isos:
                    ff['NumIsolates'] = str(num_isos)
            old_rx = ff.get('RxStatus')
            new_rx = f.get('RxStatus')
            if old_rx and new_rx and old_rx != new_rx:
                ff['RxStatus'] = 'Conflict'
            if not old_rx:
                ff['RxStatus'] = new_rx

    return fact_table_merged


def qc(gene, sequences):
    results = []
    print('In {} one-per-person aligned {} sequences:'
          .format(len(sequences), gene))
    longgaps = 0
    toomanyfs = 0
    for seq in sequences:
        if ('---' * 10) in seq['AlignedNASequence']:
            longgaps += 1
            continue
        if 'NNN' in seq['AlignedNASequence']:
            longgaps += 1
            continue
        if seq['NumFrameShifts'] > 3:
            toomanyfs += 1
            continue
        results.append(seq)
    print('- {} contained large deletions or '
          'missing nucleotides'.format(longgaps))
    print('- {} contained more than 3 frameshifts'
          .format(toomanyfs))
    return results


def attach_numbers(gene, sequences):
    for sequence in sequences:
        mutations = sequence['Mutations']
        unusuals = 0
        for m in mutations:
            isins = m['IsInsertion']
            isdel = m['IsDeletion']
            if isins:
                aa = 'i'
            elif isdel:
                aa = 'd'
            else:
                aa = m['AminoAcidText']
                cons = m['ReferenceText']
                if len(aa) > 1:
                    aa = aa.replace(cons, '')
                if len(aa) > 1:
                    continue
            prev = get_prevalence(gene, m['Position'], aa)
            if prev < 0.1:
                unusuals += 1

        sequence.update({
            'NumInsertions': len([m for m in mutations if m['IsInsertion']]),
            'NumDeletions': len([m for m in mutations if m['IsDeletion']]),
            'NumStopCodons': len([m for m in mutations
                                  if '*' in m['AminoAcidText']]),
            'NumUnusuals': unusuals
        })
    return sequences


def attach_apobecs(gene, sequences):
    apobecs = {}
    for apobec in possible_apobecs_reader(gene):
        pos = int(apobec['Position'])
        apobecs.setdefault(pos, set()).add(
            apobec['AAChange'].split('=>', 1)[1])
    apobec_cutoffs = {
        'gag': 2,
        'gp41': 2
    }

    toomanyapobecs = 0
    for sequence in sequences:
        mutations = sequence['Mutations']

        sequence['NumApobecs'] = len([
            m for m in mutations
            if set(m['AminoAcidText']) &
            apobecs.get(m['Position'], set())
        ])
        valid_apobecs = sequence['NumApobecs'] <= apobec_cutoffs[gene]
        sequence['Included'] = (
            sequence.get('Included', True) and
            valid_apobecs
        )
        if not valid_apobecs:
            toomanyapobecs += 1
    print('- {} contained {} or more APOBEC mutations'
          .format(toomanyapobecs, apobec_cutoffs[gene]))
    return sequences


def attach_rxstatus(gene, ptseqs):
    fact_table = get_fact_table(gene).values()
    pmids = {p['PMID'] for p in naive_papers(fact_table)}
    if '' in pmids:
        # sequences from unpublished papers should not be all included
        pmids.remove('')
    pubids = {p['PubID']: p['PMID']
              for p in fact_table if p['PMID'] in pmids}
    for p in naive_papers(fact_table):
        # only sequences from naive unpublished papers should be included
        if p['PubID'] not in pubids:
            pubids[p['PubID']] = p['PMID']
    for seq in ptseqs:
        if seq['_PubID'] in pubids:
            seq['RxStatus'] = 'Naive'
            seq['PubMedID'] = pubids[seq['_PubID']]
        else:
            seq['RxStatus'] = 'NA'
            seq['PubMedID'] = ''
    return ptseqs


def export_naive_sequences(gene, ptseqs):
    filename = os.path.join(
        ROOT, 'internalFiles', 'naiveStudies',
        '{}.csv'.format(gene.lower())
    )
    aligned_fasta = os.path.join(
        ROOT, 'data', 'naiveStudies',
        '{}NaiveAligned.fas'.format(gene.lower())
    )
    unaligned_fasta = os.path.join(
        ROOT, 'data', 'naiveStudies',
        '{}NaiveOriginal.fas'.format(gene.lower())
    )
    indels_csv = os.path.join(
        ROOT, 'data', 'naiveStudies',
        '{}NaiveIndels.csv'.format(gene.lower()))

    # remove non-naive sequences
    ptseqs = [seq for seq in ptseqs
              if seq['RxStatus'] == 'Naive' and seq['Included']]
    print('- {} naive {} sequences were found'.format(len(ptseqs), gene))
    genesize = int(CONSENSUS[gene]['Size'])
    siteheaders = ['P{}'.format(i) for i in range(1, genesize + 1)]

    rows = []
    indels = []
    ptseqs = sorted(ptseqs, key=lambda s: s['Accession'])
    for seq in ptseqs:
        firstaa = seq['FirstAA']
        lastaa = seq['LastAA']
        muts = {m['Position']: m for m in seq['Mutations']}
        row = {
            'PMID': seq['PubMedID'],
            'Accession': seq['Accession'],
            'RxStatus': 'Naive',
            'lanlSubtype': seq['Subtype'],
            'NumAAChanges': len(muts),
            'NumInsertions': seq['NumInsertions'],
            'NumDeletions': seq['NumDeletions'],
            'NumStopCodons': seq['NumStopCodons'],
            'NumApobecs': seq['NumApobecs'],
            'NumUnusuals': seq['NumUnusuals'],
            'NumFrameShifts': seq['NumFrameShifts']
        }
        for pos in range(1, genesize + 1):
            pname = 'P{}'.format(pos)
            if pos < firstaa or pos > lastaa:
                row[pname] = '.'
            elif pos not in muts:
                row[pname] = '-'
            else:
                mut = muts[pos]
                if mut['IsInsertion']:
                    row[pname] = 'i'
                elif mut['IsDeletion']:
                    row[pname] = 'd'
                elif mut['IsPartial']:
                    row[pname] = 'X'
                else:
                    aas = mut['AminoAcidText']
                    if len(aas) > 4:
                        aas = 'X'
                    row[pname] = aas
        rows.append(row)
        for mut in seq['Mutations']:
            if mut['IsPartial'] or not (mut['IsInsertion'] or
                                        mut['IsDeletion']):
                continue
            isins = mut['IsInsertion']
            indels.append({
                'Accession': seq['Accession'],
                'Gene': gene,
                'Position': mut['Position'],
                'IndelType': 'ins' if isins else 'del',
                'Codon': mut['CodonText'] if isins else '',
                'InsertedCodons': mut['InsertedCodonsText'] if isins else ''
            })

    csv_writer(
        filename, rows,
        ['PMID', 'Accession', 'RxStatus', 'lanlSubtype',
         'NumAAChanges', 'NumInsertions', 'NumDeletions',
         'NumStopCodons', 'NumApobecs', 'NumUnusuals',
         'NumFrameShifts'] + siteheaders
    )

    csv_writer(
        indels_csv, indels,
        ['Accession', 'Gene', 'Position',
         'IndelType', 'Codon', 'InsertedCodons']
    )

    data_writer(
        aligned_fasta,
        '\n'.join(
            '>{Accession}|{Subtype}\n{AlignedNASequence}'.format(**s)
            for s in ptseqs
        )
    )

    data_writer(
        unaligned_fasta,
        '\n'.join(
            '>{Accession}|{Subtype}\n{NASequence}'.format(**s)
            for s in ptseqs
        )
    )


def create_review_table(gene, ptseqs):
    fact_table = get_fact_table(gene)
    grouped = groupby(ptseqs, lambda s: s['_PubID'])
    results = {}
    for pubid, group_seqs in grouped:
        group_seqs = list(group_seqs)
        subtypes = sorted({s['Subtype'] for s in group_seqs})
        seq = group_seqs[0]
        fact = fact_table.get(pubid, {})
        if fact.get('PubIDCorrection'):
            pubid = fact['PubIDCorrection']

        if pubid not in results:
            results[pubid] = {
                'PubID': pubid,
                'PubMedID': fact.get('PMID', seq['PubMedID']),
                'PubYear': fact.get('PubYr', seq['PubYear']),
                'NumPts': fact.get('NumPts'),
                'NumIsolates': fact.get('NumIsolates'),
                'NumLANLIsolates': len(group_seqs),
                'NumLANLIsolatesQCPassed':
                len([s for s in group_seqs
                     if s['Included']]),
                'Title': seq['Title'],
                'Authors': seq['Authors'],
                'Subtypes': '; '.join(subtypes),
                'RxStatus': fact.get('RxStatus'),
                'Notes': fact.get('Notes'),
            }
        else:
            result = results[pubid]
            origsubtypes = result['Subtypes'].split('; ')
            subtypes = sorted(set(origsubtypes + subtypes))
            result['Subtypes'] = '; '.join(subtypes)
            if not result.get('PubYear'):
                result['PubYear'] = fact.get('PubYr')
            if not result.get('PubMedID'):
                result['PubMedID'] = fact.get('PMID')
            num_pts = int(result['NumPts'] or 0)
            num_pts += int(fact.get('NumPts') or 0)
            if num_pts:
                result['NumPts'] = str(num_pts)
            num_isos = int(result['NumIsolates'] or 0)
            num_isos += int(fact.get('NumIsolates') or 0)
            if num_isos:
                result['NumIsolates'] = str(num_isos)
            result['NumLANLIsolates'] += len(group_seqs)
            result['NumLANLIsolatesQCPassed'] += len([
                s for s in group_seqs if s['Included']
            ])
            if not result.get('RxStatus'):
                result['RxStatus'] = fact.get('RxStatus')

    results = sorted(
        results.values(), key=lambda r: (-r['NumLANLIsolates'], r['PubID']))
    csv_writer(
        os.path.join(
            ROOT, 'internalFiles', 'papersReview',
            '{}ReviewTable.csv'.format(gene)
        ),
        results,
        REVIEW_TABLE_HEADERS)
    export_excel_table(
        os.path.join(
            ROOT, 'internalFiles', 'papersReview',
            '{}ReviewTable.xlsx'.format(gene)
        ),
        results)
    export_naive_papers_table(
        os.path.join(
            ROOT, 'data', 'naiveStudies',
            '{}NaiveStudies.csv'.format(gene)
        ),
        results)


def export_naive_papers_table(filename, rows):
    rows = naive_papers(rows)
    results = []
    for row in rows:
        numpt = row['NumLANLIsolatesQCPassed']
        if not numpt:
            continue
        results.append({
            'PubMedID': row['PubMedID'],
            'PubYear': row['PubYear'],
            'Title': row['Title'],
            'Authors': row['Authors'],
            'NumPatients': numpt
        })
    csv_writer(filename, results, CLEAN_TABLE_HEADERS)


def export_excel_table(filename, rows):
    workbook = xw.Workbook(filename)
    worksheet = workbook.add_worksheet('main')
    worksheet_rx_status = workbook.add_worksheet('validRxStatus')
    fontsize = 14
    workbook.formats[0].set_font_size(fontsize)
    headers = REVIEW_TABLE_HEADERS
    valid_rx_status = sorted([
        'Lab', 'Rx', 'Rx-PI', 'Rx=>Naive', 'Check', 'Naive', 'Unknown',
        'Mixed', 'PINaive', 'ProbablyNaive', 'Unpublished', 'NonM',
        'Entry-Naive', 'Conflict', 'Problem'
    ])
    fmt = workbook.add_format()
    fmt.set_align('top')
    fmt.set_font_size(fontsize)

    numfmt = workbook.add_format()
    numfmt.set_align('top')
    numfmt.set_num_format('#')
    numfmt.set_font_size(fontsize)

    headerfmt = workbook.add_format()
    headerfmt.set_bold()
    headerfmt.set_align('top')
    headerfmt.set_font_size(fontsize)

    wrapfmt = workbook.add_format()
    wrapfmt.set_text_wrap()
    wrapfmt.set_align('top')
    wrapfmt.set_font_size(fontsize)

    for colnum, colname in enumerate(headers):
        worksheet.write(0, colnum, colname)
    for idx, row in enumerate(rows):
        rownum = idx + 1
        for colnum, colname in enumerate(headers):
            value = row[colname]
            if not value:
                continue
            if str(value).isdigit():
                worksheet.write_number(rownum, colnum, int(value))
            else:
                worksheet.write_string(rownum, colnum, value)

    worksheet.set_column('A:A', width=29, cell_format=fmt)  # PubID
    worksheet.set_column('B:B', width=10, cell_format=numfmt)  # PubMedID
    worksheet.set_column('C:C', width=7, cell_format=numfmt)  # PubYear
    worksheet.set_column('D:D', width=7, cell_format=numfmt)  # NumPts
    worksheet.set_column('E:E', width=11, cell_format=numfmt)  # NumIsolates
    worksheet.set_column('F:F', width=15, cell_format=numfmt)  # LANLIsolates
    # LANLIsolatesQCPassed
    worksheet.set_column('G:G', width=15, cell_format=numfmt)
    worksheet.set_column('H:H', width=60, cell_format=wrapfmt)  # Title
    worksheet.set_column('I:I', width=60, cell_format=wrapfmt)  # Authors
    worksheet.set_column('J:J', width=40, cell_format=wrapfmt)  # Subtypes
    worksheet.set_column('K:K', width=15, cell_format=fmt)  # RxStatus
    worksheet.set_row(0, None, cell_format=headerfmt)  # headers

    # add data validation
    for idx, rx_status in enumerate(valid_rx_status):
        worksheet_rx_status.write(idx, 0, rx_status)

    worksheet_rx_status.set_column('A:A', width=15, cell_format=fmt)

    rx_status_col = headers.index('RxStatus')
    worksheet.data_validation(
        1, rx_status_col, len(rows), rx_status_col,
        {'validate': 'list',
         'source': 'validRxStatus!A:A'})

    workbook.window_width = 2400
    workbook.window_height = 1600
    workbook.close()
    print('{} created'.format(filename))


def get_patient_sequences(gene):
    ptseqs = lanl_reader(gene, os.path.join(
        ROOT, 'local', 'hiv-db_{}_squeeze.fasta'.format(gene)
    ))
    ptseqs = [s for s in ptseqs if s['Subtype'] not in 'ONP']
    ptseqs = list(attach_references(gene, ptseqs))
    ptseqs = attach_numbers(gene, ptseqs)
    ptseqs = attach_rxstatus(gene, ptseqs)
    return qc(gene, ptseqs)


def attach_pubid(seq):
    hashed = base64.urlsafe_b64encode(
        hashlib.md5((seq['Title'] + seq['Authors'])
                    .encode('utf-8')).digest()
    ).decode('utf-8')
    seq['_PubID'] = (('PM' + seq['PubMedID'])
                     if seq['PubMedID'] else
                     ('HS' + hashed.rstrip('=')))
    return seq


def attach_references(gene, patient_sequences):
    limit = 300
    filename = os.path.join(
        ROOT, 'internalFiles', 'papersReview',
        '{}PatientSequences.csv'.format(gene)
    )
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    known_sequences = {}
    if os.path.exists(filename):
        with open(filename, 'r') as fp:
            known_sequences = {seq['Accession']: seq
                               for seq in csv.DictReader(fp)}

    header = ['Accession', 'Subtype', 'SamplingYear',
              'PubMedID', 'Title', 'Authors', 'PubYear']
    with open(filename + '.tmp', 'w') as fp:
        writer = csv.DictWriter(fp, header, extrasaction='ignore')
        writer.writeheader()

        print('Fetching references of {} from Genbank...'.format(gene))
        for offset in range(0, len(patient_sequences), limit):
            if offset:
                print('.', end='', file=sys.stderr)
                sys.stderr.flush()
            partial = OrderedDict()
            new_accessions = []
            for seq in patient_sequences[offset:offset + limit]:
                acc = seq['Accession']
                if acc in known_sequences:
                    partial[acc] = dict(known_sequences[acc], **seq)
                else:
                    partial[acc] = seq
                    new_accessions.append(acc)
            if not new_accessions:
                for seq in partial.values():
                    seq = attach_pubid(seq)
                    writer.writerow(seq)
                    yield seq
                continue
            resp = requests.post(
                EFETCH_BASE,
                data={
                    'db': 'nuccore',
                    'id': ','.join(new_accessions),
                    'rettype': 'gb',
                    'retmode': 'xml'
                },
                timeout=None
            )
            root = ET.fromstring(resp.text)
            for gbseq in root.getchildren():
                accession = gbseq.find('GBSeq_primary-accession')
                accession = accession.text
                references = gbseq.findall('GBSeq_references/GBReference')
                pubmed_id = None
                pubyear = None

                ref = ([ref for ref in references
                        if ref.find('GBReference_pubmed')
                        is not None] + [references[0]])[0]

                # find the reference with most information provided
                for ref in references:
                    # reference with a PubMed ID is preferred
                    if ref.find('GBReference_pubmed') is not None:
                        break
                else:
                    for ref in references:
                        # reference with title and authors information
                        # provided is also preferred
                        title = ref.find('GBReference_title')
                        authors = ref.findall(
                            'GBReference_authors/GBAuthor')
                        has_title = (
                            title is not None and
                            (title.text or '').strip())
                        has_authors = \
                            authors and (authors[0].text or '').strip()
                        if has_title and has_authors:
                            break
                    else:
                        # fallback to the first reference
                        ref = references[0]

                title = ref.find('GBReference_title')
                title = title.text if title is not None else ''
                authors = ref.findall('GBReference_authors/GBAuthor')
                if len(authors) > 1:
                    authors = (
                        ', '.join(au.text or '' for au in authors[:-1]) +
                        ' and ' + authors[-1].text)
                elif authors:
                    authors = authors[0].text
                else:
                    authors = ''
                journal = ref.find('GBReference_journal')
                if journal is not None:
                    matches = re.findall(
                        r'\b(?:19|20)\d\d\b', journal.text)
                    if matches:
                        pubyear = matches[0]
                pubmed = ref.find('GBReference_pubmed')
                if pubmed is not None:
                    pubmed_id = pubmed.text
                partial[accession].update({
                    'PubMedID': pubmed_id,
                    'Title': title,
                    'Authors': authors,
                    'PubYear': pubyear
                })
            for seq in partial.values():
                seq = attach_pubid(seq)
                writer.writerow(seq)
                yield seq
    os.rename(filename + '.tmp', filename)
    print('\n{} created'.format(filename))


def find_possible_apobecs(gene, ptseqs):
    filename = os.path.join(
        ROOT, 'data', 'naiveStudies', 'apobec',
        '{}PossibleApobecs.csv'.format(gene)
    )
    apobecs = Counter()
    # g2as = {
    #     'AG': '3G (GG=>AG)',
    #     'AA': '3F (GA=>AA)'
    # }
    for seq in ptseqs:
        naseq = seq['AlignedNASequence']
        # search for all positions has GG=>AG or GG=>AA change
        # deletion gaps should be also considered
        matches = re.finditer('A-*(?=[AG])', naseq)
        muts = {m['Position']: m for m in seq['Mutations']}
        start_codon_apobec_changed = False
        if 1 in muts:
            first = muts[1]
            start_codon_apobec_changed = (
                first['ReferenceText'] == 'M' and
                'I' in first['AminoAcidText'])
        if not start_codon_apobec_changed and not seq['NumStopCodons']:
            # no M=>I and no W=>* changes
            continue

        if seq['NumStopCodons']:
            for mut in seq['Mutations']:
                if '*' in mut['AminoAcidText']:
                    cons = mut['ReferenceText']
                    if cons != 'W':
                        continue
                    aa_pos = mut['Position']
                    cons_prev = get_prevalence(gene, aa_pos, cons)
                    if cons_prev < 98.5:
                        # skip non-conserved position
                        continue
                    apobecs[(aa_pos, 'W=>*')] += 1

        for match in matches:
            start, end = match.span(0)
            # na2 = naseq[end]
            aa_pos = start // 3 + 1
            na_offset = start % 3
            if aa_pos not in muts:
                continue

            mut = muts[aa_pos]
            cons = mut['ReferenceText']

            if mut['IsPartial'] or mut['IsInsertion'] or mut['IsDeletion']:
                continue

            if '*' in mut['AminoAcidText']:
                continue

            cons_prev = get_prevalence(gene, aa_pos, cons)
            if cons_prev < 98.5:
                # skip non-conserved position
                continue

            codon = mut['CodonText']
            for source in get_codons(cons):
                # find G=>A hypermutation
                if source[na_offset] != 'G':
                    continue
                target = source[:na_offset] + 'A' + source[na_offset + 1:]
                if compare_codon(target, codon):
                    target_aa = translate_codon(target)
                    if target_aa == cons:
                        # do not add things like "E=>E"
                        break
                    apobecs[(
                        aa_pos,
                        '{}=>{}'.format(cons, target_aa)
                    )] += 1
                    break

    possible_apobecs = []
    for (pos, mut), count in apobecs.most_common():
        possible_apobecs.append({
            'Position': pos,
            'AAChange': mut,
            'Consensus %': Decimal(
                get_prevalence(gene, pos, mut.split('=>', 1)[0])
            ).quantize(PREC1),
            '# with Stop': count
        })

    for seq in ptseqs:
        muts = {m['Position']: m for m in seq['Mutations']}
        for apobec in possible_apobecs:
            pos = apobec['Position']
            aa = apobec['AAChange'].split('=>', 1)[1]
            if pos in muts and aa in muts[pos]['AminoAcidText']:
                apobec['# Sequence'] = apobec.get('# Sequence', 0) + 1
    for apobec in possible_apobecs:
        apobec['% with Stop'] = Decimal(
            apobec['# with Stop'] * 100 / apobec['# Sequence']
        ).quantize(PREC0)
    possible_apobecs = sorted(
        possible_apobecs,
        key=lambda a: (a['Position'], a['AAChange']))
    possible_apobecs = [a for a in possible_apobecs
                        if a['% with Stop'] > 50]

    csv_writer(
        filename, possible_apobecs,
        ['Position', 'AAChange', 'Consensus %',
         '% with Stop', '# Sequence'],
        writer_options={
            'extrasaction': 'ignore'
        }
    )


if __name__ == '__main__':
    for gene in ('gag', 'gp41'):
        ptseqs = get_patient_sequences(gene)
        find_possible_apobecs(gene, ptseqs)
        ptseqs = attach_apobecs(gene, ptseqs)
        create_review_table(gene, ptseqs)
        export_naive_sequences(gene, ptseqs)

        csv_writer(
            os.path.join(
                ROOT, 'data', 'naiveStudies', 'apobec',
                '{}NaiveADIndex.csv'.format(gene)),
            aggregate_naiveseqs_adindex(gene, ptseqs),
            ['Accession', 'QCIssue', 'NumAPOBECs',
             'NumConservedAPOBECSites', 'ADIndex'])
