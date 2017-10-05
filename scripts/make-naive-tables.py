#! /usr/bin/env python

import os
import re
import csv
import sys
import json
import base64
import hashlib
import requests
import xlsxwriter as xw
from itertools import groupby
import xml.etree.ElementTree as ET
from collections import OrderedDict
from subprocess import Popen, PIPE

from data_reader import data_reader, ROOT, CONSENSUS as _CONS

SEARCH_INTERFACE = ('https://www.hiv.lanl.gov/components/'
                    'sequence/HIV/search/search.html')
SEARCH_TARGET = ('https://www.hiv.lanl.gov/components/'
                 'sequence/HIV/search/search.comp')
EFETCH_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
NAIVES = ('Naive', 'PINaive', 'ProbablyNaive')
GENE_RANGE = {
    'gag': (790, 2289),
    'gp41': (7758, 8792)
}

CONSENSUS = {
    'gag': _CONS['Gag'],
    'gp41': _CONS['gp41']
}

REVIEW_TABLE_HEADERS = [
    'PubID', 'PubMedID', 'PubYear', 'NumPts', 'NumIsolates',
    'NumLANLIsolates', 'Title', 'Authors', 'Subtypes',
    'RxStatus', 'Notes'
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


def nucamino(gene, sequences):
    cmd = os.path.join(ROOT, 'scripts', 'nucamino')
    fasta = '\n'.join(
        '>{Accession}\n{NASequence}'.format(**seq)
        for seq in sequences
    )
    proc = Popen(
        [cmd, 'hiv1b', '-g', gene.upper(), '--output-format', 'json'],
        stdin=PIPE, stdout=PIPE, stderr=sys.stdout)
    return proc.communicate(fasta.encode('UTF-8'))[0]


def attach_alignments(gene, sequences):
    genesize = int(CONSENSUS[gene]['Size'])

    # load cached aligned results first
    filename = os.path.join(
        ROOT, 'local',
        '{}NaiveAlignedSequences.json'.format(gene))
    cached_aligned_results = []
    if os.path.exists(filename):
        print('Local file {} found'.format(filename))
        with open(filename) as fp:
            try:
                cached_aligned_results = json.load(fp)[gene.upper()]
            except (KeyError, ValueError):
                pass
    known_seqs = {r['Name'] for r in cached_aligned_results}

    # align sequences if not aligned before
    input_sequences = [
        s for s in sequences if s['Accession'] not in known_seqs]
    if input_sequences:
        stdout = nucamino(gene, input_sequences)
        aligned_results = (
            cached_aligned_results +
            json.loads(stdout.decode('UTF-8'))[gene.upper()]
        )
        with open(filename, 'w') as fp:
            json.dump({gene.upper(): aligned_results}, fp)
            print('{} created'.format(filename))
    else:
        aligned_results = cached_aligned_results

    sequences = OrderedDict([(seq['Accession'], seq) for seq in sequences])
    for result in cached_aligned_results + aligned_results:
        acc = result['Name']
        if acc not in sequences:
            continue
        sequence = sequences[acc]
        report = result['Report']
        control = report['ControlLine']
        naseq = report['NucleicAcidsLine']
        firstaa = report['FirstAA']
        lastaa = report['LastAA']
        firstna = report['FirstNA']
        lastna = report['LastNA']
        mutations = report['Mutations']
        frameshifts = report['FrameShifts']
        # trim not aligned seq
        trimed_naseq = sequence['NASequence'][firstna - 1:lastna].upper()
        # remove insertions
        aligned_naseq = ''.join(
            na for na, c in zip(naseq, control) if c != '+')
        aligned_naseq = aligned_naseq.replace(' ', '-')
        if firstaa > 1:
            aligned_naseq = '...' * firstaa + aligned_naseq
        if lastaa < genesize:
            aligned_naseq += '...' * (genesize - lastaa)
        sequence.update({
            'TrimedNASequence': trimed_naseq,
            'AlignedNASequence': aligned_naseq,
            'FirstAA': firstaa,
            'LastAA': lastaa,
            'NumInsertions': len([m for m in mutations if m['IsInsertion']]),
            'NumDeletions': len([m for m in mutations if m['IsDeletion']]),
            'NumStopCodons': len([m for m in mutations
                                  if '*' in m['AminoAcidText']]),
            'NumFrameShifts': len(frameshifts),
            'Mutations': mutations
        })
    return list(sequences.values())


def export_naive_sequences(gene, ptseqs):
    filename = os.path.join(
        ROOT, 'result_data',
        '{}NaiveSequences.csv'.format(gene.lower())
    )
    aligned_fasta = os.path.join(
        ROOT, 'result_data',
        '{}AlignedNaiveSequences.fasta'.format(gene.lower())
    )
    unaligned_fasta = os.path.join(
        ROOT, 'result_data',
        '{}UnalignedNaiveSequences.fasta'.format(gene.lower())
    )

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
    # remove non-naive sequences
    ptseqs = [seq for seq in ptseqs if seq['_PubID'] in pubids]
    ptseqs = attach_alignments(gene, ptseqs)
    genesize = int(CONSENSUS[gene]['Size'])
    siteheaders = ['P{}'.format(i) for i in range(1, genesize + 1)]

    with open(filename, 'w') as fp, \
            open(aligned_fasta, 'w') as af, \
            open(unaligned_fasta, 'w') as uf:
        writer = csv.DictWriter(
            fp, ['PMID', 'Accession', 'RxStatus', 'lanlSubtype'] + siteheaders)
        writer.writeheader()
        ptseqs = sorted(ptseqs, key=lambda s: s['Accession'])
        af.write('\n'.join(
            '>{Accession}|{Subtype}\n{AlignedNASequence}'.format(**s)
            for s in ptseqs
        ))
        print('{} created'.format(aligned_fasta))
        uf.write('\n'.join(
            '>{Accession}|{Subtype}\n{TrimedNASequence}'.format(**s)
            for s in ptseqs
        ))
        print('{} created'.format(unaligned_fasta))
        for seq in ptseqs:
            firstaa = seq['FirstAA']
            lastaa = seq['LastAA']
            row = {
                'PMID': pubids[seq['_PubID']],
                'Accession': seq['Accession'],
                'RxStatus': 'Naive',
                'lanlSubtype': seq['Subtype']
            }
            muts = {m['Position']: m for m in seq['Mutations']}
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
            writer.writerow(row)
        print('{} created'.format(filename))


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
            if not result.get('RxStatus'):
                result['RxStatus'] = fact.get('RxStatus')

    results = sorted(
        results.values(), key=lambda r: (-r['NumLANLIsolates'], r['PubID']))
    with open(os.path.join(
        ROOT, 'internalFiles', 'papersReview',
        '{}ReviewTable.csv'.format(gene)
    ), 'w') as fp:
        # fp.write('\ufeff')  # BOM for Excel
        writer = csv.DictWriter(fp, REVIEW_TABLE_HEADERS)
        writer.writeheader()
        writer.writerows(results)
        print('{} created'.format(fp.name))
    export_excel_table(
        os.path.join(
            ROOT, 'internalFiles', 'papersReview',
            '{}ReviewTable.xlsx'.format(gene)
        ),
        results)
    export_naive_papers_table(
        os.path.join(
            ROOT, 'result_data',
            '{}NaivePapers.csv'.format(gene)
        ),
        results)


def export_naive_papers_table(filename, rows):
    rows = naive_papers(rows)
    with open(filename, 'w') as fp:
        writer = csv.DictWriter(fp, CLEAN_TABLE_HEADERS)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                'PubMedID': row['PubMedID'],
                'PubYear': row['PubYear'],
                'Title': row['Title'],
                'Authors': row['Authors'],
                'NumPatients': row['NumLANLIsolates']
            })
        print('{} created'.format(fp.name))


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
        'Entry-Naive', 'Conflict'
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
    worksheet.set_column('F:F', width=15, cell_format=numfmt)  # #LANLIsolates
    worksheet.set_column('G:G', width=60, cell_format=wrapfmt)  # Title
    worksheet.set_column('H:H', width=60, cell_format=wrapfmt)  # Authors
    worksheet.set_column('I:I', width=40, cell_format=wrapfmt)  # Subtypes
    worksheet.set_column('J:J', width=15, cell_format=fmt)  # RxStatus
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


def download_from_lanl(gene):
    filename = os.path.join(
        ROOT, 'local', 'comp-{}-lanl.txt'.format(gene)
    )
    if os.path.exists(filename):
        print('Local file {} found'.format(filename))
        return filename

    session = requests.Session()
    session.get(SEARCH_INTERFACE)  # fetch cookies first

    start, stop = GENE_RANGE[gene.lower()]
    resp = session.post(
        SEARCH_TARGET,
        data={
            'master': 'HIV-1',  # only HIV-1 seqs
            'value SequenceMap SM_start 2': start,
            'value SequenceMap SM_stop 2': stop,
            'max_rec': 100,
            'value SEQ_SAMple SSAM_badseq 4': 'on',
            'action': 'search'
        },
        timeout=None
    )
    result = resp.text
    idx, = re.findall(r'<input .*\bname="id" value="([^"]+)"', result)
    resp = session.post(
        SEARCH_TARGET,
        data={
            'incl_seq': 'on',
            'save_tbl': 'OK',
            'id': idx
        },
        timeout=None
    )
    with open(filename, 'wb') as fp:
        fp.write(resp.content)
        print('{} downloaded'.format(filename))
    return filename


def get_patient_sequences(gene):
    filename = download_from_lanl(gene)
    ptseqs = {}
    with open(filename, 'r') as fp:
        next(fp)  # skip the first line
        sequences = csv.DictReader(fp, delimiter='\t')
        for seq in sequences:
            patid = seq['PAT id(SSAM)']
            accession = seq['Accession']
            subtype = seq['Subtype']
            sampyear = seq['Sampling Year']
            naseq = seq['Sequence']
            if subtype in ('O', 'N', 'P'):
                continue
            if not patid:
                continue
            if patid not in ptseqs or \
                    accession < ptseqs[patid]['Accession']:
                # always use the first accession to be consistent
                ptseqs[patid] = {
                    'Accession': accession,
                    'Subtype': subtype,
                    'SamplingYear': sampyear,
                    'NASequence': naseq
                }
    ptseqs = sorted(ptseqs.values(), key=lambda seq: seq['Accession'])
    return list(attach_references(gene, ptseqs))


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
                    partial[acc] = known_sequences[acc]
                    partial[acc]['NASequence'] = seq['NASequence']
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


if __name__ == '__main__':
    for gene in ('gag', 'gp41'):
        ptseqs = get_patient_sequences(gene)
        create_review_table(gene, ptseqs)
        export_naive_sequences(gene, ptseqs)
