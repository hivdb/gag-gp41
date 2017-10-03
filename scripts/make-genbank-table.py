#! /usr/bin/env python

import os
import re
import csv
import base64
import hashlib
import requests
import xlsxwriter as xw
from itertools import groupby
from codonutils import translate_codon
from data_reader import data_reader, ROOT, CONSENSUS as _CONS

SEARCH_INTERFACE = ('https://www.hiv.lanl.gov/components/'
                    'sequence/HIV/search/search.html')
SEARCH_TARGET = ('https://www.hiv.lanl.gov/components/'
                 'sequence/HIV/search/search.comp')
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
    'NumLANLIsolates', 'Title', 'Authors', 'RxStatus', 'Notes'
]

CLEAN_TABLE_HEADERS = [
    'PubMedID', 'PubYear', 'Title', 'Authors', 'NumPatients'
]


def uniq_accessions_filename(gene):
    return os.path.join(
        ROOT, 'internalFiles', 'papersReview',
        '{}AccessionPerPatient.txt'.format(gene)
    )


def genbank_sequences_reader(gene):
    return data_reader(
        os.path.join(ROOT, 'internalFiles', 'genbankSequences',
                     'Comp.{}.txt'.format(gene.lower())),
        delimiter='\t')


def get_accessions(gene):
    seqs = genbank_sequences_reader(gene)
    return [seq['Accession'] for seq in seqs]


def get_fact_table(gene):
    fact_table = data_reader(
        os.path.join(ROOT, 'data', '{}GenbankFact.csv'
                     .format(gene.lower())))
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


def get_sequences_per_patients(gene):
    seqs = genbank_sequences_reader(gene)
    uniq_accs = {acc: subtype for acc, subtype in filter_lanl(gene)}
    uniq_seqs = []
    for seq in seqs:
        hashed = base64.urlsafe_b64encode(
            hashlib.md5((seq['Title'] + seq['Authors'])
                        .encode('utf-8')).digest()
        ).decode('utf-8')
        acc = seq['Accession']
        if acc in uniq_accs:
            seq['Subtype'] = uniq_accs[acc]
            seq['_PubID'] = (('PM' + seq['PubMedID'])
                             if seq['PubMedID'] else
                             ('HS' + hashed.rstrip('=')))
            uniq_seqs.append(seq)
    return sorted(uniq_seqs, key=lambda s: s['_PubID'])


def iter_aas(consensus, naseq):
    for i, cons in enumerate(consensus):
        codon = naseq[i * 3:i * 3 + 3]
        if codon == '---':
            aas = 'd'
        elif '-' in codon:
            aas = 'X'
        else:
            aas = translate_codon(codon)
            if aas == cons:
                aas = '-'
            elif len(aas) > 4:
                aas = 'X'
        yield aas


def create_naive_sequences_table(gene):
    filename = os.path.join(
        ROOT, 'result_data',
        '{}NaiveSequences.csv'.format(gene.lower())
    )
    fact_table = get_fact_table(gene)
    uniq_seqs = get_sequences_per_patients(gene)
    pubids = {pubid for pubid, f in fact_table.items()
              if f['RxStatus'] in NAIVES}
    genesize = int(CONSENSUS[gene]['Size'])
    siteheaders = ['P{}'.format(i) for i in range(1, genesize + 1)]
    with open(filename, 'w') as fp:
        writer = csv.DictWriter(
            fp, ['PMID', 'Accession', 'RxStatus', 'lanlSubtype'] + siteheaders)
        writer.writeheader()
        for seq in sorted(uniq_seqs, key=lambda s: s['Accession']):
            if seq['_PubID'] not in pubids:
                continue
            row = {
                'PMID': seq['PubMedID'],
                'Accession': seq['Accession'],
                'RxStatus': 'Naive',
                'lanlSubtype': seq['Subtype']
            }
            naseq = seq['NASeq']
            aaseq = iter_aas(CONSENSUS[gene]['AASeq'], naseq)
            for i, aas in enumerate(aaseq):
                pos = i + 1
                row['P{}'.format(pos)] = aas
            writer.writerow(row)


def create_review_table(gene):
    fact_table = get_fact_table(gene)
    uniq_seqs = get_sequences_per_patients(gene)
    grouped = groupby(uniq_seqs, lambda s: s['_PubID'])
    results = {}
    for pubid, group_seqs in grouped:
        group_seqs = list(group_seqs)
        seq = group_seqs[0]
        fact = fact_table.get(pubid, {})
        if fact.get('PubIDCorrection'):
            pubid = fact['PubIDCorrection']

        if pubid not in results:
            results[pubid] = {
                'PubID': pubid,
                'PubMedID': fact.get('PMID', seq['PubMedID']),
                'PubYear': fact.get('PubYr'),
                'NumPts': fact.get('NumPts'),
                'NumIsolates': fact.get('NumIsolates'),
                'NumLANLIsolates': len(group_seqs),
                'Title': seq['Title'],
                'Authors': seq['Authors'],
                'RxStatus': fact.get('RxStatus'),
                'Notes': fact.get('Notes'),
            }
        else:
            result = results[pubid]
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

    results = sorted(results.values(), key=lambda r: -r['NumLANLIsolates'])
    with open(os.path.join(
        ROOT, 'internalFiles', 'papersReview',
        '{}ReviewTable.csv'.format(gene)
    ), 'w') as fp:
        # fp.write('\ufeff')  # BOM for Excel
        writer = csv.DictWriter(fp, REVIEW_TABLE_HEADERS)
        writer.writeheader()
        writer.writerows(results)
    export_excel_table(
        os.path.join(
            ROOT, 'internalFiles', 'papersReview',
            '{}ReviewTable.xlsx'.format(gene)
        ),
        results)
    export_clean_table(
        os.path.join(
            ROOT, 'result_data',
            '{}ReferencesTable.csv'.format(gene)
        ),
        results)


def export_clean_table(filename, rows):
    rows = [row for row in rows if row['RxStatus'] in NAIVES]
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


def export_excel_table(filename, rows):
    workbook = xw.Workbook(filename)
    worksheet = workbook.add_worksheet('main')
    worksheet_rx_status = workbook.add_worksheet('validRxStatus')
    fontsize = 14
    workbook.formats[0].set_font_size(fontsize)
    headers = REVIEW_TABLE_HEADERS
    valid_rx_status = sorted([
        'Lab', 'Rx', 'Rx-PI', 'Rx=>Naive', 'Check', 'Naive', 'Unknown',
        'Mixed', 'PINaive', 'ProbablyNaive', 'Unpublished', 'NonM', 'Conflict'
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
    worksheet.set_column('I:I', width=15, cell_format=fmt)  # RxStatus
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

    workbook.set_size(2400, 1600)
    workbook.close()


def filter_lanl(gene):
    result_filename = uniq_accessions_filename(gene)
    if os.path.exists(result_filename):
        with open(result_filename) as fp:
            return list(csv.reader(fp))
    session = requests.Session()
    # fetch cookies first
    session.get(SEARCH_INTERFACE)

    accessions = get_accessions(gene)

    limit = 4000
    uniq_accessions = {}
    for offset in range(0, len(accessions), limit):
        # only HIV-1 seqs
        start, stop = GENE_RANGE[gene.lower()]
        resp = session.post(
            SEARCH_TARGET,
            data={
                'master': 'HIV-1',
                'value SequenceMap SM_start 2': start,
                'value SequenceMap SM_stop 2': stop,
                'max_rec': 100,
                'action': 'search'
            },
            files={
                'user_acc_file': '\n'.join(accessions[offset:offset+limit])
            }
        )
        result = resp.text
        idx, = re.findall(r'<input .*\bname="id" value="([^"]+)"', result)
        resp = session.post(
            SEARCH_TARGET,
            data={
                'save_tbl': 'OK',
                'id': idx
            })
        out = resp.text.strip().splitlines()
        out.pop(0)
        reader = csv.DictReader(out, delimiter='\t')
        for row in reader:
            uniq_accessions.setdefault(
                row['PAT id(SSAM)'] or row['Accession'],
                (row['Accession'], row['Subtype']))
    uniq_accessions = list(uniq_accessions.values())
    with open(result_filename, 'w') as fp:
        writer = csv.writer(fp)
        uniq_accessions = sorted(uniq_accessions)
        writer.writerows(uniq_accessions)
    return uniq_accessions


if __name__ == '__main__':
    create_review_table('gag')
    create_review_table('gp41')
    create_naive_sequences_table('gag')
    create_naive_sequences_table('gp41')
