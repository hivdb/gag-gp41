#! /usr/bin/env python

import os
import re
import csv
import base64
import hashlib
import requests
import xlsxwriter as xw
from itertools import groupby
from data_reader import data_reader, ROOT

SEARCH_INTERFACE = ('https://www.hiv.lanl.gov/components/'
                    'sequence/HIV/search/search.html')
SEARCH_TARGET = ('https://www.hiv.lanl.gov/components/'
                 'sequence/HIV/search/search.comp')
GENE_RANGE = {
    'gag': (790, 2289),
    'gp41': (7758, 8792)
}


def uniq_accessions_filename(gene):
    return os.path.join(
        ROOT, 'result_data',
        '{}AccessionPerPatient.txt'.format(gene)
    )


def get_accessions(gene):
    seqs = data_reader(
        os.path.join(ROOT, 'data', 'genbankSequences',
                     'Comp.{}.txt'.format(gene.lower())),
        delimiter='\t')
    return [seq['Accession'] for seq in seqs]


def create_review_table(gene):
    seqs = data_reader(
        os.path.join(ROOT, 'data', 'genbankSequences',
                     'Comp.{}.txt'.format(gene.lower())),
        delimiter='\t')
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

    fact_table = fact_table_merged
    uniq_accs = set(filter_lanl(gene))
    uniq_seqs = []
    for seq in seqs:
        hashed = base64.urlsafe_b64encode(
            hashlib.md5((seq['Title'] + seq['Authors'])
                        .encode('utf-8')).digest()
        ).decode('utf-8')
        if seq['Accession'] in uniq_accs:
            seq['_PubID'] = (('PM' + seq['PubMedID'])
                             if seq['PubMedID'] else
                             ('HS' + hashed.rstrip('=')))
            uniq_seqs.append(seq)
    uniq_seqs = sorted(uniq_seqs, key=lambda s: s['_PubID'])
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
        ROOT, 'result_data',
        '{}ReviewTable.csv'.format(gene)
    ), 'w') as fp:
        # fp.write('\ufeff')  # BOM for Excel
        writer = csv.DictWriter(fp, [
            'PubID', 'PubMedID', 'PubYear', 'NumPts', 'NumIsolates',
            'NumLANLIsolates', 'Title', 'Authors', 'RxStatus', 'Notes'])
        writer.writeheader()
        writer.writerows(results)
    export_excel_table(
        os.path.join(
            ROOT, 'result_data',
            '{}ReviewTable.xlsx'.format(gene)
        ),
        results)


def export_excel_table(filename, rows):
    workbook = xw.Workbook(filename)
    worksheet = workbook.add_worksheet('main')
    worksheet_rx_status = workbook.add_worksheet('validRxStatus')
    fontsize = 14
    workbook.formats[0].set_font_size(fontsize)
    headers = [
        'PubID', 'PubMedID', 'PubYear', 'NumPts', 'NumIsolates',
        'NumLANLIsolates', 'Title', 'Authors', 'RxStatus', 'Notes'
    ]
    valid_rx_status = sorted([
        'Lab', 'Rx', 'Rx-PI', 'Rx=>Naive', 'Check', 'Naive', 'Unknown',
        'Mixed', 'ProbablyNaive', 'Unpublished', 'NonM', 'Conflict'
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
            return [r.strip() for r in fp.readlines()]
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
                row['PAT id(SSAM)'] or row['Accession'], row['Accession'])
    uniq_accessions = list(uniq_accessions.values())
    with open(result_filename, 'w') as fp:
        fp.write('\n'.join(sorted(uniq_accessions)))
    return uniq_accessions


if __name__ == '__main__':
    create_review_table('gag')
    create_review_table('gp41')
