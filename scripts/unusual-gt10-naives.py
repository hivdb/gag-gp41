#! /usr/bin/env python

import os
import csv
from data_reader import ROOT, CONSENSUS, naive_sequence_reader, get_prevalence

for gene, gene_size in (('gag', 500), ('gp41', 345)):
    result_rows = []
    stats = naive_sequence_reader(gene)
    stats = [s for s in stats if int(s['NumUnusuals']) > 10]
    for seq in stats:
        for pos in range(1, gene_size + 1):
            aa = seq['P{}'.format(pos)]
            cons = CONSENSUS[gene]['AASeq'][pos - 1]
            if len(aa) > 1:
                aa = aa.replace(cons, '')
            if aa == '-':
                aa = cons
            prev = get_prevalence(gene, pos, aa)
            if len(aa) > 1:
                prev = 'NA'
            result_rows.append({
                'Accession': seq['Accession'],
                'Subtype': seq['lanlSubtype'],
                'NumUnusuals': seq['NumUnusuals'],
                'Position': pos,
                'AA': aa,
                'Prevalence': prev
            })
    with open(os.path.join(
            ROOT, 'local', '{}UnusualGt10NaiveAllSites.csv'.format(gene)
    ), 'w') as fp:
        writer = csv.DictWriter(fp, ['Accession', 'Subtype', 'NumUnusuals',
                                     'Position', 'AA', 'Prevalence'])
        writer.writeheader()
        writer.writerows(result_rows)
