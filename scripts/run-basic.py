import os
from itertools import chain
from data_reader import ROOT
from data_writer import csv_writer
from analysis_functions import (
    codon_changes_per_person,
    aggregate_aa_changes_by_pos)

MAX_UNUSUALS = {
    'gag': 15,
    'gp41': 10
}


def main():
    def cc_keyfunc(c):
        return int(c['PID']), c['Rx'], c['Pos']

    for gene in ('gag', 'gp41'):

        csv_writer(
            os.path.join(ROOT, 'resultData', 'aaChangesByPosWPrev',
                         '{}.csv'.format(gene)),
            chain(
                aggregate_aa_changes_by_pos(gene, 'PIs', 'PIs'),
                aggregate_aa_changes_by_pos(gene, 'NNRTIs', 'NNRTIs'),
            ),
            ['Group', 'Pos', 'PreAA', 'PostAA',
             'NumPts', 'PrePrev', 'PostPrev', 'Fold', 'LogFold'])

        csv_writer(
            os.path.join(ROOT, 'resultData', 'codonChangesByPt',
                         '{}.csv'.format(gene)),
            sorted(
                codon_changes_per_person(gene, ('PIs', 'NNRTIs')),
                key=cc_keyfunc),
            ['PID', 'Rx', 'Pos', 'Type', 'Codons', 'NumNAChanges', 'AAs'])


if __name__ == '__main__':
    main()
