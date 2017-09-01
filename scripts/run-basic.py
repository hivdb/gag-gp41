import os
import sys
import csv
from itertools import chain
from data_reader import ROOT
from analysis_functions import codon_changes_per_person
from analysis_functions import aggregate_aa_changes_by_pos


def save_csv(path, data, headers=None):
    with open(path, 'w') as fp:
        writer = csv.DictWriter(fp, headers)
        writer.writeheader()
        writer.writerows(data)
    print('{} created'.format(path))


def main():
    save_csv(
        os.path.join(ROOT, 'result_data', 'GagAAChangesByPosWPrev.csv'),
        chain(
            aggregate_aa_changes_by_pos('Gag', 'PIs', 'PIs'),
            aggregate_aa_changes_by_pos('Gag', 'NNRTIs', 'NNRTIs'),
        ),
        ['Group', 'Pos', 'PreAA', 'PostAA',
         'NumPts', 'PrePrev', 'PostPrev', 'Fold', 'LogFold'])

    save_csv(
        os.path.join(ROOT, 'result_data', 'Gp41AAChangesByPosWPrev.csv'),
        chain(
            aggregate_aa_changes_by_pos('gp41', 'PIs', 'PIs'),
            aggregate_aa_changes_by_pos('gp41', 'NNRTIs', 'NNRTIs'),
        ),
        ['Group', 'Pos', 'PreAA', 'PostAA',
         'NumPts', 'PrePrev', 'PostPrev', 'Fold', 'LogFold'])

    def cc_keyfunc(c):
        return int(c['PID']), c['Rx'], c['Pos']

    save_csv(
        os.path.join(ROOT, 'result_data', 'GagCodonChangesByPt.csv'),
        sorted(
            codon_changes_per_person('Gag', ('PIs', 'NNRTIs')),
            key=cc_keyfunc),
        ['PID', 'Rx', 'Pos', 'Type', 'Codons', 'NumNAChanges', 'AAs'])

    save_csv(
        os.path.join(ROOT, 'result_data', 'Gp41CodonChangesByPt.csv'),
        sorted(
            codon_changes_per_person('gp41', ('PIs', 'NNRTIs')),
            key=cc_keyfunc),
        ['PID', 'Rx', 'Pos', 'Type', 'Codons', 'NumNAChanges', 'AAs'])

    save_csv(
        os.path.join(ROOT, 'result_data', 'PRCodonChangesByPt.csv'),
        sorted(
            codon_changes_per_person('PR', ('PIs', 'NNRTIs')),
            key=cc_keyfunc),
        ['PID', 'Rx', 'Pos', 'Type', 'Codons', 'NumNAChanges', 'AAs'])

    save_csv(
        os.path.join(ROOT, 'result_data', 'RTCodonChangesByPt.csv'),
        sorted(
            codon_changes_per_person('RT', ('PIs', 'NNRTIs')),
            key=cc_keyfunc),
        ['PID', 'Rx', 'Pos', 'Type', 'Codons', 'NumNAChanges', 'AAs'])


if __name__ == '__main__':
    main()