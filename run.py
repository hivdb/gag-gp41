import os
import csv
from itertools import chain
from analysis_functions import codon_changes_per_person
from analysis_functions import aggregate_aa_changes_by_pos

PWD = os.path.dirname(__file__)


def save_csv(path, data, headers=None):
    with open(path, 'w') as fp:
        writer = csv.DictWriter(fp, headers, delimiter='\t')
        writer.writeheader()
        writer.writerows(data)


def main():
    save_csv(
        os.path.join(PWD, 'result_data', 'GagAAChangesByPosWPrev.txt'),
        chain(
            aggregate_aa_changes_by_pos('Gag', 'PIs', 'PIs'),
            aggregate_aa_changes_by_pos(
                'Gag', ('None', 'NNRTIs'), 'Controls'),
        ),
        ['Group', 'Pos', 'PreAA', 'PostAA',
         'NumPts', 'PrePrev', 'PostPrev', 'Fold', 'LogFold'])

    save_csv(
        os.path.join(PWD, 'result_data', 'Gp41AAChangesByPosWPrev.txt'),
        chain(
            aggregate_aa_changes_by_pos('gp41', 'PIs', 'PIs'),
            aggregate_aa_changes_by_pos(
                'gp41', ('None', 'NNRTIs'), 'Controls'),
        ),
        ['Group', 'Pos', 'PreAA', 'PostAA',
         'NumPts', 'PrePrev', 'PostPrev', 'Fold', 'LogFold'])

    def cc_keyfunc(c):
        return int(c['PID']), c['Rx'], c['Pos']

    save_csv(
        os.path.join(PWD, 'result_data', 'GagCodonChangesByNSS_And_Pt.txt'),
        sorted(codon_changes_per_person('Gag'), key=cc_keyfunc),
        ['PID', 'Rx', 'Pos', 'Type', 'Codons', 'NumNAChanges', 'AAs'])

    save_csv(
        os.path.join(PWD, 'result_data', 'Gp41CodonChangesByNSS_And_Pt.txt'),
        sorted(codon_changes_per_person('gp41'), key=cc_keyfunc),
        ['PID', 'Rx', 'Pos', 'Type', 'Codons', 'NumNAChanges', 'AAs'])

    save_csv(
        os.path.join(PWD, 'result_data', 'PRCodonChangesByNSS_And_Pt.txt'),
        sorted(codon_changes_per_person('PR'), key=cc_keyfunc),
        ['PID', 'Rx', 'Pos', 'Type', 'Codons', 'NumNAChanges', 'AAs'])

    save_csv(
        os.path.join(PWD, 'result_data', 'RTCodonChangesByNSS_And_Pt.txt'),
        sorted(codon_changes_per_person('RT'), key=cc_keyfunc),
        ['PID', 'Rx', 'Pos', 'Type', 'Codons', 'NumNAChanges', 'AAs'])


if __name__ == '__main__':
    main()
