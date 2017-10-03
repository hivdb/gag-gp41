import os
import csv
from itertools import chain
from data_reader import ROOT
from analysis_functions import (
    codon_changes_per_person,
    aggregate_aa_changes_by_pos,
    aggregate_naiveseqs_distribution,
    aggregate_mut_prevalence)


def save_csv(path, data, headers=None):
    with open(path, 'w') as fp:
        writer = csv.DictWriter(fp, headers)
        writer.writeheader()
        writer.writerows(data)
    print('{} created'.format(path))


def main():
    save_csv(
        os.path.join(ROOT, 'result_data', 'MutPrevalence.csv'),
        chain(
            aggregate_mut_prevalence('Gag'),
            aggregate_mut_prevalence('gp41'),
        ),
        ['Gene', 'Pos', 'AA', 'Pcnt', 'Count', 'PosTotal'])

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

    for gene in ('Gag', 'gp41'):
        diffs, stops = aggregate_naiveseqs_distribution(gene)
        save_csv(
            os.path.join(
                ROOT, 'result_data',
                '{}NaiveSeqsChangesDistribution.csv'.format(gene)
            ),
            diffs,
            ['NumChanges', 'NumSequences'])
        save_csv(
            os.path.join(
                ROOT, 'result_data',
                '{}NaiveSeqsStopCodonsDistribution.csv'.format(gene)
            ),
            stops,
            ['NumStopCodons', 'NumSequences'])


if __name__ == '__main__':
    main()
