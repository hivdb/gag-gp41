import os
from itertools import chain
from data_reader import ROOT
from data_writer import csv_writer
from analysis_functions import (
    get_most_common_subtypes,
    codon_changes_per_person,
    aggregate_aa_changes_by_pos,
    aggregate_naiveseqs_stat,
    aggregate_naiveseqs_posstat,
    aggregate_mut_prevalence)


def main():
    def cc_keyfunc(c):
        return int(c['PID']), c['Rx'], c['Pos']

    for gene in ('gag', 'gp41'):

        major_subtypes = [None] + get_most_common_subtypes(gene)
        header = None
        all_prevalence = []
        for subtype in major_subtypes:
            header, prevs = aggregate_mut_prevalence(gene, subtype)
            all_prevalence.append(prevs)
        csv_writer(
            os.path.join(ROOT, 'data', 'naiveStudies',
                         '{}AAPrevalence.csv'.format(gene)),
            chain(*all_prevalence),
            header)

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

        csv_writer(
            os.path.join(ROOT, 'data', 'naiveStudies',
                         '{}StatBySeq.csv'.format(gene)),
            aggregate_naiveseqs_stat(gene),
            ['Accession', 'PMID', 'Gene', 'Subtype',
             'NumAAChanges', 'NumInsertions', 'NumDeletions',
             'NumStopCodons', 'NumApobecs', 'NumUnusuals', 'NumFrameShifts'])

        csv_writer(
            os.path.join(ROOT, 'data', 'naiveStudies',
                         '{}StatByPos.csv'.format(gene)),
            aggregate_naiveseqs_posstat(gene),
            ['Gene', 'AAPosition', 'NumAAChanges',
             'NumStopCodons', 'NumInsertions', 'NumDeletions'])


if __name__ == '__main__':
    main()
