#! /usr/bin/env python

import os
from multiprocessing import Pool, cpu_count

from data_writer import csv_writer
from data_reader import ROOT, any_fasta_reader

ACGT = set('ACGT')


def get_distance(args):
    acc1, seq1, acc2, seq2 = args
    diff = 0
    total = len(seq1)
    for na1, na2 in zip(seq1, seq2):
        if na1 in ACGT and na2 in ACGT:
            diff += na1 != na2
    return {
        'Sequence1': acc1,
        'Sequence2': acc2,
        'Distance': diff / total
    }


def make_sequence_pairs(sequences):
    sequences = list(sequences)
    for acc1, seq1 in sequences:
        for acc2, seq2 in sequences:
            if acc1 >= acc2:
                continue
            yield acc1, seq1, acc2, seq2


def calc_distances(gene, pool):
    print('Calculating distances for {} naive sequences...'.format(gene))
    sequences = any_fasta_reader(os.path.join(
        ROOT, 'data', 'naiveStudies', '{}NaiveAligned.fas'.format(gene)
    ))
    sequences = list(sequences)
    sequence_pairs = make_sequence_pairs(sequences)
    return pool.imap(get_distance, sequence_pairs, chunksize=1000)


def main():
    with Pool(max(1, cpu_count() - 2)) as pool:
        for gene in ('gag', 'gp41'):
            result = calc_distances(gene, pool)
            csv_writer(os.path.join(
                ROOT, 'local', 'naiveStudies',
                '{}NaiveDistance.csv'.format(gene)
            ), result, ['Sequence1', 'Sequence2', 'Distance'])


if __name__ == '__main__':
    main()
