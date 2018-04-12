#! /usr/bin/env python
import csv

from data_reader import sequence_reader, any_fasta_reader

FASTA_DIR = '/app/data/fasta'
INTERNAL_FASTA_DIR = '/app/internalFiles/fasta'
TARGET = '/app/internalFiles/sequences.csv'
HEADER = ['PID', 'TimePoint', 'Accession', 'Date',
          'Gene', 'FirstAA', 'LastAA', 'NASequence']


sequences = list(sequence_reader())

new_sequences = {}
for gene in ('gag', 'gp41'):
    filename = '{}/{}Aligned.fas'.format(FASTA_DIR, gene)
    pis_only = open(
        '{}/{}PIs.aln.fasta.txt'.format(INTERNAL_FASTA_DIR, gene), 'w')
    nnrtis_only = open(
        '{}/{}NNRTIs.aln.fasta.txt'.format(INTERNAL_FASTA_DIR, gene), 'w')
    for header, seq in any_fasta_reader(filename):
        ptid, dc, ts = header.split('|', 1)[1].split('_', 2)
        key = (gene, ptid, ts)
        new_sequences[key] = seq
        if dc == 'PIs':
            pis_only.write('>{}\n{}\n'.format(header, seq))
        elif dc == 'NNRTIs':
            nnrtis_only.write('>{}\n{}\n'.format(header, seq))

with open(TARGET, 'w') as fp:
    writer = csv.DictWriter(fp, HEADER)
    writer.writeheader()
    for seq in sequences:
        seqdata = seq._data
        key = (seq.gene, seq.pid, seq.time_point)
        first_na0 = seq.first_aa * 3 - 3
        last_na0 = seq.last_aa * 3
        if key in new_sequences:
            seqdata['NASequence'] = new_sequences[key][first_na0:last_na0]
        writer.writerow(seqdata)
