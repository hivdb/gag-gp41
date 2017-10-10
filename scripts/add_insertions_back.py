#! /usr/bin/env python
import sys

from itertools import groupby
from data_reader import data_reader, any_fasta_reader


def sortkey(ins):
    return (
        ins['Gene'],
        int(ins['PID']),
        ins['Timepoint'],
        -int(ins['Pos'])
    )


def groupkey(ins):
    return (
        ins['Gene'],
        int(ins['PID']),
        ins['Timepoint']
    )


gene = sys.argv[1]

all_insertions = data_reader(sys.argv[3])
all_insertions = sorted(all_insertions, key=sortkey)
all_insertions = groupby(all_insertions, groupkey)
all_insertions = {key: list(val) for key, val in all_insertions}

for header, seq in any_fasta_reader(sys.argv[2]):
    ptid, _, ts = header.split('|', 1)[1].split('_', 2)
    key = (gene, int(ptid), ts)
    if key in all_insertions:
        insertions = all_insertions[key]
        for ins in insertions:
            inspos = int(ins['Pos']) * 3
            insnas = ins['NA']
            seq = seq[:inspos] + insnas + seq[inspos:]
    print('>{}\n{}'.format(header, seq))
