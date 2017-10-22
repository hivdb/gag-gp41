#! /usr/bin/env python

import os
from collections import Counter
from data_reader import data_reader, ROOT


distmt = data_reader(os.path.join(
    ROOT, 'local', 'hyphyOutput', 'gp41NaiveDistanceMatrix.txt'
), delimiter='\t')
counter = Counter()
for p in distmt:
    if float(p['Distance']) < 0.3:
        continue
    counter[p['Sequence1']] += 1
    counter[p['Sequence2']] += 1

for seq, n in counter.most_common(30):
    if n > 5:
        print(seq, n)
