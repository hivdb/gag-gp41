#! /usr/bin/env python

import sys
import json

GENES = {
    'GAG': 500,
    'GP41': 345
}

gene = sys.argv[1]
input_filename = sys.argv[2]
output_filename = sys.argv[3]
genelen = GENES[gene]
with open(input_filename) as infp, open(output_filename, 'w') as outfp:
    data = json.load(infp)
    all_nas = []
    for seq in data[gene]:
        leftpad = (seq['Report']['FirstAA'] - 1) * 3
        rightpad = (genelen - seq['Report']['LastAA']) * 3
        ctrl = seq['Report']['ControlLine']
        naseq = seq['Report']['NucleicAcidsLine'].replace(' ', '-')
        ctrl = '-' * leftpad + ctrl + '-' * rightpad
        naseq = '-' * leftpad + naseq + '-' * rightpad
        nas = [''] * (genelen * 3 + 1)
        pos = 1
        for c, n in zip(ctrl, naseq):
            if c == '+':
                nas[pos - 1] += n
            else:
                nas[pos] = n
                pos += 1
        all_nas.append(nas)
    result_naseqs = [''] * len(all_nas)
    for pos_nas in zip(*all_nas):
        pos_nas = list(pos_nas)
        maxlen = len(max(pos_nas, key=lambda nas: len(nas)))
        for idx, nas in enumerate(pos_nas):
            result_naseqs[idx] += nas + '-' * (maxlen - len(nas))
    for seq, naseq in zip(data[gene], result_naseqs):
        outfp.write('>{}\n'.format(seq['Name']))
        outfp.write(naseq + '\n')
