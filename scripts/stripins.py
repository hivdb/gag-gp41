#! /usr/bin/env python
from __future__ import print_function

import sys
import csv
from codonutils import translate_codons

INSHEADER = ['Header', 'Gene', 'Position', 'InsertedAAs', 'InsertedNAs']


def fasta_reader(filename):
    with open(filename) as fp:
        header = None
        seq = []
        for line in fp:
            if line.startswith('#'):
                continue
            elif line.startswith('>'):
                if seq:
                    yield header, ''.join(seq)
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if seq:
            yield header, ''.join(seq)


def stripins(gene, refheader, infile):
    sequences = list(fasta_reader(infile))
    for header, seq in sequences:
        if header.endswith(refheader):
            refseq = seq
            break
    for header, seq in sequences:
        # skip consensus
        if header.endswith(refheader):
            continue
        result_seq = []
        insertions = []
        i = 0
        insnas = []
        for r, n in zip(refseq, seq):
            if r == '-' and i % 3 == 0:
                if n != '-':
                    insnas.append(n)
            else:
                if insnas:
                    while len(insnas) % 3 != 0:
                        # remove unwanted frameshifts
                        insnas.pop()
                    insnas = ''.join(insnas)
                    insaas = translate_codons(insnas)

                    insertions.append({
                        'Header': header,
                        'Gene': gene,
                        'Position': i // 3,
                        'InsertedAAs': insaas,
                        'InsertedNAs': insnas
                    })
                    insnas = []
                i += 1
                result_seq.append(n)
        yield header, ''.join(result_seq), insertions


def main():
    if len(sys.argv) != 6:
        print('Usage: {} <GENE_NAME> <REF_HEADER> '
              '<INPUT FASTA> <OUTPUT FASTA> <OUTPUT INSERTIONS>'
              .format(sys.argv[0]), file=sys.stderr)
        exit(1)
    gene, refheader, infn, outfn, insfn = sys.argv[1:]
    with open(outfn, 'w') as outfile, open(insfn, 'w') as insfile:
        writer = csv.DictWriter(insfile, INSHEADER)
        writer.writeheader()
        for header, seq, insertions in stripins(gene, refheader, infn):
            outfile.write('>{}\n{}\n'.format(header, seq))
            writer.writerows(insertions)


if __name__ == '__main__':
    main()
