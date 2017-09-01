#! /usr/bin/env python
import sys

with open(sys.argv[1]) as fp:
    skip = True
    for line in fp:
        if skip:
            if line.startswith('Pair, dN/dS'):
                skip = False
            else:
                continue
        print('\t'.join(s.strip() for s in line.split(',')))
