#! /usr/bin/env python
import sys

with open(sys.argv[1]) as fp:
    skip = True
    for line in fp:
        if skip:
            if line.startswith('|'):
                skip = False
            else:
                continue
        elif line.startswith('|:'):
            continue
        elif not line.strip():
            skip = True
            continue
        line = line.strip('|\r\n\t ')
        print('\t'.join(s.strip() for s in line.split('|')))
