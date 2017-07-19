import os
import csv
from datetime import datetime
from collections import namedtuple

from codonutils import translate_codon

PWD = os.path.dirname(__file__)

with open(os.path.join(PWD, 'data', 'consensus.csv')) as fp:
    CONSENSUS = {c['Gene']: c for c in csv.DictReader(fp)}


class Codon(namedtuple('Codon', ['position', 'gene',
                                 'triplet', 'cons_aa'])):

    @property
    def aa(self):
        if self.triplet == '---':
            return '-'
        aa = translate_codon(self.triplet)
        if aa == self.cons_aa:
            return '.'
        return aa

    @property
    def is_mutation(self):
        return self.aa != '.'

    @property
    def is_deletion(self):
        return self.aa != '-'


class Sequence:

    def __init__(self, seqdata):
        self._data = seqdata

    @property
    def pid(self):
        return self._data['PID']

    @property
    def time_point(self):
        return self._data['TimePoint']

    @property
    def date(self):
        return datetime.strptime(self._data['Date'], '%Y-%m-%d')

    @property
    def gene(self):
        return self._data['Gene']

    @property
    def first_aa(self):
        return int(self._data['FirstAA'])

    @property
    def last_aa(self):
        return int(self._data['LastAA'])

    @property
    def na_sequence(self):
        return self._data['NASequence']

    def iter_codons(self):
        cons = CONSENSUS[self.gene]['AASeq']
        for i in range(0, self.last_aa - self.first_aa + 1):
            triplet = self.na_sequence[i * 3:i * 3 + 3]
            cons_aa = cons[i + self.first_aa - 1]
            yield Codon(
                position=i + self.first_aa,
                gene=self.gene,
                triplet=triplet,
                cons_aa=cons_aa)


def sequence_reader(filter_func=None):
    with open(os.path.join(PWD, 'data', 'sequences.csv')) as fp:
        sequences = csv.DictReader(fp)
        for seq in sequences:
            seq = Sequence(seq)
            if not filter_func or filter_func(seq):
                yield seq
