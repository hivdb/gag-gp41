import os
import csv
from datetime import datetime
from collections import namedtuple

from codonutils import translate_codon

PWD = os.path.dirname(__file__)

with open(os.path.join(PWD, 'data', 'consensus.csv')) as fp:
    CONSENSUS = {c['Gene']: c for c in csv.DictReader(fp)}


class MutPrevalence:

    def __init__(self, data):
        self._data = data

    @property
    def gene(self):
        return self._data['Gene']

    @property
    def position(self):
        return int(self._data['Pos'])

    @property
    def aa(self):
        return self._data['aa']

    @property
    def percent(self):
        return float(self._data['Pcnt'])

    @property
    def count(self):
        return int(self._data['Count'])

    @property
    def total(self):
        return int(self._data['PosTotal'])


with open(os.path.join(PWD, 'data', 'mut_prevalence.csv')) as fp:
    PREVALENCE = {(p['Gene'], int(p['Pos']), p['AA']):
                  MutPrevalence(p) for p in csv.DictReader(fp)}


def get_prevalence(gene, pos, aa):
    prev = PREVALENCE.get((gene, pos, aa))
    if prev:
        return prev.percent
    else:
        return .0


class Codon(namedtuple('Codon', ['position', 'gene',
                                 'triplet', 'cons_aa'])):

    @property
    def aa(self):
        if self.triplet == '---':
            return '-'
        return translate_codon(self.triplet)

    @property
    def is_mutation(self):
        return self.aa != self.cons_aa

    @property
    def is_deletion(self):
        return self.aa != '-'

    @property
    def prevalence(self):
        return get_prevalence(self.gene, self.position, self.aa)


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
        return datetime.strptime(self._data['Date'], '%Y-%m-%d').date()

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

    def iter_codons(self, start_aa=None, end_aa=None):
        cons = CONSENSUS[self.gene]['AASeq']
        for i in range(0, self.last_aa - self.first_aa + 1):
            aa_pos = i + self.first_aa
            if (start_aa and aa_pos < start_aa) or \
                    (end_aa and aa_pos > end_aa):
                continue
            triplet = self.na_sequence[i * 3:i * 3 + 3]
            cons_aa = cons[aa_pos - 1]
            yield Codon(
                position=aa_pos,
                gene=self.gene,
                triplet=triplet,
                cons_aa=cons_aa)


class Sample:

    def __init__(self, sample_data):
        self._data = sample_data

    @property
    def pid(self):
        return self._data['PID']

    @property
    def time_point(self):
        return self._data['TimePoint']

    @property
    def date(self):
        return datetime.strptime(self._data['Date'], '%Y-%m-%d').date()

    @property
    def category(self):
        return self._data['Category']


def data_reader(filepath, decorator, filter_func=None):
    with open(filepath) as fp:
        items = csv.DictReader(fp)
        for item in items:
            item = decorator(item)
            if not filter_func or filter_func(item):
                yield item


def sequence_reader(filter_func=None):
    return data_reader(
        os.path.join(PWD, 'data', 'sequences.csv'),
        Sequence, filter_func)


def sample_reader(filter_func=None):
    return data_reader(
        os.path.join(PWD, 'data', 'samples.csv'),
        Sample, filter_func)
