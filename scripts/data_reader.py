import os
import csv
from datetime import datetime
from collections import namedtuple

from codonutils import translate_codon

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

with open(os.path.join(ROOT, 'data', 'consensus.csv')) as fp:
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


PREVALENCE = None


def get_prevalence(gene, pos, aa):
    global PREVALENCE
    if not PREVALENCE:
        PREVALENCE = {}
        for _gene in ('gag', 'gp41'):
            with open(os.path.join(
                    ROOT, 'data', 'naiveStudies',
                    '{}AAPrevalence.csv'.format(_gene))) as fp:
                for p in csv.DictReader(fp):
                    if p['Subtype']:
                        continue
                    _aa = (p['AA']
                           .replace('ins', 'i')
                           .replace('del', 'd'))
                    PREVALENCE[
                        (p['Gene'], int(p['Pos']), _aa)] = MutPrevalence(p)
    prev = PREVALENCE.get((gene, pos, aa))
    if prev is None:
        return .0
    else:
        return prev.percent


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

    def iter_codons(self, start_aa=None, end_aa=None, frameshift=0):
        cons = CONSENSUS[self.gene]['AASeq']
        fs = frameshift
        for i in range(0, self.last_aa - self.first_aa + 1):
            aa_pos = i + self.first_aa
            if (start_aa and aa_pos < start_aa) or \
                    (end_aa and aa_pos > end_aa):
                continue
            triplet = self.na_sequence[i * 3 + fs:i * 3 + 3 + fs]
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


def data_reader(filepath, decorator=dict, filter_func=None, delimiter=','):
    with open(filepath) as fp:
        bom = fp.read(1)
        if bom != '\ufeff':
            fp.seek(0)
        items = csv.DictReader(fp, delimiter=delimiter)
        for item in items:
            item = decorator(item)
            if not filter_func or filter_func(item):
                yield item


def sequence_reader(filter_func=None):
    return data_reader(
        os.path.join(ROOT, 'data', 'sequences.csv'),
        Sequence, filter_func)


def possible_apobecs_reader(gene, filter_func=None):
    filename = os.path.join(
        ROOT, 'data', 'naiveStudies', 'apobec',
        '{}PossibleApobecs.csv'.format(gene.lower())
    )
    if not os.path.exists(filename):
        return []
    return data_reader(
        filename,
        filter_func=filter_func
    )


def sample_reader(filter_func=None):
    return data_reader(
        os.path.join(ROOT, 'data', 'samples.csv'),
        Sample, filter_func)


def naive_sequence_reader(gene, filter_func=None):
    return data_reader(
        os.path.join(
            ROOT, 'internalFiles', 'naiveStudies',
            '{}.csv'.format(gene.lower())),
        filter_func=filter_func)


def fasta_reader(gene, rx):
    filename = os.path.join(
        ROOT, 'data', 'fasta', '{}{}.aln.fasta.txt'.format(gene, rx))
    return any_fasta_reader(filename)


def any_fasta_reader(filename):
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
