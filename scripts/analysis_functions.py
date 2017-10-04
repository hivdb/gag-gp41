import math
from itertools import groupby
from decimal import Decimal
from collections import Counter, OrderedDict

from data_reader import get_prevalence
from data_reader import (
    CONSENSUS, sample_reader,
    sequence_reader, naive_sequence_reader
)


PREC2 = Decimal('1.00')
PREC3 = Decimal('1.000')


def calcfold(left, right):
    fold = -999.0
    if right == 0:
        if left < .1:
            fold = 1.0
        elif left < 1:
            fold = 10.0
        elif left < 10:
            fold = 100.0
        else:
            fold = 1000.0
    else:
        fold = left / right
    fold = max(fold, 0.001)
    fold = min(fold, 1000)
    return fold


def compare_aa_change(gene, prev_codon, post_codon):
    prev_aa = set(prev_codon.aa)
    post_aa = set(post_codon.aa)
    if prev_aa == post_aa or post_aa.issubset(prev_aa):
        return
    post_aa -= prev_aa
    if not post_aa or post_aa == {'*'}:
        return
    prev_aa = ''.join(sorted(prev_aa))
    post_aa = ''.join(sorted(post_aa))

    pos = prev_codon.position
    if len(prev_aa) == 1 and len(post_aa) == 1:
        pre_prev = get_prevalence(gene, pos, prev_aa)
        post_prev = get_prevalence(gene, pos, post_aa)
        fold = calcfold(pre_prev, post_prev)
        logfold = math.log10(fold)
        pre_prev = Decimal(pre_prev).quantize(PREC3)
        post_prev = Decimal(post_prev).quantize(PREC3)
        fold = Decimal(fold).quantize(PREC3)
        logfold = Decimal(logfold).quantize(PREC2)
    else:
        pre_prev = post_prev = fold = logfold = 'NA'

    return {
        'Pos': pos,
        'PreAA': prev_aa,
        'PostAA': post_aa,
        'PrePrev': pre_prev,
        'PostPrev': post_prev,
        'Fold': fold,
        'LogFold': logfold,
    }


def compare_codon_change(gene, prev_codon, post_codon):
    prev_triplet = prev_codon.triplet
    post_triplet = post_codon.triplet
    if prev_triplet == post_triplet or \
            prev_triplet in ('---', '...') or \
            post_triplet in ('---', '...'):
        return
    type_ = 'syn'
    if prev_codon.aa != post_codon.aa:
        type_ = 'non'
    num_changes = sum([bp0 != bp1 for bp0, bp1
                       in zip(prev_triplet, post_triplet)])

    return {
        'Pos': prev_codon.position,
        'Type': type_,
        'Codons': '{}-->{}'.format(prev_triplet, post_triplet),
        'NumNAChanges': num_changes,
        'AAs': (prev_codon.aa if type_ == 'syn'
                else '{}-->{}'.format(prev_codon.aa, post_codon.aa)),
    }


def aggregate_mut_prevalence(gene):
    result = OrderedDict()
    total = Counter()
    all_aas = 'ACDEFGHIKLMNPQRSTVWY~X#*'
    consensus = CONSENSUS[gene]['AASeq']
    sequences = list(naive_sequence_reader(gene))
    for pos in range(1, len(consensus) + 1):
        cons = consensus[pos - 1]
        for aa in all_aas:
            result[(pos, aa)] = 0
        for sequence in sequences:
            aa = sequence['P{}'.format(pos)]
            if aa == '.':
                continue
            aa = (aa
                  .replace(cons, '')
                  .replace('d', '~')
                  .replace('-', cons))
            if len(aa) == 1:
                result[(pos, aa)] += 1
            else:
                result[(pos, 'X')] += 1
            total[pos] += 1
    return [{
        'Gene': gene,
        'Pos': pos,
        'AA': aa,
        'Pcnt': Decimal(count / total[pos] * 100).quantize(PREC3),
        'Count': count,
        'PosTotal': total[pos]
    } for (pos, aa), count in result.items()]


def aggregate_naiveseqs_stat(gene):
    consensus = CONSENSUS[gene]['AASeq']
    sequences = list(naive_sequence_reader(gene))
    result = []

    for sequence in sequences:
        diffs = 0
        stopcodons = 0
        for pos in range(1, len(consensus) + 1):
            aa = sequence['P{}'.format(pos)]
            if aa == '.':
                continue
            if aa != '-':
                diffs += 1
            if '*' in aa:
                stopcodons += 1
        result.append({
            'Accession': sequence['Accession'],
            'PMID': sequence['PMID'],
            'Gene': gene,
            'Subtype': sequence['lanlSubtype'],
            'NumAAChanges': diffs,
            'NumStopCodons': stopcodons
        })

    return result


def aggregate_naiveseqs_posstat(gene):
    consensus = CONSENSUS[gene]['AASeq']
    sequences = list(naive_sequence_reader(gene))
    result = []

    for pos in range(1, len(consensus) + 1):
        diffs = 0
        stopcodons = 0
        for sequence in sequences:
            aa = sequence['P{}'.format(pos)]
            if aa == '.':
                continue
            if aa != '-':
                diffs += 1
            if '*' in aa:
                stopcodons += 1
        result.append({
            'AAPosition': pos,
            'Gene': gene,
            'NumAAChanges': diffs,
            'NumStopCodons': stopcodons
        })

    return result


def iter_sequence_pairs(gene, category=None):
    if gene == 'Gp41':
        # a stupid fix
        gene = 'gp41'

    if isinstance(category, str):
        samples = list(
            sample_reader(lambda spl: spl.category == category))
    elif hasattr(category, '__contains__'):
        samples = list(
            sample_reader(lambda spl: spl.category in category))
    else:
        samples = list(sample_reader())

    pids = {spl.pid for spl in samples}
    sequences = sequence_reader(
        lambda seq: seq.gene == gene and seq.pid in pids)

    def sortkeyfunc(item):
        return item.pid, item.time_point

    def groupkeyfunc(item):
        return item.pid

    samples = sorted(samples, key=sortkeyfunc)
    samples = {k: list(l) for k, l in groupby(samples, groupkeyfunc)}
    sequences = sorted(sequences, key=sortkeyfunc)
    sequences = {k: list(l) for k, l in groupby(sequences, groupkeyfunc)}

    for pid in sorted(sequences.keys()):
        post_spl, prev_spl = samples[pid]
        post_seq, prev_seq = sequences[pid]
        yield pid, prev_spl.category, prev_seq, post_seq


def iter_codon_pairs(gene, category=None):

    for pid, cat, prev_seq, post_seq in iter_sequence_pairs(gene, category):
        start_aa = max(prev_seq.first_aa, post_seq.first_aa)
        end_aa = min(prev_seq.last_aa, post_seq.last_aa)
        prev_codons = prev_seq.iter_codons(start_aa, end_aa)
        post_codons = post_seq.iter_codons(start_aa, end_aa)

        for prev_codon, post_codon in zip(prev_codons, post_codons):
            yield pid, cat, prev_codon, post_codon


def aa_changes_per_person(gene, category, group):
    for pid, _, prev_codon, post_codon in iter_codon_pairs(gene, category):
        result = compare_aa_change(gene, prev_codon, post_codon)
        if not result:
            continue
        yield {
            'PID': pid,
            'Group': group,
            **result
        }


def codon_changes_per_person(gene, category):
    for pid, rx, prev_codon, post_codon \
            in iter_codon_pairs(gene, category):
        result = compare_codon_change(gene, prev_codon, post_codon)
        if not result:
            continue
        yield {
            'PID': pid,
            'Rx': rx,
            **result
        }


def aggregate_aa_changes_by_pos(gene, category, group):
    aa_changes = aa_changes_per_person(gene, category, group)

    def keyfunc(c):
        return c['Pos'], c['PostAA'], c['PreAA']

    aa_changes = sorted(aa_changes, key=keyfunc)
    aa_changes = groupby(aa_changes, keyfunc)

    for _, changes in aa_changes:
        changes = list(changes)
        num_pts = len(changes)
        item = {
            **changes[0],
            'NumPts': num_pts
        }
        item.pop('PID')
        yield item
