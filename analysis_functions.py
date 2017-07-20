import math
from itertools import groupby
from decimal import Decimal

from data_reader import get_prevalence
from data_reader import sample_reader, sequence_reader


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


def compare_codons(gene, prev_codon, post_codon):
    prev_aa = set(prev_codon.aa)
    post_aa = set(post_codon.aa)
    if prev_aa == post_aa or post_aa.issubset(prev_aa):
        return
    post_aa -= prev_aa
    if not post_aa or post_aa == {'*'}:
        return
    prev_aa = ''.join(prev_aa)
    post_aa = ''.join(post_aa)

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


def aa_changes_per_person(gene, category, group):
    if isinstance(category, str):
        samples = list(
            sample_reader(lambda spl: spl.category == category))
    else:
        samples = list(
            sample_reader(lambda spl: spl.category in category))
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
        start_aa = max(prev_seq.first_aa, post_seq.first_aa)
        end_aa = min(prev_seq.last_aa, post_seq.last_aa)
        prev_codons = prev_seq.iter_codons(start_aa, end_aa)
        post_codons = post_seq.iter_codons(start_aa, end_aa)

        for prev_codon, post_codon in zip(prev_codons, post_codons):
            result = compare_codons(gene, prev_codon, post_codon)
            if not result:
                continue
            yield {
                'PID': pid,
                'Group': group,
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
