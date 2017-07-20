import math
from itertools import groupby
from decimal import Decimal

from data_reader import get_prevalence
from data_reader import sample_reader, sequence_reader


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

    prec2 = Decimal('1.00')
    prec3 = Decimal('1.000')

    for pid in sorted(sequences.keys()):
        post_spl, prev_spl = samples[pid]
        post_seq, prev_seq = sequences[pid]
        start_aa = max(prev_seq.first_aa, post_seq.first_aa)
        end_aa = min(prev_seq.last_aa, post_seq.last_aa)
        prev_codons = prev_seq.iter_codons(start_aa, end_aa)
        post_codons = post_seq.iter_codons(start_aa, end_aa)

        for prev_codon, post_codon in zip(prev_codons, post_codons):
            prev_aa = set(prev_codon.aa)
            post_aa = set(post_codon.aa)
            if prev_aa == post_aa or post_aa.issubset(prev_aa):
                continue
            post_aa -= prev_aa
            if not post_aa or post_aa == {'*'}:
                continue
            prev_aa = ''.join(prev_aa)
            post_aa = ''.join(post_aa)

            pos = prev_codon.position
            if len(prev_aa) == 1 and len(post_aa) == 1:
                pre_prev = get_prevalence(gene, pos, prev_aa)
                post_prev = get_prevalence(gene, pos, post_aa)
                fold = calcfold(pre_prev, post_prev)
                logfold = math.log10(fold)
                pre_prev = Decimal(pre_prev).quantize(prec3)
                post_prev = Decimal(post_prev).quantize(prec3)
                fold = Decimal(fold).quantize(prec3)
                logfold = Decimal(logfold).quantize(prec2)
            else:
                pre_prev = post_prev = fold = logfold = 'NA'
            yield {
                'PID': pid,
                'Group': group,
                'Pos': pos,
                'PreAA': prev_aa,
                'PostAA': post_aa,
                'PrePrev': pre_prev,
                'PostPrev': post_prev,
                'Fold': fold,
                'LogFold': logfold,
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
