import os
import csv
import math
from itertools import groupby
from collections import namedtuple

from numpy import median
from tabulate import tabulate
from numpy import percentile
from scipy.stats import mannwhitneyu

from data_reader import fasta_reader

APPDIR = '/app'
FASTAS = os.path.join(APPDIR, 'data/fasta')
HYPHYOUT = os.path.join(APPDIR, 'result_data/hyphy_output')
CLEANOUT = os.path.join(APPDIR, 'result_data/hyphy_cleaned_output')
NAN = float('nan')
GENE_LENGTHS = {'Gag': 500, 'Gp41': 345}


PairwiseResult = namedtuple(
    'PairwiseResult',
    ('pvalue', 'leftmin', 'leftmed', 'leftmax',
     'rightmin', 'rightmed', 'rightmax')
)


def dnds_result(title, leftfname, rightfname):
    with open(leftfname) as leftfp, open(rightfname) as rightfp:
        left = csv.DictReader(leftfp, delimiter='\t')
        right = csv.DictReader(rightfp, delimiter='\t')
        left = [float(one['dN/dS']) for one in left]
        right = [float(one['dN/dS']) for one in right]

        # remove NaNs
        left = [n for n in left if not math.isnan(n)]
        right = [n for n in right if not math.isnan(n)]

        # perform two-sample Wilcoxon test, aka ‘Mann-Whitney’ test
        # see: https://stackoverflow.com/a/33890615/2644759
        pvalue = mannwhitneyu(left, right).pvalue

        result = PairwiseResult(
            pvalue=pvalue,
            leftmin=min(left),
            leftmed=median(left),
            leftmax=max(left),
            rightmin=min(right),
            rightmed=median(right),
            rightmax=max(right)
        )

    print('### {}\n'.format(title))
    print('- PIs (treated): {:.2f}, range {:.2f} to {:.2f}'.format(
        result.leftmed, result.leftmin, result.leftmax
    ))
    print('- NNRTIs (treated): {:.2f}, range {:.2f} to {:.2f}'.format(
        result.rightmed, result.rightmin, result.rightmax
    ))
    print('- P value: {:.1f}'.format(result.pvalue))
    print()


def summarize_dnds(gene, rx, domain=None):
    if domain:
        domain = '-' + domain
    else:
        domain = ''
    filename = os.path.join(
        CLEANOUT, '{}{}{}.pairwise.tsv'.format(gene, rx, domain))
    with open(filename) as fp:
        dnds = csv.DictReader(fp, delimiter='\t')
        dnds = [float(one['dN/dS']) for one in dnds]

        # remove NaNs
        dnds = [n for n in dnds if not math.isnan(n)]

        q25, q50, q75 = percentile(dnds, (25, 50, 75))
        return '{:.2f} ({:.2f}-{:.2f})'.format(q50, q25, q75), dnds


def summarize_na_diffs(gene, rx):
    cchanges = get_codon_changes(gene)
    cchanges = [cc for cc in cchanges if cc['Rx'] == rx]
    na_length = GENE_LENGTHS[gene] * 3
    perpt = groupby(cchanges, lambda cc: cc['PID'])
    perpt = [sum(int(cc['NumNAChanges'])
                 for cc in ccs) * 100 / na_length for _, ccs in perpt]
    q25, q50, q75 = percentile(perpt, (25, 50, 75))
    return '{:.1f} ({:.1f}-{:.1f}%)'.format(q50, q25, q75)


def summarize_aa_diffs(gene, rx):
    cchanges = get_codon_changes(gene)
    cchanges = [
        cc for cc in cchanges if cc['Rx'] == rx and cc['Type'] == 'syn']
    aa_length = GENE_LENGTHS[gene]
    perpt = groupby(cchanges, lambda cc: cc['PID'])
    perpt = [len(list(ccs)) * 100 / aa_length for _, ccs in perpt]
    q25, q50, q75 = percentile(perpt, (25, 50, 75))
    return '{:.1f} ({:.1f}-{:.1f}%)'.format(q50, q25, q75)


def summarize_ambiguities(gene, rx):
    perpt = {}
    na_length = GENE_LENGTHS[gene] * 3
    for header, sequence in fasta_reader(gene, rx):
        pt = header.split('|', 1)[1].rsplit('_', 1)[0]
        ts = header.rsplit('_', 1)[1]
        count = len([na for na in sequence.upper() if na not in 'ACTG-.'])
        perpt.setdefault(pt, [NAN, NAN])
        perpt[pt][('Pre', 'Post').index(ts)] = count * 100 / na_length
    delta = [pre - post for pre, post in perpt.values()]
    return (
        '{:.2f} ({:.2f}-{:.2f}%)'.format(
            *percentile([pre for pre, _ in perpt.values()], (50, 25, 75))),
        '{:.2f} ({:.2f}-{:.2f}%)'.format(
            *percentile([post for _, post in perpt.values()], (50, 25, 75))),
        '{:.2f} ({:.2f}-{:.2f}%)'.format(
            *percentile(delta, (50, 25, 75)))
    )


def summarize_gag_cleavage_sites():

    sites = (
        132, 373, 374, 375, 376, 378, 380, 381,
        429, 436, 451, 453, 484, 485, 490
    )
    result = []
    with open(os.path.join(
            APPDIR, 'result_data',
            'GagAAChangesByPosWPrev.csv'.format(gene))
    ) as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            if int(row['Pos']) in sites:
                result.append([
                    row['Group'], row['Pos'],
                    '{PreAA}=&gt;{PostAA}'.format(**row),
                    row['NumPts']])
    return result


def fel_result(title, fname):
    with open(fname) as fp:
        data = csv.DictReader(fp, delimiter='\t')
        r = [['~', one['Codon'],
              float(one['Selection detected?'].split(' = ')[-1])]
             for one in data
             if one['Selection detected?'].startswith('Pos.')]
        r = [one for one in r if one[2] < 0.05]
    if r:
        r[0][0] = title
    return r


def get_codon_changes(gene):
    with open(os.path.join(
            APPDIR, 'result_data',
            '{}CodonChangesByPt.csv'.format(gene))
    ) as fp:
        reader = csv.DictReader(fp)
        return list(reader)


def meds_result(gene, rx):
    fname = os.path.join(HYPHYOUT, '{}{}.meds.result.csv'.format(gene, rx))
    cchanges = get_codon_changes(gene)
    medsdata = []
    with open(fname) as fp:
        skip = True
        for line in fp:
            if skip:
                if line.startswith('MEDS'):
                    skip = False
                continue
            elif not line.strip():
                skip = True
            medsdata.append(line)
    reader = csv.DictReader(medsdata)
    results = []
    for pos, rows in groupby(reader, lambda r: r['Site']):
        aas = [r['AA'] for r in rows]
        r = []
        for c in cchanges:
            if c['Rx'] != rx or c['Type'] == 'syn' or c['Pos'] != pos:
                continue
            aapre, aapost = c['AAs'].split('-->')
            aapost = set(aapost) - set(aapre)
            if not aapost:
                continue
            aapost = ''.join(sorted(aapost))
            if set(aas) & set(aapre + aapost):
                r.append((aapre + '=&gt;' + aapost, c['PID']))
        if not r:
            continue
        r = groupby(sorted(r, key=lambda i: i[0]), lambda i: i[0])
        r = [['~', '~', m, len(list(p))] for m, p in r]
        r[0][1] = '{}{}'.format(pos, '/'.join(aas))
        results.extend(r)
    if results:
        results[0][0] = '{} - {}'.format(gene, rx)
    else:
        results = [['{} - {}'.format(gene, rx), 'None', '-', '-']]
    return results


if __name__ == '__main__':
    print(
        '## Graphical summary of gag/gp41 sites '
        'at which amino acid mutations developed during therapy\n'
    )
    print(
        '### Gag sites\n\n'
        '![Gag sites](https://github.com/hivdb/'
        'gag-gp41/raw/master/report/gag-mutations.png)\n'
    )
    print(
        '### Gp41 sites\n\n'
        '![Gp41 sites](https://github.com/hivdb/'
        'gag-gp41/raw/master/report/gp41-mutations.png)\n'
    )
    print('## Pairwise comparison\n')
    for gene in ('Gag', 'Gp41'):
        pairwise = []
        print('### {}\n'.format(gene))

        for rx in ('NNRTIs', 'PIs'):
            pairwise.append([
                '{} NA differences'.format(rx),
                summarize_na_diffs(gene, rx)])

        for rx in ('NNRTIs', 'PIs'):
            pairwise.append([
                '{} AA differences'.format(rx),
                summarize_aa_diffs(gene, rx)])

        all_dndslist = []
        for rx in ('NNRTIs', 'PIs'):
            ratio, dndslist = summarize_dnds(gene, rx)
            all_dndslist.append(dndslist)
            pairwise.append(['{} dN/dS ratio'.format(rx), ratio])
        # perform two-sample Wilcoxon test, aka ‘Mann-Whitney’ test
        # see: https://stackoverflow.com/a/33890615/2644759
        pvalue = mannwhitneyu(*all_dndslist).pvalue
        pairwise.append(['P value of dN/dS ratio', '{:.1f}'.format(pvalue)])

        for rx in ('NNRTIs', 'PIs'):
            ab_pre, ab_post, ab_delta = summarize_ambiguities(gene, rx)
            pairwise.append(['{} ambiguities (pre)'.format(rx), ab_pre])
            pairwise.append(['{} ambiguities (post)'.format(rx), ab_post])
            pairwise.append(['{} ambiguities (delta)'.format(rx), ab_delta])

        print(tabulate(pairwise, ['Title', 'Value'], tablefmt='pipe'))
        print()

    print('\n## Gag cleavage sites\n')
    print(tabulate(
        summarize_gag_cleavage_sites(),
        ['Rx', 'Position', 'AA change', '# patients'], tablefmt='pipe'))

    print('\n## Positions with evidence for diversifying selection (FEL)\n')
    fel = []
    for gene in ('Gag', 'Gp41'):
        for rx in ('PIs', 'NNRTIs'):
            fel.extend(fel_result(
                '{} - {}'.format(gene, rx),
                os.path.join(CLEANOUT, '{}{}.fel.tsv'.format(gene, rx))
            ))
    print(tabulate(fel, ['Group', 'Position', 'P value'], tablefmt='pipe'))
    print('\n\n## Positions with evidence for directional selection (MEDS)\n')
    meds = []
    for gene in ('Gag', 'Gp41'):
        for rx in ('PIs', 'NNRTIs'):
            meds.extend(meds_result(gene, rx))

    print(tabulate(
        meds, ['Group', 'Mutation', 'Detail', '# patients'], tablefmt='pipe'))
    print()
