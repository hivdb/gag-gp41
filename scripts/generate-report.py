import os
import csv
import math
from itertools import groupby
from collections import namedtuple

from numpy import median
from tabulate import tabulate
from scipy.stats import mannwhitneyu

APPDIR = '/app'
HYPHYOUT = os.path.join(APPDIR, 'result_data/hyphy_output')
CLEANOUT = os.path.join(APPDIR, 'result_data/hyphy_cleaned_output')
NAN = float('nan')


PairwiseResult = namedtuple(
    'PairwiseResult',
    ('pvalue', 'leftmin', 'leftmed', 'leftmax',
     'rightmin', 'rightmed', 'rightmax')
)


def pairwise_result(title, leftfname, rightfname):
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
        '![Gag sites](https://github.com/hivdb/'
        'gag-gp41/raw/master/report/gag-mutations.png)\n'
    )
    print(
        '![Gp41 sites](https://github.com/hivdb/'
        'gag-gp41/raw/master/report/gp41-mutations.png)\n'
    )
    print('## Pairwise dN/dS ratio comparison\n')
    for gene in ('Gag', 'Gp41'):
        pairwise_result(
            'Complete {}'.format(gene),
            os.path.join(CLEANOUT, '{}PIs.pairwise.tsv'.format(gene)),
            os.path.join(CLEANOUT, '{}NNRTIs.pairwise.tsv'.format(gene))
        )
        if gene == 'Gag':
            pairwise_result(
                '{} MA Domain'.format(gene),
                os.path.join(CLEANOUT, '{}PIs-MA.pairwise.tsv'.format(gene)),
                os.path.join(CLEANOUT, '{}NNRTIs-MA.pairwise.tsv'.format(gene))
            )
            pairwise_result(
                '{} CA Domain'.format(gene),
                os.path.join(CLEANOUT, '{}PIs-CA.pairwise.tsv'.format(gene)),
                os.path.join(CLEANOUT, '{}NNRTIs-CA.pairwise.tsv'.format(gene))
            )
        else:
            pairwise_result(
                '{} CD Domain'.format(gene),
                os.path.join(CLEANOUT, '{}PIs-CD.pairwise.tsv'.format(gene)),
                os.path.join(CLEANOUT, '{}NNRTIs-CD.pairwise.tsv'.format(gene))
            )
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
