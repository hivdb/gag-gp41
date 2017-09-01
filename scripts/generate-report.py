import os
import csv
import math
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

    print('- {}:'.format(title))
    print('  - PIs (treated): {:.2f}, range {:.2f} to {:.2f}'.format(
        result.leftmed, result.leftmin, result.leftmax
    ))
    print('  - NNRTIs (treated): {:.2f}, range {:.2f} to {:.2f}'.format(
        result.rightmed, result.rightmin, result.rightmax
    ))
    print('  - P value: {:.1f}'.format(result.pvalue))


def fel_result(title, fname):
    with open(fname) as fp:
        data = csv.DictReader(fp, delimiter='\t')
        r = [(one['Codon'], float(one['Selection detected?'].split(' = ')[-1]))
             for one in data
             if one['Selection detected?'].startswith('Pos.')]
        r = [one for one in r if one[1] < 0.05]
    print('### {} (P <= 0.05)\n'.format(title))
    print(tabulate(r, ['Position', 'P value'], tablefmt='pipe'))
    print()


if __name__ == '__main__':
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
    for gene in ('Gag', 'Gp41'):
        for rx in ('PIs', 'NNRTIs'):
            fel_result(
                '{} - {}'.format(gene, rx),
                os.path.join(CLEANOUT, '{}{}.fel.tsv'.format(gene, rx))
            )
