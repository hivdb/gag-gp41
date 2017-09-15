import os
import csv
import math
from itertools import groupby

from tabulate import tabulate
from numpy import percentile
from scipy.stats import mannwhitneyu

APPDIR = '/app'
FASTAS = os.path.join(APPDIR, 'data/fasta')
HYPHYOUT = os.path.join(APPDIR, 'result_data/hyphy_output')
CLEANOUT = os.path.join(APPDIR, 'result_data/hyphy_cleaned_output')
NAN = float('nan')
GENE_LENGTHS = {'Gag': 500, 'Gp41': 345}


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


def get_patients_number(gene):
    cchanges = get_codon_changes(gene)
    return (
        len(set(cc['PID'] for cc in cchanges if cc['Rx'] == 'NNRTIs')),
        len(set(cc['PID'] for cc in cchanges if cc['Rx'] == 'PIs')),
    )


def summarize_na_diffs(gene, rx):
    cchanges = get_codon_changes(gene)
    cchanges = [cc for cc in cchanges if cc['Rx'] == rx]
    na_length = GENE_LENGTHS[gene] * 3
    perpt = groupby(cchanges, lambda cc: cc['PID'])
    perpt = [sum(int(cc['NumNAChanges'])
                 for cc in ccs) * 100 / na_length for _, ccs in perpt]
    q25, q50, q75 = percentile(perpt, (25, 50, 75))
    return '{:.1f} ({:.1f}-{:.1f}%)'.format(q50, q25, q75), perpt


def summarize_aa_diffs(gene, rx):
    cchanges = get_codon_changes(gene)

    cchanges = [
        cc for cc in cchanges if cc['Rx'] == rx and cc['Type'] == 'non']
    aa_length = GENE_LENGTHS[gene]
    perpt = groupby(cchanges, lambda cc: cc['PID'])
    perpt = [len(list(ccs)) * 100 / aa_length for _, ccs in perpt]
    q25, q50, q75 = percentile(perpt, (25, 50, 75))
    return '{:.1f} ({:.1f}-{:.1f}%)'.format(q50, q25, q75), perpt


def summarize_diversity(gene, rx):
    cchanges = get_codon_changes(gene)
    cchanges = [
        [cc['PID']] + cc['AAs'].split('-->')
        for cc in cchanges
        if cc['Rx'] == rx and cc['Type'] == 'non']
    perpt = groupby(cchanges, lambda cc: cc[0])
    shrinks = 0
    total = 0
    for _, ccs in perpt:
        ccs = list(ccs)
        shrinks += len([cc for cc in ccs if len(cc[1]) > len(cc[2])])
        total += len(ccs)
    return '{}/{} ({:.1f}%)'.format(shrinks, total, shrinks / total * 100)


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


def pvalue(left, right):
    # perform two-sample Wilcoxon test, aka ‘Mann-Whitney’ test
    # see: https://stackoverflow.com/a/33890615/2644759
    p = mannwhitneyu(left, right).pvalue
    return '{:.2f}'.format(p)


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

        row = ['NA differences']
        all_nadiffs = []
        for rx in ('NNRTIs', 'PIs'):
            num, nadiffs = summarize_na_diffs(gene, rx)
            all_nadiffs.append(nadiffs)
            row.append(num)
        row.append(pvalue(*all_nadiffs))
        pairwise.append(row)

        row = ['AA differences']
        all_aadiffs = []
        for rx in ('NNRTIs', 'PIs'):
            num, aadiffs = summarize_aa_diffs(gene, rx)
            all_aadiffs.append(aadiffs)
            row.append(num)
        row.append(pvalue(*all_aadiffs))
        pairwise.append(row)

        row = ['dN/dS ratio']
        all_dndslist = []
        for rx in ('NNRTIs', 'PIs'):
            ratio, dndslist = summarize_dnds(gene, rx)
            all_dndslist.append(dndslist)
            row.append(ratio)
        row.append(pvalue(*all_dndslist))
        pairwise.append(row)

        if gene == 'Gag':
            row = ['dN/dS ratio (MA)']
            all_dndslist = []
            for rx in ('NNRTIs', 'PIs'):
                ratio, dndslist = summarize_dnds(gene, rx, 'MA')
                all_dndslist.append(dndslist)
                row.append(ratio)
            row.append(pvalue(*all_dndslist))
            pairwise.append(row)

            row = ['dN/dS ratio (C terminal)']
            all_dndslist = []
            for rx in ('NNRTIs', 'PIs'):
                ratio, dndslist = summarize_dnds(gene, rx, 'CTerminal')
                all_dndslist.append(dndslist)
                row.append(ratio)
            row.append(pvalue(*all_dndslist))
            pairwise.append(row)

        else:
            row = ['dN/dS ratio (CD)']
            all_dndslist = []
            for rx in ('NNRTIs', 'PIs'):
                ratio, dndslist = summarize_dnds(gene, rx, 'CD')
                all_dndslist.append(dndslist)
                row.append(ratio)
            row.append(pvalue(*all_dndslist))
            pairwise.append(row)

        row = ['Diversity']
        for rx in ('NNRTIs', 'PIs'):
            row.append(summarize_diversity(gene, rx))
        row.append('-')
        pairwise.append(row)

        numpt_nnrtis, numpt_pis = get_patients_number(gene)
        print(tabulate(
            pairwise, [
                'Title', 'NNRTIs ({})'.format(numpt_nnrtis),
                'PIs ({})'.format(numpt_pis), 'P value'
            ], tablefmt='pipe'))
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
