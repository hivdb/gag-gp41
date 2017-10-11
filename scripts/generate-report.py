import os
import csv
import math
from itertools import groupby

from tabulate import tabulate
from numpy import percentile
from scipy.stats import mannwhitneyu
from analysis_functions import iter_sequence_pairs

APPDIR = '/app'
FASTAS = os.path.join(APPDIR, 'data/fasta')
HYPHYOUT = os.path.join(APPDIR, 'result_data/hyphy_output')
CLEANOUT = os.path.join(APPDIR, 'result_data/hyphy_cleaned_output')
NAN = float('nan')
GENE_LENGTHS = {'gag': 500, 'gp41': 345}
DOMAINS = (
    ('All', None, None),
    ('MA', (1, 132), 'gag'),
    ('CTerminal', (364, 500), 'gag'),
    ('CD', (195, 345), 'gp41'))

GAG_CLEAVAGE_SITES_DEFINE = [
    ('MA/CA', (132, 133)),   # 132/133+-5
    ('CA/SP1', (363, 364)),  # 363/364+-5
    ('SP1/NC', (377, 378)),  # 377/378+-5
    ('NC/SP2', (432, 433)),  # 432/433+-5
    ('SP2/p6', (448, 449)),  # 448/449+-5
    ('p6/PR', (488, 489))    # 488/489+-5
]


def summarize_dnds(gene, rx, domain=None):
    if domain == 'All':
        domain = ''
    elif domain:
        domain = '-' + domain
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


def summarize_na_diffs(gene, rx, domain_range=None):
    cchanges = get_codon_changes(gene, domain_range)
    cchanges = [cc for cc in cchanges if cc['Rx'] == rx]
    na_length = GENE_LENGTHS[gene] * 3
    perpt = groupby(cchanges, lambda cc: cc['PID'])
    perpt = [sum(int(cc['NumNAChanges'])
                 for cc in ccs) * 100 / na_length for _, ccs in perpt]
    q25, q50, q75 = percentile(perpt, (25, 50, 75))
    return '{:.1f} ({:.1f}-{:.1f}%)'.format(q50, q25, q75), perpt


def summarize_aa_diffs(gene, rx, domain_range=None):
    cchanges = get_codon_changes(gene, domain_range)
    cchanges = [
        cc for cc in cchanges if cc['Rx'] == rx and cc['Type'] == 'non']

    aa_length = GENE_LENGTHS[gene]
    perpt = groupby(cchanges, lambda cc: cc['PID'])
    perpt = [len(list(ccs)) * 100 / aa_length for _, ccs in perpt]
    q25, q50, q75 = percentile(perpt, (25, 50, 75))
    return '{:.1f} ({:.1f}-{:.1f}%)'.format(q50, q25, q75), perpt


def summarize_ambiguities(gene, rx):
    baselines = []
    followups = []
    for pid, _, prev_seq, post_seq in iter_sequence_pairs(gene, rx):
        baseline = sum([prev_seq.na_sequence.count(na)
                        for na in 'ACTG-.']) / len(prev_seq.na_sequence)
        followup = sum([post_seq.na_sequence.count(na)
                        for na in 'ACTG-.']) / len(post_seq.na_sequence)
        baselines.append(100 - baseline * 100)
        followups.append(100 - followup * 100)
    return (
        ('{:.1f} ({:.1f}-{:.1f}%)'
         .format(*percentile(baselines, (50, 25, 75))),
         baselines),
        ('{:.1f} ({:.1f}-{:.1f}%)'
         .format(*percentile(followups, (50, 25, 75))),
         followups)
    )


def summarize_gag_cleavage_sites():
    result = []
    for cat in ('PIs', 'NNRTIs'):
        for pid, cat, prev_seq, post_seq in iter_sequence_pairs('gag', cat):
            for name, (start, end) in GAG_CLEAVAGE_SITES_DEFINE:
                start -= 4
                end += 4
                prev_codons = list(prev_seq.iter_codons(start, end))
                post_codons = list(post_seq.iter_codons(start, end))
                if not prev_codons or not post_codons:
                    continue
                pos = []
                prev_aas = []
                post_aas = []
                for i, (codon0, codon1) in \
                        enumerate(zip(prev_codons, post_codons)):
                    aa0, aa1 = [codon0.aa, codon1.aa]
                    if len(aa0) > 1:
                        aa0 = aa0.replace(codon0.cons_aa, '')
                    if len(aa1) > 1:
                        aa1 = aa1.replace(codon1.cons_aa, '')
                    prev_aas.append(aa0 if len(aa0) == 1 else 'X')
                    post_aas.append(aa1 if len(aa1) == 1 else 'X')
                    if aa0 != aa1:
                        pos.append(str(start + i))
                prev_aas = ''.join(prev_aas)
                post_aas = ''.join(post_aas)
                if prev_aas == post_aas:
                    continue
                result.append([cat, name, ', '.join(pos),
                               prev_aas, post_aas])

    def sortkey(r):
        return (r[0], r[1], r[3], r[4])

    result = sorted(result, key=sortkey)
    result = groupby(result, sortkey)

    agg_result = []
    for _, rows in result:
        rows = list(rows)
        agg_result.append(rows[0] + [len(rows)])

    return sorted(
        agg_result, key=lambda r: ({'PIs': 0, 'NNRTIs': 1}[r[0]],
                                   int(r[2].split(', ', 1)[0]),
                                   r[5], r[3], r[4]))


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


def get_codon_changes(gene, domain_range=None):
    with open(os.path.join(
            APPDIR, 'resultData', 'codonChangesByPt',
            '{}.csv'.format(gene))
    ) as fp:
        reader = csv.DictReader(fp)
        cchanges = list(reader)
    if domain_range:
        left, right = domain_range
        cchanges = [cc for cc in cchanges if left <= int(cc['Pos']) <= right]
    return cchanges


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
    print(
        '### Distribution figures for Naive Sequences\n\n'
        '![Gag AA changes](https://github.com/hivdb/gag-gp41'
        '/raw/master/report/gag-naive-aachanges-dist.png)\n'
        '![Gag stop codons](https://github.com/hivdb/gag-gp41'
        '/raw/master/report/gag-naive-stopcodons-dist.png)\n'
        '![Gp41 AA changes](https://github.com/hivdb/gag-gp41'
        '/raw/master/report/gp41-naive-aachanges-dist.png)\n'
        '![Gp41 stop codons](https://github.com/hivdb/gag-gp41'
        '/raw/master/report/gp41-naive-stopcodons-dist.png)\n'
    )
    print('## Pairwise comparison\n')
    for gene in ('gag', 'gp41'):
        pairwise = []
        print('### {}\n'.format(gene))

        pairwise.append(['NA differences', '', '', ''])

        for domain, domain_range, geneOnly in DOMAINS:
            if geneOnly and geneOnly != gene:
                continue
            row = ['- {}'.format(domain)]
            all_nadiffs = []
            for rx in ('PIs', 'NNRTIs'):
                num, nadiffs = summarize_na_diffs(gene, rx, domain_range)
                all_nadiffs.append(nadiffs)
                row.append(num)
            row.append(pvalue(*all_nadiffs))
            pairwise.append(row)

        pairwise.append(['AA differences', '', '', ''])

        for domain, domain_range, geneOnly in DOMAINS:
            if geneOnly and geneOnly != gene:
                continue
            row = ['- {}'.format(domain)]
            all_aadiffs = []
            for rx in ('PIs', 'NNRTIs'):
                num, aadiffs = summarize_aa_diffs(gene, rx, domain_range)
                all_aadiffs.append(aadiffs)
                row.append(num)
            row.append(pvalue(*all_aadiffs))
            pairwise.append(row)

        pairwise.append(['dN/dS ratio', '', '', ''])

        for domain, domain_range, geneOnly in DOMAINS:
            if geneOnly and geneOnly != gene:
                continue
            row = ['- {}'.format(domain)]
            all_dndslist = []
            for rx in ('PIs', 'NNRTIs'):
                ratio, dndslist = summarize_dnds(gene, rx, domain)
                all_dndslist.append(dndslist)
                row.append(ratio)
            row.append(pvalue(*all_dndslist))
            pairwise.append(row)

        pairwise.append(['% Ambiguities', '', '', ''])

        all_baselines = []
        all_followups = []
        all_deltas = []

        baseline_row = ['Baseline']
        followup_row = ['Follow-up']
        pvalue_row = ['P value']

        for rx in ('PIs', 'NNRTIs'):
            ((numbl, baselines),
             (numfu, followups)) = summarize_ambiguities(gene, rx)
            all_baselines.append(baselines)
            all_followups.append(followups)
            baseline_row.append(numbl)
            followup_row.append(numfu)
            pvalue_row.append(pvalue(baselines, followups))
        baseline_row.append(pvalue(*all_baselines))
        followup_row.append(pvalue(*all_followups))
        pvalue_row.append('-')

        pairwise.extend([baseline_row, followup_row, pvalue_row])

        numpt_nnrtis, numpt_pis = get_patients_number(gene)
        print(tabulate(
            pairwise, [
                '', 'PIs ({})'.format(numpt_pis),
                'NNRTIs ({})'.format(numpt_nnrtis), 'P value'
            ], tablefmt='pipe'))
        print()

    print('\n## Gag cleavage sites\n')
    print(tabulate(
        summarize_gag_cleavage_sites(),
        ['Rx', 'Cleavage site', 'Position',
         'Baseline AAs', 'Follow-up AAs', '# patients'], tablefmt='pipe'))

    print('\n## Positions with evidence for diversifying selection (FEL)\n')
    fel = []
    for gene in ('gag', 'gp41'):
        for rx in ('PIs', 'NNRTIs'):
            fel.extend(fel_result(
                '{} - {}'.format(gene, rx),
                os.path.join(CLEANOUT, '{}{}.fel.tsv'.format(gene, rx))
            ))
    print(tabulate(fel, ['Group', 'Position', 'P value'], tablefmt='pipe'))
    print('\n\n## Positions with evidence for directional selection (MEDS)\n')
    meds = []
    for gene in ('gag', 'gp41'):
        for rx in ('PIs', 'NNRTIs'):
            meds.extend(meds_result(gene, rx))

    print(tabulate(
        meds, ['Group', 'Mutation', 'Detail', '# patients'], tablefmt='pipe'))
    # print('\n\n## Insertion table\n')
    # insertions = insertions_table()
    # print(tabulate(
    #     insertions,
    #     ['Accession', 'Gene', 'PID', 'Timepoint',
    #      'Position', 'AminoAcid', 'NucleicAcid'],
    #     tablefmt='pipe'
    # ))
    print()
