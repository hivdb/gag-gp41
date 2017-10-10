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
    deltas = []
    for pid, _, prev_seq, post_seq in iter_sequence_pairs(gene, rx):
        baseline = sum([prev_seq.na_sequence.count(na)
                        for na in 'ACTG-.']) / len(prev_seq.na_sequence)
        followup = sum([post_seq.na_sequence.count(na)
                        for na in 'ACTG-.']) / len(post_seq.na_sequence)
        delta = baseline - followup
        baselines.append(100 - baseline * 100)
        followups.append(100 - followup * 100)
        deltas.append(delta * 100)
    return (
        ('{:.1f} ({:.1f}-{:.1f}%)'
         .format(*percentile(baselines, (50, 25, 75))),
         baselines),
        ('{:.1f} ({:.1f}-{:.1f}%)'
         .format(*percentile(followups, (50, 25, 75))),
         followups),
        ('{:.1f} ({:.1f}-{:.1f}%)'
         .format(*percentile(deltas, (50, 25, 75))),
         deltas)
    )


def summarize_gag_cleavage_sites():

    sites = (
        132, 373, 374, 375, 376, 378, 380, 381,
        429, 436, 451, 453, 484, 485, 490
    )
    result = []
    with open(os.path.join(
            APPDIR, 'resultData', 'aaChangesByPosWPrev',
            'gag.csv'.format(gene))
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


def insertions_table():
    rows = []
    with open(os.path.join(APPDIR, 'data', 'insertions.csv')) as fp, \
            open(os.path.join(APPDIR, 'data', 'accessions.csv')) as fp2:
        insertions = csv.DictReader(fp)
        accessions = csv.DictReader(fp2)
        accessions = {
            (int(a['PID']), a['Gene'], a['Timepoint']): a['Accession']
            for a in accessions
        }
        for ins in insertions:
            ptid = int(ins['PID'])
            key = (ptid, ins['Gene'], ins['Timepoint'])
            if key not in accessions:
                continue
            accession = accessions[key]
            rows.append([
                accession, ins['Gene'], ptid, ins['Timepoint'],
                int(ins['Pos']), ins['AA'], ins['NA']
            ])
    return rows


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
        delta_row = ['Delta']

        for rx in ('PIs', 'NNRTIs'):
            ((numbl, baselines),
             (numfu, followups),
             (numdt, deltas)) = summarize_ambiguities(gene, rx)
            all_baselines.append(baselines)
            all_followups.append(followups)
            all_deltas.append(deltas)
            baseline_row.append(numbl)
            followup_row.append(numfu)
            delta_row.append(numdt)
        baseline_row.append(pvalue(*all_baselines))
        followup_row.append(pvalue(*all_followups))
        delta_row.append(pvalue(*all_deltas))

        pairwise.extend([baseline_row, followup_row, delta_row])

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
        ['Rx', 'Position', 'AA change', '# patients'], tablefmt='pipe'))

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
    print('\n\n## Insertion table\n')
    insertions = insertions_table()
    print(tabulate(
        insertions,
        ['Accession', 'Gene', 'PID', 'Timepoint',
         'Position', 'AminoAcid', 'NucleicAcid'],
        tablefmt='pipe'
    ))
    print()
