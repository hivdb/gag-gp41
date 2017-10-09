#! /usr/bin/env python

from collections import OrderedDict

from codonutils import translate_codon
from data_reader import (
    CONSENSUS, any_fasta_reader)

MUTINS = 0
MUTDEL = 1
FSINS = 2
FSDEL = 3


def build_mutation(pos, cons, codon, indel_type=-1, indel=None):
    codon = codon if indel_type != MUTDEL else ''
    aa = translate_codon(codon) if codon else ''
    if aa == cons:
        return
    return {
        'Position': pos,
        'ReferenceText': cons,
        'CodonText': codon,
        'AminoAcidText': aa,
        'InsertedCodonsText': indel if indel_type == MUTINS else '',
        'IsInsertion': indel_type == MUTINS,
        'IsDeletion': indel_type == MUTDEL,
        'IsPartial': '-' in codon
    }


def lanl_reader(gene, filename):
    consensus = CONSENSUS[gene]['AASeq']
    sequences = list(any_fasta_reader(filename))

    patients = OrderedDict()
    for header, naseq in sequences:
        accession, ptid, subtype = header.split('.', 2)
        if accession == 'K03455':
            pass
        elif not ptid or ptid == '-':
            continue
        elif ptid in patients:
            continue
        patients[ptid] = (accession, naseq, subtype)
    headers = [header for header, _, _ in patients.values()]
    sequences = [seq.upper() for _, seq, _ in patients.values()]
    subtypes = [subtype for _, _, subtype in patients.values()]
    hxb2idx = headers.index('K03455')

    insertions = set()
    for i, nas in enumerate(zip(*sequences)):
        # use HXB2 to identify indels
        is_insertion = nas[hxb2idx] == '-'
        if is_insertion:
            insertions.add(i)

    header = ['Accession']
    results = []
    for j, sequence in enumerate(sequences):
        all_codons = []
        all_indels = {}
        cur_codon = ''
        cur_ins = ''
        cur_pos = 0
        for i, na in enumerate(sequence):
            if len(cur_codon) == 3:
                cur_pos += 1
                numdel = cur_codon.count('-')
                if numdel == 3:
                    all_indels[(cur_pos, MUTDEL)] = True
                elif numdel > 0:
                    all_indels[(cur_pos, FSDEL)] = numdel
                all_codons.append(cur_codon)
                cur_codon = ''
            if i in insertions:
                if na != '-':
                    cur_ins += na
            else:
                if cur_ins:
                    mut_ins = cur_ins[:len(cur_ins) // 3 * 3]
                    fs_ins = cur_ins[-len(cur_ins) % 3:]
                    if mut_ins:
                        all_indels[(cur_pos, MUTINS)] = mut_ins
                    if fs_ins:
                        all_indels[(cur_pos, FSINS)] = fs_ins
                    cur_ins = ''
                cur_codon += na
        else:
            if len(cur_codon) != 3:
                raise ValueError(
                    'Sequence {} is not fully aligned'.format(headers[j]))
            cur_pos += 1
            numdel = cur_codon.count('-')
            if numdel == 3:
                all_indels[(cur_pos, MUTDEL)] = True
            elif numdel > 0:
                all_indels[(cur_pos, FSDEL)] = numdel
            all_codons.append(cur_codon)
        seqresult = {
            'Accession': headers[j],
            'FirstAA': 1,
            'LastAA': len(consensus),
            'Subtype': subtypes[j],
            'Mutations': [],
            'NumFrameShifts': len([p for p, t in all_indels.keys()
                                   if t in (FSINS, FSDEL)]),
            'NASequence': sequence.replace('-', ''),
            'AlignedNASequence': ''.join(all_codons)
        }
        for posm1, codon in enumerate(all_codons):
            pos = posm1 + 1
            cons = consensus[posm1]
            mut = None
            if (pos, MUTINS) in all_indels:
                # ins
                mut = build_mutation(pos, cons, codon, MUTINS,
                                     all_indels[(pos, MUTINS)])
            elif (pos, MUTDEL) in all_indels:
                # del
                mut = build_mutation(pos, cons, codon, MUTDEL,
                                     all_indels[(pos, MUTDEL)])
            else:
                mut = build_mutation(pos, cons, codon)
            if mut:
                seqresult['Mutations'].append(mut)
        results.append(seqresult)
    return results


if __name__ == '__main__':
    lanl_reader('Gag', '/app/local/hiv-db_gag_squeeze.fasta')
