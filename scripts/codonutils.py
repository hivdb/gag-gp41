CODON_TABLE = {
    'TTT': 'F',
    'TTC': 'F',
    'TTA': 'L',
    'TTG': 'L',

    'CTT': 'L',
    'CTC': 'L',
    'CTA': 'L',
    'CTG': 'L',

    'ATT': 'I',
    'ATC': 'I',
    'ATA': 'I',
    'ATG': 'M',

    'GTT': 'V',
    'GTC': 'V',
    'GTA': 'V',
    'GTG': 'V',

    'TCT': 'S',
    'TCC': 'S',
    'TCA': 'S',
    'TCG': 'S',

    'CCT': 'P',
    'CCC': 'P',
    'CCA': 'P',
    'CCG': 'P',

    'ACT': 'T',
    'ACC': 'T',
    'ACA': 'T',
    'ACG': 'T',

    'GCT': 'A',
    'GCC': 'A',
    'GCA': 'A',
    'GCG': 'A',

    'TAT': 'Y',
    'TAC': 'Y',

    'CAT': 'H',
    'CAC': 'H',
    'CAA': 'Q',
    'CAG': 'Q',

    'AAT': 'N',
    'AAC': 'N',
    'AAA': 'K',
    'AAG': 'K',

    'GAT': 'D',
    'GAC': 'D',
    'GAA': 'E',
    'GAG': 'E',

    'TGT': 'C',
    'TGC': 'C',
    'TGG': 'W',

    'CGT': 'R',
    'CGC': 'R',
    'CGA': 'R',
    'CGG': 'R',

    'AGT': 'S',
    'AGC': 'S',
    'AGA': 'R',
    'AGG': 'R',

    'GGT': 'G',
    'GGC': 'G',
    'GGA': 'G',
    'GGG': 'G',

    'TAA': '*',
    'TGA': '*',
    'TAG': '*',
}

AMBIBUOUS_NAS = {
    'W': 'AT',
    'S': 'CG',
    'M': 'AC',
    'K': 'GT',
    'R': 'AG',
    'Y': 'CT',
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ACT',
    'V': 'ACG',
    'N': 'ACGT'
}


def translate_codon(nas):
    if nas in CODON_TABLE:
        return CODON_TABLE[nas]
    nas = nas.replace('-', 'N')[:3]
    aas = set()
    for na0 in AMBIBUOUS_NAS.get(nas[0], nas[0]):
        for na1 in AMBIBUOUS_NAS.get(nas[1], nas[1]):
            for na2 in AMBIBUOUS_NAS.get(nas[2], nas[2]):
                aas.add(CODON_TABLE[na0 + na1 + na2])
    CODON_TABLE[nas] = aas = ''.join(sorted(aas))
    return aas
