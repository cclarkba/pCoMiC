####### GLOBAL VARIABLES #######

## TODO Convert the domain information into a GenBank-like file, want to avoid having these gloval variables

CODONS = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V',
        'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V',
        'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V',
        'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V',
        'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',
        'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
        'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
        'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
        'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D',
        'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
        'TAA': 'Stop', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
        'TAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
        'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G',
        'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
        'TGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
        'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}

PROPERTIES = {
    # 'ONE-LETTER CODON': (PI, HYDROPATHY SCORE, H BONDS)
    'A': (6.11,1.8,1), 'R': (10.76,-4.5,2), 'N': (5.41,-3.5,2), 'D': (2.87,-3.5,2), 'C': (5.02,2.5,1),
    'Q': (5.65,-3.5,2), 'E': (3.08,-3.5,2), 'G': (6.06,-0.4,1), 'H': (7.64,-3.2,1), 'I': (6.04,4.5,1),
    'L': (6.04,3.8,1), 'K': (9.47,-3.9,2), 'M': (5.74,1.9,1), 'F': (5.91,2.8,1), 'P': (6.30,-1.6,0), 
    'S': (5.68,-0.8,2), 'T': (5.60,-0.7,2), 'W': (5.88,-0.9,1), 'Y': (5.63,-1.3,2), 'V': (6.02,4.2,1)
}

# Nucleotide Sequence
CDS = (70,4513)

# Protein Sequence
NBD1 = (372,683)
NBD2 = (1206,1480)
TM1 = (77,98)
TM6 = (339,358)
TM12 = (1130,1152)
REG = (653,831)


