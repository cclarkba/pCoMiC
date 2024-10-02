####### CONVERTING FUNCTIONS #######

from globals import CODONS

def fasta_to_dict(fasta):
    """
    Converts fasta file to dictionary where key is sequence info following '>' and value is sequence
    """
    seqs = {}
    with open(fasta) as f:
        curr_seq = ""
        curr_id = ""
        line = f.readline()
        while line:
            if line[0] == ">":
                if curr_seq != "":
                    seqs[curr_id] = curr_seq
                    curr_seq = ""
                curr_id = "".join([x for x in line if x!='\n' and x!='>'])
            else:
                curr_seq += "".join([x for x in line if x!='\n'])
            line = f.readline()
        seqs[curr_id] = curr_seq
    return seqs

def dna_to_protein(s): 
    out = ""
    for i in range(0,len(s) - 3,3):
        if CODONS[s[i:i+3]] == 'Stop':
            break
        out += CODONS[s[i:i+3]]
    return out