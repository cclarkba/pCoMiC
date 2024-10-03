####### TESTING #######

import pandas as pd
from conversions import fasta_to_dict


def make_mut_list(db):
    df = pd.read_excel(db)
    out = []
    i = 0
    while i < len(df['Variant cDNA name']):
        out.append(df['Variant cDNA name'][i])
        i += 1
    return out

# If you want to run this, change these to whatever file path each is located in
control = "C:/Users/musta/Desktop/4A/Biol469/Final/CFTR_mrna.fasta"
db = "C:/Users/musta/Desktop/4A/Biol469/Final/database.xlsx"
out = "C:/Users/musta/Desktop/4A/Biol469/Final/"

seqs = fasta_to_dict(control)
ctrl = seqs['NM_000492.4 Homo sapiens CF transmembrane conductance regulator (CFTR), mRNA']

# Testing database

# mut_list = make_mut_list(db)
# for elem in mut_list:
#     mut = induce_muts(elem,ctrl,1)
#     find_muts(mut,ctrl,db,out)
# mut = induce_muts('c.1519_1521del',ctrl,1)
# find_muts(mut,ctrl,db,out)