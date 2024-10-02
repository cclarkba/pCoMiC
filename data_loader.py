####### DATA LOADER #######

import pandas as pd
import numpy as np
import torch
from torch.nn.utils.rnn import pad_sequence
from induce_mutations import induce_muts
from conversions import dna_to_protein
from globals import PROPERTIES
from sklearn.model_selection import train_test_split

# class MutDataset(torch.utils.data.Dataset):
#     def __init__(self, data, targets, test_size):

#         train_data, test_data, train_labels, test_labels = train_test_split(
#             data,
#             targets,
#             test_size=test_size,
#             random_state=42
#         )
#         self.train_dataset = pad_sequence([torch.tensor(d) for d in train_data], batch_first=True, padding_value=0)
#         self.train_labels = torch.tensor(np.array(train_labels))
#         self.test_dataset = pad_sequence([torch.tensor(d) for d in test_data], batch_first=True, padding_value=0)
#         self.test_labels = torch.tensor(np.array(test_labels))
    
#     def __len__(self):
#         return len(self.train_dataset), len(self.test_dataset)
    
#     def __getitem__(self, idx, dataset):
#         if dataset == 'train':
#             return self.train_dataset[idx], self.train_labels[idx]
#         else:
#             return self.test_dataset[idx], self.test_labels[idx]
        
class MutDataset(torch.utils.data.Dataset):
    def __init__(self, data, targets):
        self.dataset = pad_sequence([torch.tensor(d) for d in data], batch_first=True, padding_value=0)
        self.labels = torch.tensor(np.array(targets))
    
    def __len__(self):
        return len(self.dataset)
    
    def __getitem__(self, idx, ):
        return self.dataset[idx], self.labels[idx]
    
class MutDataset_Unsuper(torch.utils.data.Dataset):
    def __init__(self, data):
        self.dataset = pad_sequence([torch.tensor(d) for d in data], batch_first=True, padding_value=0)
    
    def __len__(self):
        return len(self.dataset)
    
    def __getitem__(self, idx):
        return self.dataset[idx]

def excel_to_dataset(db, ref_seq, mut_cats, out, supervised=True):
    '''
    Takes in an excel database formatted as the example and outputs data and targets.
    Data is an array of mutated sequences represented as a Nx20 matrix where each row is a one-hot vector corresponding to the amino acid as position n
    Targets is a multi-hot vector corresponding to the correct mutant classification [Production, Processing, Gating, Conducting, Insufficient]
    data and targets passed through MutDataset class and returned
    '''
    df = pd.read_excel(db)
    mut_cdna = df['Variant cDNA name'].tolist()

    data = []
    aa = sorted(PROPERTIES.keys())
    if supervised:
        classes = df['Expected'].tolist()
        targets = []
        for mut,c in zip(mut_cdna, classes):
            # Apply mutation to nucleotide seq, convert to aaseq then to one-hot matrix
            mut_prot = dna_to_protein(induce_muts(mut, ref_seq))
            mat = np.zeros((len(mut_prot),len(aa)))
            for i,a in enumerate(mut_prot):
                idx = aa.index(a)
                mat[i][idx] = 1
            data.append(mat)

            # Run through targets for said mutant and update corresponding multi-hot vector
            targs = np.zeros(5)
            c = c.split(',')
            c = [x.strip() for x in c]
            for j in c:
                try:
                    targs[mut_cats.index(j)] = 1
                except:
                    pass
            targets.append(targs)
        dataset = MutDataset(data, targets)

    else:
        # unsupervised dataset with no classes
        for mut in mut_cdna:
            mut_prot = dna_to_protein(induce_muts(mut, ref_seq))
            mat = np.zeros((len(mut_prot), len(aa)))
            for i,a in enumerate(mut_prot):
                idx = aa.index(a)
                mat[i][idx] = 1
            data.append(mat)
        dataset = MutDataset_Unsuper(data)


    torch.save(dataset,out)
    return dataset

import math
from globals import NBD1,NBD2,TM1,TM12,TM6,PROPERTIES,REG
from align import align
from conversions import dna_to_protein 

def char_scores(ref, mut):
    """"
    Calculates predicted probabilities for a mutation to be classed into each category individually based on changes in relative charge, hydropathy and hydrogen bonding capability from WT
    """
    ref_prot = dna_to_protein(ref)
    mut_prot = dna_to_protein(mut)
    scores = {
        'hydro': 0,
        'charge': 0,
        'hbond': 0,
        'size': 0
    }
    alignment = align(ref_prot,mut_prot)
    scores['size'] += len(mut_prot) / len(ref_prot)
    ref_aligned = alignment[0]
    mut_aligned = alignment[1]
    if abs(len(ref_prot) - len(mut_prot)) > 15:
        return scores
    else:
        i = 0
        while i < len(ref_aligned):
            if ref_aligned[i] != mut_aligned[i]: # Mutation Identified
                # Calculating average surrounding hydropathy and charge
                j,k = i-1,i+1
                avg_hydro = 0
                ref_avg_hydro = 0
                avg_charge = 0
                ref_avg_charge = 0
                aa_seen = 0
                while j > i - 4:
                    hydro_tot = 0 
                    charge_tot = 0
                    if j >= 0:
                        val = PROPERTIES.get(mut_aligned[j])
                        if val:
                            hydro_tot += val[1]
                            charge_tot += 1 - (val[2]/7.4)
                            aa_seen += 1
                        j -= 1
                    if k < len(mut_prot):
                        val = PROPERTIES.get(mut_aligned[k])
                        if val:
                            hydro_tot += val[1]
                            charge_tot += 1 - (val[2]/7.4)
                            aa_seen += 1
                        k += 1
                    avg_hydro += hydro_tot
                    avg_charge += charge_tot

                # Defining parameters for H Bonding
                lost_hbonds = 0
                if ref_aligned[i] == '-': # Insertion
                    # Hydropathy
                    ref_avg_hydro = avg_hydro
                    avg_hydro += PROPERTIES[mut_aligned[i]][1]
                    ref_avg_hydro /= aa_seen
                    avg_hydro /= (aa_seen + 1)
                    # No H bonds lost in insertion
                    # Charge
                    ref_avg_charge = avg_charge
                    avg_charge += PROPERTIES[mut_aligned[i]][0]
                    ref_avg_charge /= aa_seen
                    avg_charge /= (aa_seen + 1)
                elif mut_aligned[i] == '-': # Deletion
                    # Hydropathy
                    ref_avg_hydro = avg_hydro + PROPERTIES[ref_aligned[i]][1]
                    ref_avg_hydro /= (aa_seen + 1)
                    avg_hydro /= aa_seen
                    # H Bonding
                    lost_hbonds = PROPERTIES[ref_aligned[i]][2]
                    # Charge
                    ref_avg_charge = avg_charge + PROPERTIES[ref_aligned[i]][0]
                    ref_avg_charge /= (aa_seen + 1)
                    avg_charge /= aa_seen
                else: # Substitution
                    # Hydropathy
                    aa_seen += 1
                    ref_avg_hydro = (avg_hydro + PROPERTIES[ref_aligned[i]][1]) / aa_seen
                    if ref_avg_hydro == 0:
                        ref_avg_hydro = avg_hydro/(aa_seen-1)
                    avg_hydro = (avg_hydro + PROPERTIES[mut_aligned[i]][1]) / aa_seen
                    # H Bonding
                    lost_hbonds = PROPERTIES[ref_aligned[i]][2] - PROPERTIES[mut_aligned[i]][2]
                    if lost_hbonds < 0: lost_hbonds = 0
                    # Charge
                    ref_avg_charge = (avg_charge + PROPERTIES[ref_aligned[i]][0]) / aa_seen
                    avg_charge = (avg_charge + PROPERTIES[mut_aligned[i]][0]) / aa_seen

                hydro_change = abs((ref_avg_hydro-avg_hydro)/ref_avg_hydro)
                # Hydropathy being weighed too heavily, need to normalize this by running through sigmoidal normalizer
                scale = 15
                if hydro_change > 1:
                    if hydro_change > 3:
                        scale = 85
                    elif hydro_change > 2:
                        scale = 55
                    ### TODO gotta change this additional sigmoidal normalizing function to its own function
                    hydro_change = 5/(5+(scale*math.e**-hydro_change))
                charge_change = abs((ref_avg_charge-avg_charge)/ref_avg_charge)
                lost_hbonds_score = 0
                if lost_hbonds != 0:
                    lost_hbonds_score = 1/(1+(5*math.e**(lost_hbonds)))

                # Calling scoring functions and adding to dict
                # Will probably have to add additional thing for when there are more than one mutations
                scores['hydro'] += hydro_change
                scores['charge'] += charge_change
                scores['hbond'] += lost_hbonds_score
            i += 1
    return scores

def char_dataset(db, ref_seq, mut_cats, out, supervised=True):
    '''
    Takes in an excel database formatted as the example and outputs data and targets.
    Uses only 4 inputs, ratio of normal/mutated seq length and three scores
    '''
    df = pd.read_excel(db)
    mut_cdna = df['Variant cDNA name'].tolist()

    data = []
    if supervised:
        classes = df['Expected'].tolist()
        targets = []
        for mut,c in zip(mut_cdna, classes):
            # Apply mutation to nucleotide seq, convert to aaseq then to one-hot matrix
            mut_seq = induce_muts(mut, ref_seq)
            scores = char_scores(ref_seq, mut_seq)
            hydro = scores['hydro']
            charge = scores['charge']
            hbond = scores['hbond']
            size = scores['size']
            input = np.array([hydro, charge, hbond, size])
            data.append(input)

            # Run through targets for said mutant and update corresponding multi-hot vector
            targs = np.zeros(5)
            c = c.split(',')
            c = [x.strip() for x in c]
            for j in c:
                try:
                    targs[mut_cats.index(j)] = 1
                except:
                    pass
            targets.append(targs)
        dataset = MutDataset(data, targets)

    else:
        # unsupervised dataset with no classes
        for mut in mut_cdna:
            mut_seq = induce_muts(mut, ref_seq)
            scores = char_scores(ref_seq, mut_seq)
            hydro = scores['hydro']
            charge = scores['charge']
            hbond = scores['hbond']
            size = scores['size']
            input = np.array([hydro, charge, hbond, size])
            data.append(input)
        dataset = MutDataset_Unsuper(data)


    torch.save(dataset,out)
    return dataset


if __name__ == '__main__':
    from conversions import fasta_to_dict
    control = "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/CFTR_mrna.fasta"
    seqs = fasta_to_dict(control)
    ctrl = seqs['NM_000492.4 Homo sapiens CF transmembrane conductance regulator (CFTR), mRNA']
    mut_cats = ['Production', 'Processing', 'Gating', 'Conducting', 'Insufficient']

    # db = "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/database.xlsx"
    # out = "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/dataset_char_scores_super.pt"
    # dataset = char_dataset(db, ctrl, mut_cats, out, supervised=True)
    db = "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/database.xlsx"
    out = "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/dataset_split.pt"
    dataset = excel_to_dataset(db, ctrl, mut_cats, out, supervised=True)

