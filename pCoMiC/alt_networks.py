####### ALTERNATIVE NEURAL NETWORKS TO TRY IF TIME #######

import torch.nn as nn
import torch.nn.functional as F
import torch
import numpy as np
from data_loader import MutDataset, MutDataset_Unsuper


class AASeq(nn.Module):
    # A 1480x20 matrix of one-hot vectors denoting which amino acid is at each position
    # May change this to be nucleotide sequence to encompass class V mutations as well, will have to see
    # Looks like nucleotide sequence is way too large, will have to discuss this in the paper
    def __init__(self):
        super(AASeq, self).__init__()

    def forward(self, x):
        return x
    
class SeqLength(nn.Module):
    # A single scalar value representing the length of the sequence
    def __init__(self):
        super(SeqLength, self).__init__()
    
    def forward(self, x):
        return x
    
class Domain(nn.Module):
    # An Rx1 one-hot vector where R represents the number of domains identified by separate GenBank-like file
    # Whatever value is 1 in the vector denotes which domain the mutation occurs in
    def __init__(self):
        super(Domain, self).__init__()

    def forward(self, x):
        return x

class Properties(nn.Module):
    # A 3x1 vector of scalar values representing changes in hydropathy, relative charge and hydrogen bonding capability
    # between the wild-type and mutant proteins
    def __init__(self):
        super(Properties, self).__init__()

    def forward(self, x):
        return x
    
class Network(nn.Module):
    def __init__(self):
        super(Network, self).__init__()
        self.AASeq = AASeq()
        self.SeqLength = SeqLength()
        self.Domain = Domain()
        self.Properties = Properties()
    
    def forward(self, x1, x2, x3, x4):
        x1_proc = self.AASeq(x1)
        x2_proc = self.SeqLength(x2)
        x3_proc = self.Domain(x3)
        x4_proc = self.Properties(x4)
        
        x_combined = torch.cat([x1_proc, x2_proc, x3_proc, x4_proc], dim=1)
        return x_combined