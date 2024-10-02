####### ALIGNER #######

from Bio import Align

def align(ref,mut):
    """
    Simple pairwise alignment between reference seq and mutant seq
    """
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    aligner.mismatch_score = -1
    aligner.match_score = 2
    alignments = aligner.align(ref,mut)
    alignment = alignments[0]
    return alignment
