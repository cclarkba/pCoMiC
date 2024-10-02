####### MAIN FUNCTIONS #######
from align import align
import sys
from cleaning import clean_muts
from output import cleaned_output, scores_figure
from scoring import predict

def find_muts(mut_seq,ref_seq,db,out):
    """
    Identifies all single basepair mutations in mut_seq compared to ref_seq and labels as either insertion, deletion or substitution event
    """
    try:
        alignment = align(ref_seq,mut_seq)
    except:
        print("Error: Alignment likely timed out. Ensure sequences correspond to CFTR CDS and device has sufficient memory and computing power.")
        sys.exit()
    ref_aligned = alignment[0]
    mut_aligned = alignment[1]
    i = 0
    mutations = {}
    while i < len(ref_aligned):
        if ref_aligned[i] != mut_aligned[i]:
            if ref_aligned[i] == '-':
                # insertion
                mutation = f'ins{mut_seq[i]}'
            elif mut_aligned[i] == '-':
                # deletion
                mutation = f'del{ref_seq[i]}'
            else:
                # substituion
                mutation = f'{ref_seq[i]}>{mut_seq[i]}'
            mutations[i+1] = mutation
        i += 1
    cleaned = clean_muts(mutations,ref_aligned,mut_aligned)

    # Call for output of cleaned muts
    scores = predict(ref_seq,mut_seq)
    cleaned_output(cleaned,db,out,scores)
    # Call for creation of figure
    scores_figure(cleaned,db,out,scores)
