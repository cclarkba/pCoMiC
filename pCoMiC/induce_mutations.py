####### INDUCE MUTATIONS #######

def induce_muts(mutation, reference):
    """if not re.match('c.+[0-9]_+[0-9](del|dup|ins+[A,C,T,G]|delins+[A,C,T,G])|c.+[0-9](dup|[A,C,T,G]>[A,C,T,G])',mutation):
        print("Please enter correct formatting for mutation. \nRegEx: 'c.+[0-9]_+[0-9](del|dup|ins+[A,C,T,G]|delins+[A,C,T,G])|c.+[0-9](dup|[A,C,T,G]>[A,C,T,G])'")
        sys.exit()"""
    # seq = fasta_to_dict(reference)
    seq = reference
    posn = ''
    mut = ''
    second_int = ''
    multiple_posn = False
    for c in mutation:
        if c.isdigit():
            if multiple_posn:
                second_int += c
            else:
                posn += c
        elif c == '_':
            multiple_posn = True
        elif c in ('c','.'):
            continue
        else:
            mut += c
    insertion = 'ins' in mut
    deletion = 'del' in mut
    if '>' in mut:
        # substitution
        seq = seq[:int(posn)-1] + mut[-1] + seq[int(posn):]
    elif insertion and not deletion:
        # insertion
        seq = seq[:int(posn)] + mut[3:] + seq[int(posn):]
    elif deletion and not insertion:
        # deletion
        if second_int == '':
            seq = seq[:int(posn)-1] + seq[int(posn):]
        else:
            seq = seq[:int(posn)-1] + seq[int(second_int):]
    elif 'dup' in mut:
        # duplication
        if second_int == '':
            seq = seq[:int(posn)] + seq[int(posn)-1] + seq[int(posn):]
        else:
            seq = seq[:int(second_int)] + seq[int(posn)-1:int(second_int)] + seq[int(second_int):]
    else:
        #delins
        if second_int == '':
            seq = seq[:int(posn)-1] + mut[6:] + seq[int(posn):]
        else:
            seq = seq[:int(posn)-1] + mut[6:] + seq[int(second_int):]
    return seq

if __name__ == "__main__":
    from conversions import fasta_to_dict
    control = "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/CFTR_mrna.fasta"
    db = "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/database.xlsx"
    out = "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/output"

    seqs = fasta_to_dict(control)
    ctrl = seqs['NM_000492.4 Homo sapiens CF transmembrane conductance regulator (CFTR), mRNA']
    mut = induce_muts('c.1519_1521del', ctrl, out)