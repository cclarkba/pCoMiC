####### CLEANING #######

def clean_muts(mutations,ref_aligned,mut_aligned):
    """
    Takes single basepair mutation output from find_muts and cleans them. 
    Converts list of single basepair mutations to variant cDNA form for later comparison of database.
    Assumes that mutations will likely not simultaneously occur within 50bp of one another and a single mutation will not exceed 50bp.
    For total deletion and delins events, accounts for one basepair of overhang that may have naturally aligned to genome
        and reports 4 potential mutations.
    Ex. {6056: 'delT', 6057: 'delA'}, will convert to c.6056_6057del, c.6055_6057delinsX, c.6056_6058delinsX, c.6055_6058delinsXX
        {6056: 'delT', 6057: 'G>C'} will convert to c.6056_6057delinsC, c.6055_6057delinsXC, c.6056_6058delinsCX, c.6055_6058delinsXCX
        where X represents the inserted basepair that just so happened to align to the genome in the deleted region
    """

    if not mutations: return 'No mutations found.'

    keys = sorted(list(mutations.keys()))
    intervals = []
    if keys[-1] - keys[0] <= 50:
        intervals.append([keys[0],keys[-1]])
    else:
        i = 1
        start = keys[i]
        while i < len(keys):
            if keys[i] - keys[i-1] > 50:
                intervals.append([start,keys[i-1]])
                start = keys[i]
            i += 1
        intervals.append([start,keys[i-1]])
    
    cleaned = {}
    for interval in intervals:
        if interval[1] - interval[0] == 0:
            # If only single basepair mutation observed
            if (mut_aligned[interval[0]-1] == ref_aligned[interval[0]-2] or mut_aligned[interval[0]-1] == ref_aligned[interval[0]]) and \
                ref_aligned[interval[0]-1] == '-':
                # Accounting for duplication which would appear as ins from find_muts
                cleaned[interval[0]] = [f'c.{interval[0]}dup']
            else:
                cleaned[interval[0]] = [f'c.{interval[0]}{mutations[interval[0]]}']
        else:
            ref_interval = ref_aligned[interval[0]-1:interval[1]]
            mut_interval = mut_aligned[interval[0]-1:interval[1]]
            key = str(interval[0]+1) + '_' + str(interval[1]+1)
            if mut_interval == '-'*len(mut_interval):
                # deletion of entire interval
                cleaned[key] = [f'c.{key}del']
                # Need to account for delins one bp to either side
                cleaned[key].append(f'c.{key}delins{mut_aligned[interval[0]-2]}')
                cleaned[key].append(f'c.{key}delins{mut_aligned[interval[1]]}')
            elif ref_interval == '-'*len(ref_interval):
                if mut_interval == ref_aligned[interval[0]-1+len(ref_interval):interval[1]+len(ref_interval)] or mut_interval == ref_aligned[interval[0]-1-len(ref_interval):interval[1]-len(ref_interval)]:
                    # duplication of entire interval
                    cleaned[key] = [f'c.{interval[0]}_{interval[1]}dup']
                else:
                    # complete insertion event, insertion occurs between two bps 
                    cleaned[key] = [f'c.{interval[0]-1}_{interval[0]}ins{mut_interval}']
            else:
                # delins
                # Could include substitutions as well, delins behaves very oddly and will always have random binding throughout
                ins = ''
                for nuc in mut_interval:
                    if nuc != '-':
                        ins += nuc
                cleaned[key] = [f'c.{key}delins{ins}']
                # Account for ins one bp to either side that aligned naturally with ref gene
                ins_1 = mut_aligned[interval[0]-2] + ins
                ins_2 = ins + mut_aligned[interval[1]] 
                cleaned[key].append(f'c.{key}delins{ins_1}')
                cleaned[key].append(f'c.{key}delins{ins_2}')
    return cleaned