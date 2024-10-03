####### SCORING #######

import math
from globals import NBD1,NBD2,TM1,TM12,TM6,PROPERTIES,REG
from align import align
from conversions import dna_to_protein 

def sig_calc(score,weight,shift):
    """"
    Sigmoidal normalizing function
    """
    # Previous shift was 0.35, higher makes smaller values more insignificant, shifts midpoint
    # Previous weight was 2, higher also makes smaller values more insignificant, less impact on steepness
    if score > 1:
        score = 1
    if score < 0:
        score = 0
    return 1/(1+((score**(math.log(weight)/math.log(shift))-1))**2)

def proc_score(pos,hydro,charge,hbond):
    """
    Calculates probability of categorizing mutation into processing category
    """
    # Will weigh everything equally, if values are below certain threshold, will add to see their combined effects
    # Else, will only take highest weighted factor since others are likely insignificant
    if pos in range(NBD1[0],NBD1[1]) or pos in range(NBD2[0],NBD2[1]):
            # Reweighing for important region, hydro matters most then charge then hbond
            hydro = sig_calc(hydro,2,0.25)
            charge = sig_calc(charge,2,0.45)
            if hbond != 0:
                hbond = sig_calc(charge,2,0.5)
    if pos == 507:
        hydro = sig_calc(hydro,1.5,0.1)
    if hydro+charge+hbond < 0.6:
        return sig_calc(hydro+charge+hbond,2,0.3)
    else:
        return sig_calc(max(hydro,charge,hbond),2.5,0.3)

def gat_score(pos,hydro,charge,ref_aligned,mut_aligned):
    """
    Calculates probability of categorizing mutation into gating category
    """
    # Weigh NBD regions and specific residues higher (T,Y,S), weigh hydro higher than charge
    charge = sig_calc(charge,2,0.7)
    weight = 2
    shift = 0.25
    if pos in range(NBD1[0],NBD1[1]) or pos in range(NBD2[0],NBD2[1]):
        # Mutation in important NBD1 domain, critical for gating
        shift = 0.3
    if ref_aligned[pos] in ('Y','T','S') and mut_aligned not in ('Y','T','S'):
        # Mutated residue may have been involved in phosphorylation before but now is not
        weight = 1.65
    if hydro+charge<0.5:
        return sig_calc(hydro+charge,weight,shift)
    else:
        if pos == 507:
            shift = 0.6
        return sig_calc(max(hydro,charge),weight+0.5,shift+0.1)

def cond_score(pos,hydro,charge,hbond,ref_aligned,mut_aligned):
    """
    Calculates probability of categorizing mutation into conducting category
    """
    # Weigh charge at transmembrane regions highest, hbond and hydro will be weighed to account for tightening or loosening of channel disallowing only chloride from entering
    # Really only need to take into account channel regions
    weight = 2.5
    shift = 0.55
    if pos in range(TM1[0],TM1[1]) or pos in range(TM6[0],TM6[1]) or pos in range(TM12[0],TM12[1]):
        # Mutation in transmembrane region comprising channel, increase weight
        shift = 0.25
        weight = 2
    try: mut_charge = (PROPERTIES[mut_aligned[pos]][0]/7.4) - 1
    except: mut_charge = 0
    try: ref_charge = (PROPERTIES[ref_aligned[pos]][0]/7.4) - 1
    except: ref_charge = 0
    if mut_charge < 0 and ref_charge >= 0:
        # Mutated residue has negative relative charge, original residue has positive relative charge, will affect conductance
        return sig_calc(charge,weight-0.5,shift-0.1)
    elif hydro+charge+hbond < 0.6:
        return sig_calc(hydro+charge+hbond,weight,shift)
    else:
        # Hydro and hbond insignificant, only account for most important factor being charge
        return sig_calc(charge,weight,shift)

def insuf_score(pos,hydro,hbond,charge):
    """
    Calculates probability of categorizing mutation into insufficient category
    """
    # Can't do much here, weigh everything pretty low since they don't really contribute to this, maybe score regulatory domain higher, look into this
    weight = 4
    shift = 0.75
    if pos in range(REG[0],REG[1]):
        shift = 0.35
    if all(i<0.25 for i in (hydro,charge,hbond)):
        return sig_calc(hydro+charge+hbond,weight,shift)
    else:
        return sig_calc(max(hydro,charge,hbond),weight,shift)

def predict(ref, mut):
    """"
    Calculates predicted probabilities for a mutation to be classed into each category individually based on changes in relative charge, hydropathy and hydrogen bonding capability from WT
    """
    ref_prot = dna_to_protein(ref)
    mut_prot = dna_to_protein(mut)
    scores = {
        'production': 0,
        'processing': 0,
        'gating': 0,
        'conducting': 0,
        'insufficient': 0
    }
    alignment = align(ref_prot,mut_prot)
    ref_aligned = alignment[0]
    mut_aligned = alignment[1]
    diff = len(ref_prot) - len(mut_prot)
    if diff > 15:
        scores['production'] = 1
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
                scores['processing'] += proc_score(i,hydro_change,charge_change,lost_hbonds_score)
                scores['gating'] += gat_score(i,hydro_change,charge_change,ref_aligned,mut_aligned)
                scores['conducting'] += cond_score(i,hydro_change,charge_change,lost_hbonds_score,ref_aligned,mut_aligned)
                scores['insufficient'] += insuf_score(i,hydro_change,charge_change,lost_hbonds_score)
            i += 1
    return scores
