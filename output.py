####### OUTPUT ########

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def cleaned_output(cleaned,db,out,scores):
    df = pd.read_excel(db)
    known = []
    unknown = []
    for interval in cleaned.values():
        is_known = False
        for mut in interval:
            i = 0
            while i < len(df['Variant cDNA name']):
                if mut == df['Variant cDNA name'][i]:
                    known.append(i)
                    is_known = True
                    break
                i += 1
        if not is_known:
            unknown.extend([interval[0]])
    new_data = {
        'Variant cDNA name': [],
        'Variant Protein Name': [],
        '# alleles in database': [],
        'Allele frequency': [],
        'Influence': [],
        'Expected': [],
        'Predicted': [],
        'Production': [],
        'Processing': [],
        'Gating': [],
        'Conducting': [],
        'Insufficient': []
    }
    for var in unknown:
        new_data['Variant cDNA name'].append(var)
        new_data['Variant Protein Name'].append('NA')
        new_data['# alleles in database'].append('NA')
        new_data['Allele frequency'].append('NA')
        new_data['Influence'].append('NA')
        new_data['Expected'].append('NA')
    for var in known:
        values = df.loc[var].values
        new_data['Variant cDNA name'].append(values[0])
        new_data['Variant Protein Name'].append(values[1])
        new_data['# alleles in database'].append(values[2])
        new_data['Allele frequency'].append(values[3])
        new_data['Influence'].append(values[4])
        new_data['Expected'].append(values[5])
    new_data['Production'].append(scores['production'])
    new_data['Processing'].append(scores['processing'])
    new_data['Gating'].append(scores['gating'])
    new_data['Conducting'].append(scores['conducting'])
    new_data['Insufficient'].append(scores['insufficient'])
    predicted = []
    for s in ('Production','Processing','Gating','Conducting','Insufficient'):
        if new_data[s][0] > 0.6:
            predicted.append(s)
    for elem in predicted:
        if new_data['Predicted'] == []:
            new_data['Predicted'].append(elem)
        else:
            new_data['Predicted'][0] += f', {elem}'
    if not predicted:
        new_data['Predicted'].append('Benign')
    new_df = pd.DataFrame(new_data)
    if not os.path.isdir(out+'output'):
        os.mkdir(out+'output')
    if os.path.isfile(out+'output/mutation_info.csv'):
        new_df.to_csv(out+'output/mutation_info.csv',mode='a',index=False,header=None)
    else:
        new_df.to_csv(out+'output/mutation_info.csv',index=False)

def scores_figure(cleaned,db,out,scores):
    df = pd.read_excel(db)
    mutation = ''
    for interval in cleaned.values():
        is_known = False
        for mut in interval:
            i = 0
            while i < len(df['Variant cDNA name']):
                if mut == df['Variant cDNA name'][i]:
                    mutation = mut
                    is_known = True
                    break
                i += 1
        if not is_known:
            mutation = interval[0]
    mutations = list(scores.keys())
    mut_scores = list(scores.values())
    x_pos = np.arange(len(mutations))
    plt.bar(x_pos,mut_scores)
    plt.title(mutation)
    plt.xlabel('Mutation Type')
    plt.ylabel('Probability')
    plt.xticks(x_pos,mutations)
    plt.ylim(0,1)
    if not os.path.isdir(out+'output/figures'):
        os.mkdir(out+'output/figures')
    if '>' in mutation[2:]:
        index = mutation.rfind('>')
        mutation = mutation[:index] + '_' + mutation[index+1:]
    plt.savefig(out+'output/figures/'+mutation[2:]+'_scores.png')
    plt.close()
