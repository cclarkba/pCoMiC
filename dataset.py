####### CLEANING DATA #######

# This will not be included in final pipeline, just need it to clean some data up

import pandas as pd

def non_dups(df1, df2, out):
    df1_var = df1['Variant cDNA name'].to_list()
    df2_var = df2['Variant cDNA name'].to_list()

    alt_splice = ['+', '-', '|','[']

    out_data = []

    for i in range(len(df2_var)):
        if any(x in df2_var[i] for x in alt_splice):
            continue
        app = True
        for j in range(len(df1_var)):
            if df2_var[i] == df1_var[j]:
                app = False
                break
        if app:
            out_data.append(df2_var[i])
    
    out_df = pd.DataFrame(out_data, columns=['Variant cDNA name'])
    out_df.to_excel(out, index=False)


if __name__ == '__main__':
    db1 = "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/database.xlsx"
    db2 = "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/to_clean.xlsx"

    out = "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/cleaned.xlsx"

    df1 = pd.read_excel(db1)
    df2 = pd.read_excel(db2)

    non_dups(df1, df2, out)
    