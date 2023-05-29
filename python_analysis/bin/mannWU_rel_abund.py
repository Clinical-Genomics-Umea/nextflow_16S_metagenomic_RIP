import pandas as pd
import os
import sys
from pathlib import Path
from datetime import date
from scipy.stats import mannwhitneyu

'''
Perform Mann-Whithney U tests on cancer group vs normal group for all ASVs (rel abundance data)   
'''


def format_input(file):
    '''Format the input csv for MWU '''
    indata = pd.DataFrame(pd.read_csv(file, sep=',', header=0))  
    indata = indata.T
    indata.columns = indata.iloc[0]
    indata = indata[1:]
    bac_data = indata.drop(columns=['Age', 'Sex', 'Stage', 'Cancer'])
    bac_data = bac_data.astype('float')
    bac_data = bac_data.T
    return bac_data

def clin_info(meta_file, clin_var):
    ''' Output clin data for samples included '''
    metadata = pd.DataFrame(pd.read_csv(meta_file, sep=',', header=0))
    metadata['Cancer'] = metadata['Cancer'].replace(['C','N'], ['1','0'])
    clin_group = metadata[['Novogene_ID', clin_var]]
    # drop samples with Nan in clinical variable
    clin_group = clin_group[(clin_group[clin_var] != 'Nan')]
    return clin_group



def mannWU(clin_group, bac_data, out_prefix, clin_var):
    '''Runs mannWU on groups based on clinical variable'''
    stats = []
    for ASV in range(len(bac_data)):
        asv_data = bac_data.iloc[[ASV]]
        asv_data = asv_data.T.reset_index()
        asv_data.rename(columns={'index':'Novogene_ID'}, inplace=True)
        asv_groups = clin_group.merge(asv_data, on='Novogene_ID')
        means_asv = asv_groups.groupby(asv_groups.columns[1]).mean(numeric_only=True)
        zeros_mean = means_asv[:1].iloc[0][0]
        ones_mean = means_asv[1:2].iloc[0][0]
        ones_data = asv_groups[asv_groups[clin_var] == '1'].drop('Novogene_ID', axis=1)
        zeros_data = asv_groups[asv_groups[clin_var] == '0'].drop('Novogene_ID', axis=1)
        for_mann_wu = pd.concat([ones_data, zeros_data], ignore_index=True, axis=1).rename(columns={1:'ones_data', 3:'zeroes_data'})
        stat, p = mannwhitneyu(for_mann_wu['ones_data'].dropna(), for_mann_wu['zeroes_data'].dropna())
        stats.append(bac_data.index[ASV] + ' ' + str(stat) + ' ' + str(p) + ' ' + str(ones_mean) + ' ' + str(zeros_mean))
    stats = pd.DataFrame(stats)
    stats.columns = ['ASV']
    stats[['ASV', 'MannWU_stat', 'p_value', str(clin_var) + '_pos_mean', str(clin_var) + '_neg_mean',]] = stats['ASV'].str.split(' ', expand=True)
    stats = stats.astype({'MannWU_stat': float, 'p_value': float, str(clin_var) + '_pos_mean': float, str(clin_var) + '_neg_mean': float,})
    stats = stats.sort_values('p_value', ascending=False)
    stats.to_csv(str(out_prefix) + '_' + clin_var + '_' + date.today().strftime('%y%m%d') + '.csv', index=False)

def main():
    infile = Path(sys.argv[1])
    meta_file = Path(sys.argv[2])
    out_prefix = Path(sys.argv[3])
    bac_data = format_input(infile)
    # Cancer
    clin_cancer = clin_info(meta_file, 'Cancer')
    mannWU(clin_cancer, bac_data, out_prefix, 'Cancer')
    # MSI
    clin_msi = clin_info(meta_file, 'MSI')
    mannWU(clin_msi, bac_data, out_prefix, 'MSI')
   

if __name__ == "__main__":
    main()