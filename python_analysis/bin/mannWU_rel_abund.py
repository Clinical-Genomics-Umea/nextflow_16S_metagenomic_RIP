import pandas as pd
import os
import sys
from pathlib import Path
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
    indata = indata.sort_values(by='Cancer')
    sample_groups = indata[['Cancer']]
    sample_groups = sample_groups.drop(indata.index[0]).reset_index()
    sample_groups.rename(columns={'index':'SampleID'}, inplace=True)

    bac_data = indata.drop(columns=['Age', 'Sex', 'Stage', 'Cancer'])
    bac_data = bac_data.astype('float')
    bac_data = bac_data.T
    return sample_groups, bac_data

def mannWU(sample_groups, bac_data, outfile):
    '''Runs mannWU on cancer vs normal groups '''
    stats = []
    for ASV in range(len(bac_data)):
        asv_data = bac_data.iloc[[ASV]]
        asv_data = asv_data.T.reset_index()
        asv_data.rename(columns={'index':'SampleID'}, inplace=True)
        asv_groups = sample_groups.merge(asv_data, on='SampleID')
        means_asv = asv_groups.groupby(asv_groups.columns[1]).mean(numeric_only=True)
        c_mean = means_asv[:1].iloc[0][0]
        n_mean = means_asv[1:2].iloc[0][0]
        c_data = asv_groups[asv_groups['Cancer'] == 'C'].drop('SampleID', axis=1)
        n_data = asv_groups[asv_groups['Cancer'] == 'N'].drop('SampleID', axis=1)
        for_mann_wu = pd.concat([c_data, n_data], ignore_index=True, axis=1).rename(columns={1:'c_data', 3:'n_data'})
        stat, p = mannwhitneyu(for_mann_wu['c_data'].dropna(), for_mann_wu['n_data'].dropna())
        stats.append(bac_data.index[ASV] + ' ' + str(stat) + ' ' + str(p) + ' ' + str(c_mean) + ' ' + str(n_mean))
    stats = pd.DataFrame(stats)
    stats.columns = ['ASV']
    stats[['ASV', 'MannWU_stat', 'p_value', 'mean_cancer', 'mean_normal',]] = stats['ASV'].str.split(' ', expand=True)
    stats = stats.astype({'MannWU_stat': float, 'p_value': float, 'mean_cancer': float, 'mean_normal': float,})
    stats = stats.sort_values('p_value', ascending=False)
    stats.to_csv(outfile, index=False)

def main():
    infile = Path(sys.argv[1])
    outfile = Path(sys.argv[2])
    sample_groups, bac_data = format_input(infile)
    mannWU(sample_groups, bac_data, outfile)
   

if __name__ == "__main__":
    main()