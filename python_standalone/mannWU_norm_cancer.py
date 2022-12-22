import pandas as pd
import os
from scipy.stats import mannwhitneyu

'''
Take a lule output file and perform Mann-Whithney U tests on cancer group vs normal group for all ASVs in file  
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
        means_asv = asv_groups.groupby(asv_groups.columns[1]).mean()
        c_mean = means_asv[:1].iloc[0][0]
        n_mean = means_asv[1:2].iloc[0][0]
        c_data = asv_groups[asv_groups['Cancer'] == 'C'].drop('SampleID', 1)
        n_data = asv_groups[asv_groups['Cancer'] == 'N'].drop('SampleID', 1)
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
    #sample_groups, bac_data = format_input("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/lulu_out_filtered_clin_noS9_221201.csv")
    #mannWU(sample_groups, bac_data, "/home/lindak/project/nextflow_16S_metagenomic_RIP/data/mannWU_filtered_noS9_221220.csv")
    sample_groups, bac_data = format_input("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/norm_lulu_filtered_221222.csv")
    mannWU(sample_groups, bac_data, "/home/lindak/project/nextflow_16S_metagenomic_RIP/data/mannWU_norm_filtered_221222.csv")

main()