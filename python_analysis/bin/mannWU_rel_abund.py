import pandas as pd
import os
import sys
from pathlib import Path
from datetime import date
from scipy.stats import mannwhitneyu
from cliffs_delta import cliffs_delta
from math import sqrt
import matplotlib.pyplot as plt

'''
Perform Mann-Whithney U tests on cancer/normal and MSI/MSS for all ASVs (rel abundance data)   
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
    #for ASV in range(len(bac_data)):
    for ASV in range(357, 358):
        asv_data = bac_data.iloc[[ASV]]
        asv_data = asv_data.T.reset_index()
        asv_data.rename(columns={'index':'Novogene_ID'}, inplace=True)
        asv_groups = clin_group.merge(asv_data, on='Novogene_ID')
        print(asv_groups)
        asv_groups_for_cliffs = asv_groups.groupby(asv_groups.columns[1])[asv_groups.columns[2]]
        cliffs_d, _ = cliffs_delta(asv_groups_for_cliffs.get_group('0'), asv_groups_for_cliffs.get_group('1')) 
        print(cliffs_d)
        #test_groups = pd.DataFrame({'MSI': [0,0,0,0,0,1,1,1,1,1], 'bac': [12,15,21,22,23,18,11,18,23,17]})
        #print(test_groups)
        #test_grouped = test_groups.groupby('MSI')['bac']
        #test_grouped = test_groups.groupby(test_groups.columns[0])[test_groups.columns[1]]
        #cliffs_dtest, _ = cliffs_delta(test_grouped.get_group(0), test_grouped.get_group(1))
        #print(float(cliffs_dtest))
        means_asv = asv_groups.groupby(asv_groups.columns[1]).mean(numeric_only=True)
        std_asv = asv_groups.groupby(asv_groups.columns[1]).std(numeric_only=True)
        zeros_mean = means_asv[:1].iloc[0][0]
        ones_mean = means_asv[1:2].iloc[0][0]
        pooled_std = sqrt(((std_asv[:1].iloc[0][0])**2 + (std_asv[1:2].iloc[0][0])**2) / 2)
        cohens_d = (ones_mean - zeros_mean) / pooled_std
        ones_data = asv_groups[asv_groups[clin_var] == '1'].drop('Novogene_ID', axis=1)
        zeros_data = asv_groups[asv_groups[clin_var] == '0'].drop('Novogene_ID', axis=1)
        for_mann_wu = pd.concat([ones_data, zeros_data], ignore_index=True, axis=1).rename(columns={1:'ones_data', 3:'zeroes_data'})
        stat, p = mannwhitneyu(for_mann_wu['ones_data'].dropna(), for_mann_wu['zeroes_data'].dropna())
        stats.append(bac_data.index[ASV] + ' ' + str(stat) + ' ' + str(p) + ' ' + str(ones_mean) + ' ' + str(zeros_mean) + ' ' + str(pooled_std) + ' ' + str(cohens_d))   
    stats = pd.DataFrame(stats)
    stats.columns = ['ASV']
    stats[['ASV', 'MannWU_stat', 'p_value', str(clin_var) + '_pos_mean', str(clin_var) + '_neg_mean', str(clin_var) + '_pooled_stdev', str(clin_var) + '_cohens_d',]] = stats['ASV'].str.split(' ', expand=True)
    stats = stats.astype({'MannWU_stat': float, 'p_value': float, str(clin_var) + '_pos_mean': float, str(clin_var) + '_neg_mean': float, str(clin_var) + '_pooled_stdev': float,  str(clin_var) + '_cohens_d': float,})
    stats = stats.sort_values(str(clin_var) + '_cohens_d', ascending=False)
    #stats.to_csv(str(out_prefix) + '_' + clin_var + '_' + date.today().strftime('%y%m%d') + '.csv', index=False)
    print(stats)
    return stats

def plot_cohens_d(stats, out_prefix, clin_var):
    '''Plot cohen's d values sorted on cohen's d'''
    #cohens_values = pd.DataFrame(stats[['ASV', str(clin_var) + '_cohens_d']])
    cohens_values = pd.DataFrame(stats[['ASV', 'p_value', str(clin_var) + '_cohens_d']])        # for p_value sorting
    cohens_values.set_index('ASV', inplace=True)
    #cohens_values = cohens_values.nlargest(50, str(clin_var) + '_cohens_d')
    cohens_values = cohens_values.nsmallest(50, 'p_value')                                      # for p_value sorting
    cohens_values = cohens_values.sort_values(['p_value'], ascending=False)                     # for p_value sorting
    cohens_values = cohens_values.drop(columns='p_value')                                       # for p_value sorting
    cohens_values = cohens_values[str(clin_var) + '_cohens_d'].astype(float)
    #cohens_values = cohens_values.sort_values(ascending=False)   
    cohens_values.plot(kind='barh', x=str(clin_var) + '_cohens_d', y='ASV', rot=0, figsize=(10,8))
    plt.xlabel("Cohen's d")
    plt.ylabel(None)
    plt.subplots_adjust(left=0.6)
    plt.xlim([-0.8, 0.8])
    plt.yticks(fontsize=6)
    #plt.savefig(str(out_prefix) + '_cohens_d_' + clin_var + '_' + date.today().strftime('%y%m%d') + '.png', dpi=200)
    #plt.savefig(str(out_prefix) + '_cohens_d_sorted_on_p_value' + clin_var + '_' + date.today().strftime('%y%m%d') + '.png', dpi=200)    # for p_value sorting
    plt.show()
    plt.close()

def main():
    infile = Path(sys.argv[1])
    meta_file = Path(sys.argv[2])
    out_prefix = Path(sys.argv[3])
    bac_data = format_input(infile)
    # Cancer
    clin_cancer = clin_info(meta_file, 'Cancer')
    stats_cancer = mannWU(clin_cancer, bac_data, out_prefix, 'Cancer')
    plot_cohens_d(stats_cancer, out_prefix, 'Cancer')
    # MSI
    clin_msi = clin_info(meta_file, 'MSI')
    stats_msi = mannWU(clin_msi, bac_data, out_prefix, 'MSI')
    plot_cohens_d(stats_msi, out_prefix, 'MSI')
   

if __name__ == "__main__":
    main()