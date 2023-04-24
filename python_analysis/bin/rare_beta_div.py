import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from datetime import date
from pathlib import Path
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova
from skbio.stats.distance import DistanceMatrix


pandas2ri.activate()
vegan = importr('vegan')


def rare_beta(infile):
    ''' Make a rarefied (rarefaction is included in avgdist()) bray-curtis distance matrix  '''   
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    indata.drop(indata.tail(4).index,inplace=True)
    indata = indata.T
    indata.columns = indata.iloc[0]
    indata = indata.drop(indata.index[0]).astype(float)
    side_len = len(indata)
    n_seqs = indata.sum(axis=1).min()
    r_data = pandas2ri.py2rpy(indata)
    dist = vegan.avgdist(r_data, sample = n_seqs, iterations = 100, dmethod = "bray")
    
    # creating dissimilarity matrix from numpy 1D array
    a = np.tril(np.ones((side_len, side_len), dtype='int'), -1)
    b = np.cumsum(np.roll(a.sum(axis=0), 1))
    b[-1] = 0
    lower = np.cumsum(a + np.diagflat(b), axis=0)
    lower = np.tril(lower, -1)
    lower = np.tril(dist[lower - 1], -1)
    lower_df = pd.DataFrame(lower)
    beta_matrix = lower_df + lower_df.T
    beta_matrix.columns = indata.index
    return beta_matrix        


def clin_info(meta_file, clin_var):
    ''' Output clin data for samples included '''
    metadata = pd.DataFrame(pd.read_csv(meta_file, sep=',', header=0))
    clin_var = metadata[['Novogene_ID', clin_var]]
    return clin_var


def run_pcoa(beta_matrix, clin_var, out_prefix):
    ''' Calculate principal coordinates and plot '''
    dist_matrix = beta_matrix.copy()
    dist_matrix.insert(0, 'Novogene_ID', dist_matrix.columns)
    dist_matrix = dist_matrix.merge(clin_var, on='Novogene_ID')
    filt = dist_matrix.columns[-1]  

    # drop samples with Nan in clinical variable
    drop = dist_matrix[~(dist_matrix[filt] != 'Nan')]
    drop = drop['Novogene_ID']
    dist_matrix = dist_matrix[dist_matrix[filt] != 'Nan'] 
    dist_matrix = dist_matrix.drop(columns=drop)
    dist_matrix = dist_matrix.drop(columns = ['Novogene_ID', filt])

    pcoa_values = pcoa(dist_matrix.values)
    plot_df = pcoa_values.samples[['PC1', 'PC2']]
    plot_df.insert(0, 'Novogene_ID', dist_matrix.columns)
    plot_df = plot_df.merge(clin_var, on='Novogene_ID')
    plot_df = plot_df.drop(columns='Novogene_ID')
    group_codes = {k:idx for idx, k in enumerate(plot_df.iloc[:,2].unique())}
    colors = plot_df.iloc[:,2].apply(lambda x : group_codes[x])
    fig, ax = plt.subplots()
    scatter = ax.scatter(plot_df['PC1'], plot_df['PC2'], c=colors)
    handles = scatter.legend_elements(num=[0,1,2])[0]
    plt.subplots_adjust(right=0.8)
    ax.legend(handles=handles, labels=group_codes.keys(), loc='right', bbox_to_anchor=(1.2, 0.9), fontsize=8)
    percent_pc1 = pcoa_values.proportion_explained['PC1']
    percent_pc2 = pcoa_values.proportion_explained['PC2']
    plt.xlabel('PC1 (' + str(percent_pc1.round(2)) + ') %')
    plt.ylabel('PC2 (' + str(percent_pc2.round(2)) + ') %')
    plt.title(clin_var.columns[1])
    plt.savefig(str(out_prefix) + 'pcoa_' + clin_var.columns[1] + '_' + date.today().strftime('%y%m%d') + '.png', dpi=200)
    plt.close()

def run_permanova(beta_matrix, clin_var):
    ''' Run permanova test on distance matrix. Output permanova p-value and df for plotting '''
    dist_matrix = beta_matrix.copy()
    dist_matrix.insert(0, 'Novogene_ID', dist_matrix.columns)
    dist_matrix = dist_matrix.merge(clin_var, on='Novogene_ID')
    filt = dist_matrix.columns[-1]  

    # drop samples with Nan in clinical variable
    drop = dist_matrix[~(dist_matrix[filt] != 'Nan')]
    drop = drop['Novogene_ID']
    dist_matrix = dist_matrix[dist_matrix[filt] != 'Nan'] 
    dist_matrix = dist_matrix.drop(columns=drop)
    dist_matrix = dist_matrix.drop(columns = ['Novogene_ID', filt])

    dist_matrix.insert(0, 'Novogene_ID', dist_matrix.columns)
    dist_matrix = dist_matrix.merge(clin_var, on='Novogene_ID')
    dist_matrix = dist_matrix.drop(columns = ['Novogene_ID'])
    groups = dist_matrix.iloc[:,-1]
    samples = dist_matrix.columns.to_list()
    samples = samples[0:-1]
    dist_matrix[dist_matrix.columns[-1]] = pd.factorize(dist_matrix[dist_matrix.columns[-1]])[0]
    num_groups = list(dist_matrix.iloc[:,-1].unique())
    plot_df = pd.DataFrame()
    for group in range(len(num_groups)):
        group_p = str(group) + '_group'
        group_p = dist_matrix[dist_matrix.iloc[:,-1] == group]
        group_p = group_p.drop(group_p.columns[-1], axis =1)
        plot_df[group] = pd.concat((group_p[col] for col in group_p.columns), ignore_index=True).to_frame()
        
    perm_matrix = dist_matrix.drop(columns = [filt]).to_numpy()
    perm_matrix = perm_matrix.copy(order='C')
    dm = DistanceMatrix(perm_matrix, ids=samples)
    perm = permanova(dm, grouping=groups, permutations=999)
    return perm['p-value'], plot_df     


def make_boxplot(plot_df, perm_p, out_prefix, name_0, name_1):
    ''' Make boxplot for dissimilarity matrix values showing p-value from permanova '''
    plot_df = plot_df.rename(columns={0: name_0, 1: name_1})
    perm_p = str(perm_p.round(4))
    bx = sns.boxplot(data=plot_df, palette='Blues')
    plt.ylabel('Rarefied Bray-Curtis distance')
    plt.text(1, 1.05 , 'p_value: ' + perm_p, fontsize=10, horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)
    #plt.title('')
    plt.savefig(str(out_prefix) + 'boxplot_' + name_0 + '_' + name_1 + '_' + date.today().strftime('%y%m%d') + '.png', dpi=200)
    plt.close()


def main():
    infile = Path(sys.argv[1])
    meta_file = Path(sys.argv[2])
    out_prefix = Path(sys.argv[3])   
    bray_matrix = rare_beta(infile)
    # Cancer vs normal
    cancer_info = clin_info(meta_file, 'Cancer')
    run_pcoa(bray_matrix, cancer_info, out_prefix)
    perm_p_cancer, cancer_df = run_permanova(bray_matrix, cancer_info)
    make_boxplot(cancer_df, perm_p_cancer, out_prefix, 'Cancer', 'Normal')
    # MSI vs MSS
    msi_info = clin_info(meta_file, 'MSI')
    run_pcoa(bray_matrix, msi_info, out_prefix)
    perm_p_msi, msi_df = run_permanova(bray_matrix, msi_info)
    make_boxplot(msi_df, perm_p_msi, out_prefix, 'MSS', 'MSI')
    # Cancer specific death
    death_info = clin_info(meta_file, 'Can_spec_death')
    run_pcoa(bray_matrix, death_info, out_prefix)
    perm_p_death, death_df = run_permanova(bray_matrix, death_info)
    make_boxplot(death_df, perm_p_death, out_prefix, 'No_can_spec_death', 'Can_spec_death') # Är 0 = levande?
    # KRAS
    kras_info = clin_info(meta_file, 'KRAS')
    run_pcoa(bray_matrix, kras_info, out_prefix)
    perm_p_kras, kras_df = run_permanova(bray_matrix, kras_info)
    make_boxplot(kras_df, perm_p_kras, out_prefix, 'KRAS-', 'KRAS+')
    # BRAF
    braf_info = clin_info(meta_file, 'BRAF')
    run_pcoa(bray_matrix, braf_info, out_prefix)
    perm_p_braf, braf_df = run_permanova(bray_matrix, braf_info)
    make_boxplot(braf_df, perm_p_braf, out_prefix, 'BRAF-', 'BRAF+')
    

if __name__ == "__main__":
    main()       