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

np.set_printoptions(linewidth=np.inf)


'''
Create a dissimilarity matrix with Bray-Curtis distances. Data is rarefied with vegan.avgdist()
to the depth of min sum of all ASVS (iterated 100 times). For each clinical variable:
cancer/normal, cancer specific death, MSI/MSS, KRAS, BRAF:
    - Principal coordinate analysis is performed and plotted
    - A permanova test is performed
'''


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
    #dist = vegan.vegdist(r_data, sample = n_seqs, dmethod = "dbray")
    
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
 
    pcoa_plot_legend = {'Cancer': ['Cancer', 'Normal'], 
                           'MSI': ['MSS', 'MSI'],
                'Can_spec_death': ['No', 'Yes'], 
                          'KRAS': ['KRAS-', 'KRAS+'],
                          'BRAF': ['BRAF-', 'BRAF+'],
                    'Local_3grp': ['1', '2', '3']}
    group_codes = {k:idx for idx, k in enumerate(plot_df.iloc[:,2].unique())}
    colors = plot_df.iloc[:,2].apply(lambda x : group_codes[x])
    fig, ax = plt.subplots()
    scatter = ax.scatter(plot_df['PC1'], plot_df['PC2'], c=colors)
    handles = scatter.legend_elements(num=[0,1,2])[0]
    plt.subplots_adjust(right=0.8)
    labels = pcoa_plot_legend[plot_df.columns[-1]]
    ax.legend(handles=handles, labels = labels, loc='right', bbox_to_anchor=(1.2, 0.9), fontsize=8)
    percent_pc1 = (pcoa_values.proportion_explained['PC1']) * 100
    percent_pc2 = (pcoa_values.proportion_explained['PC2']) * 100
    plt.xlabel('PC1 (' + str(percent_pc1.round(1)) + ') %')
    plt.ylabel('PC2 (' + str(percent_pc2.round(1)) + ') %')
    plt.title(clin_var.columns[1])
    plt.savefig(str(out_prefix) + 'pcoa_' + clin_var.columns[1] + '_' + date.today().strftime('%y%m%d') + '.png', dpi=200)
    #plt.show()
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
    samples = dist_matrix.columns.to_list()[0:-1]
    
    if dist_matrix.columns[-1] != 'Local_3grp':
        dist_matrix[dist_matrix.columns[-1]] = pd.factorize(dist_matrix[dist_matrix.columns[-1]])[0]
    elif dist_matrix.columns[-1] == 'Local_3grp':
        dist_matrix.iloc[:, -1] = dist_matrix.iloc[:, -1].astype('int')
      
    groups_cat = pd.DataFrame(dist_matrix.iloc[:,-1])
    groups_cat.index = samples
    groups_cat_dict = groups_cat.groupby([groups_cat.columns[0]]).apply(lambda x: x.index.tolist()).to_dict()
    perm_matrix = dist_matrix.drop(columns = [filt]).to_numpy()
    plot_matrix = perm_matrix.copy()
    plot_matrix = np.triu(plot_matrix, 1)
    plot_matrix[np.tril_indices(plot_matrix.shape[0])] = np.nan
    plot_matrix = pd.DataFrame(plot_matrix)
    plot_matrix.columns = samples
    plot_matrix.index = samples
    
    plot_df = pd.DataFrame()
    for key in groups_cat_dict.keys():
        group_dists = []
        sample_group = groups_cat_dict.get(key)
        for row_name in sample_group:
            for col_name in sample_group:
               group_dists.append(plot_matrix.loc[row_name][col_name])

        clean_group_dist = []
        for dist in group_dists:
            if str(dist) != 'nan':
                clean_group_dist.append(dist)
        clean_group_dist = pd.DataFrame(clean_group_dist)
        clean_group_dist.columns = [key]
        plot_df = pd.concat([plot_df, clean_group_dist], axis=1)           
    perm_matrix = perm_matrix.copy(order='C')
    dm = DistanceMatrix(perm_matrix, ids=samples)
    perm = permanova(dm, grouping=groups, permutations=999)
    return perm['p-value'], plot_df 


def make_boxplot(plot_df, perm_p, out_prefix, name_list, clin_var):
    ''' Make boxplot for beta diversity values showing p-value from permanova '''
    plot_df.columns = name_list
    perm_p = str(perm_p.round(4))
    bx = sns.boxplot(data=plot_df, palette='Blues')
    plt.ylabel('Rarefied Bray-Curtis distance')
    plt.text(1, 1.05 , 'p_value: ' + perm_p, fontsize=10, horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)
    #plt.title('')
    plt.savefig(str(out_prefix) + 'boxplot_' + clin_var.columns[1] + '_' + date.today().strftime('%y%m%d') + '.png', dpi=200)
    #plt.show()
    plt.close()


def main():
    infile = Path(sys.argv[1])
    meta_file = Path(sys.argv[2])
    out_prefix = Path(sys.argv[3])   
    bray_matrix = rare_beta(infile)
   
    ## Cancer vs normal
    cancer_info = clin_info(meta_file, 'Cancer')
    run_pcoa(bray_matrix, cancer_info, out_prefix)
    perm_p_cancer, cancer_df = run_permanova(bray_matrix, cancer_info)
    make_boxplot(cancer_df, perm_p_cancer, out_prefix, ['Cancer', 'Normal'], cancer_info)
    # Local_3grp
    local_info = clin_info(meta_file, 'Local_3grp')
    run_pcoa(bray_matrix, local_info, out_prefix)
    perm_p_local, local_df = run_permanova(bray_matrix, local_info)
    make_boxplot(local_df, perm_p_local, out_prefix, ['1', '2', '3'], local_info)
    # MSI vs MSS
    msi_info = clin_info(meta_file, 'MSI')
    run_pcoa(bray_matrix, msi_info, out_prefix)
    perm_p_msi, msi_df = run_permanova(bray_matrix, msi_info)
    make_boxplot(msi_df, perm_p_msi, out_prefix, ['MSS', 'MSI'], msi_info)
    # Cancer specific death
    death_info = clin_info(meta_file, 'Can_spec_death')
    run_pcoa(bray_matrix, death_info, out_prefix)
    perm_p_death, death_df = run_permanova(bray_matrix, death_info)
    make_boxplot(death_df, perm_p_death, out_prefix, ['No_can_spec_death', 'Can_spec_death'], death_info)
    # KRAS
    kras_info = clin_info(meta_file, 'KRAS')
    run_pcoa(bray_matrix, kras_info, out_prefix)
    perm_p_kras, kras_df = run_permanova(bray_matrix, kras_info)
    make_boxplot(kras_df, perm_p_kras, out_prefix, ['KRAS-', 'KRAS+'], kras_info)
    # BRAF
    braf_info = clin_info(meta_file, 'BRAF')
    run_pcoa(bray_matrix, braf_info, out_prefix)
    perm_p_braf, braf_df = run_permanova(bray_matrix, braf_info)
    make_boxplot(braf_df, perm_p_braf, out_prefix, ['BRAF-', 'BRAF+'], braf_info)
    

if __name__ == "__main__":
    main()       