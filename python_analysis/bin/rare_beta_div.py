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
    dist = vegan.avgdist(r_data, sample = n_seqs, iterations = 10, dmethod = "bray")
  
    # creating dissimilarity matrix from numpy 1D array
    a = np.tril(np.ones((side_len, side_len), dtype='int'), -1)
    b = np.cumsum(np.roll(a.sum(axis=0), 1))
    b[-1] = 0
    lower = np.cumsum(a + np.diagflat(b), axis=0)
    lower = np.tril(lower, -1)
    lower = np.tril(dist[lower - 1], -1)
    lower_df = pd.DataFrame(lower)
    sym_df = lower_df + lower_df.T
    sym_df.columns = indata.index
    return sym_df        


def clin_info(meta_file, clin_var):
    ''' Output clin data for samples included '''
    metadata = pd.DataFrame(pd.read_csv(meta_file, sep=',', header=0))
    clin_var = metadata[['Novogene_ID', clin_var]]
    return clin_var


def run_pcoa(dist_matrix, clin_var, out_prefix):
    ''' Calculate principal coordinates and plot '''
    dist_matrix.insert(0, 'Novogene_ID', dist_matrix.columns)
    dist_matrix = dist_matrix.merge(clin_var, on='Novogene_ID')
    filt = dist_matrix.columns[-1]
    dist_matrix = dist_matrix[dist_matrix[filt] != 'Nan']
    print(dist_matrix)

    # Måste ta bort även de proverna (med Nan) från kolumnerna

    dist_matrix = dist_matrix.drop(columns = ['Novogene_ID', filt])
    print(dist_matrix)
    pcoa_values = pcoa(dist_matrix.values)
    plot_df = pcoa_values.samples[['PC1', 'PC2']]
    #plot_df.insert(0, 'Novogene_ID', dist_matrix.columns)
    #plot_df = plot_df.merge(clin_var, on='Novogene_ID')
    plot_df = plot_df.drop(columns='Novogene_ID')
    plot_df[clin_var.columns[1]] = pd.factorize(plot_df[clin_var.columns[1]])[0]
    colors = plot_df[plot_df.columns[2]].to_numpy()
    plt.scatter(plot_df['PC1'], plot_df['PC2'], c=colors)
    percent_pc1 = pcoa_values.proportion_explained['PC1']
    percent_pc2 = pcoa_values.proportion_explained['PC2']
    plt.xlabel('PC1 (' + str(percent_pc1.round(2)) + ') %')
    plt.ylabel('PC2 (' + str(percent_pc2.round(2)) + ') %')
    plt.title(clin_var.columns[1])
    plt.savefig(str(out_prefix) + clin_var.columns[1] + '_pcoa_' + date.today().strftime('%y%m%d') + '.png', dpi=200)



def main():
    infile = Path(sys.argv[1])
    meta_file = Path(sys.argv[2])
    out_prefix = Path(sys.argv[3])   
    bray_matrix = rare_beta(infile)
    # Cancer vs normal
    #cancer_info = clin_info(meta_file, 'Cancer')
    #run_pcoa(bray_matrix, cancer_info, out_prefix)
    # MSI vs MSS
    msi_info = clin_info(meta_file, 'MSI')
    run_pcoa(bray_matrix, msi_info, out_prefix)
    

if __name__ == "__main__":
    main()       