#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
import numpy as np
import fastcluster
import sys
from pathlib import Path

pd.options.display.max_colwidth = 10
pd.options.display.max_rows = 50
pd.set_option('display.max_columns', 12)
pd.set_option('display.width', 0)

def make_clust_hm(infile):
    """ Subset data for plotting and make a clustered heatmap """
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    indata = indata.T  
    indata.columns = indata.iloc[0]
    indata = indata.drop(indata.index[0])
    indata = indata.reset_index()
    indata.rename(columns = {'index': 'Sample_ID'}, inplace=True)
    indata = indata.sort_values(['Cancer', 'Sample_ID'])  
    indata = indata.T
    indata.columns = indata.iloc[0]
    indata = indata.drop(indata.index[0])
    plotdata = indata.drop(['Cancer', 'Sex', 'Age', 'Stage'])
    log_plotdata = plotdata.copy()    
    log_plotdata = log_plotdata.astype('float')  
    log_plotdata += 1
    log_plotdata = np.log10(log_plotdata)

    

    # Create heatmap using log values with pseudo count +1

    # Cluster ASVs:
    sns.set(font_scale=0.5)
    cols = (['grey'] * indata.loc['Cancer'].value_counts()['C']) + (['black'] * indata.loc['Cancer'].value_counts()['N'])
  
    hm = sns.clustermap(log_plotdata, standard_scale=None, col_cluster=False, row_cluster=True,
                        method='ward', metric='euclidean',
                        cmap='vlag',
                        figsize=(20, 30),
                        col_colors=cols,
                       )
    plt.savefig("/home/lindak/project/nextflow_16S_metagenomic_RIP/plots/hm_norm_clust_ASV_noS9_NoNormASVs_221222.png", dpi=1000)

    # # Cluster ASVs and samples:
    # sns.set(font_scale=0.5)
    # cols = (['grey'] * indata.loc['Cancer'].value_counts()['C']) + (['black'] * indata.loc['Cancer'].value_counts()['N'])
    # hm = sns.clustermap(log_plotdata, standard_scale=0, col_cluster=True,row_cluster=True,
    #                     method='ward', metric='euclidean',
    #                     cmap='vlag',
    #                     figsize=(20, 30),
    #                     col_colors=cols,
    #                    )
    # plt.savefig("/home/lindak/project/nextflow_16S_metagenomic_RIP/plots/clust_hm_ASV_sample_clust_noS9_221201.png", dpi=1000)

    # # Cluster samples:
    # sns.set(font_scale=0.5)
    # cols = (['grey'] * indata.loc['Cancer'].value_counts()['C']) + (['black'] * indata.loc['Cancer'].value_counts()['N'])
    # hm = sns.clustermap(log_plotdata, standard_scale=0, col_cluster=True, row_cluster=False,
    #                     method='ward', metric='euclidean',
    #                     cmap='vlag',
    #                     figsize=(20, 30),
    #                     col_colors=cols,
    #                    )

    # plt.savefig("/home/lindak/project/nextflow_16S_metagenomic_RIP/plots/clust_hm_sample_clust_noS9_221201.png", dpi=1000)


def main():
    indata = Path(sys.argv[1])
    make_clust_hm(indata)

if __name__ == "__main__":
    main()