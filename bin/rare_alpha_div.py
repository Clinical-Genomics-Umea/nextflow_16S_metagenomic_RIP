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
from scipy.stats import mannwhitneyu
from skbio.diversity.alpha import shannon, simpson, ace, chao1

pandas2ri.activate()
vegan = importr('vegan')

''' 
Rarefy data (vegan) to the depth of min sum of all ASVS (iterate 100 times). 
For each iteration calculate alpha diversity indices (skbio.divversity.alpha) shannon, simpson, chao1 and ace
and calc MWU between cancer and normal groups.
'''

def rare_diversity(infile, index_type, dtype):
    ''' Rarefy and calculate diversity index for each rarefied data set. Sort indices into cancer and normal group '''   
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    clin = indata.tail(4)
    indata.drop(indata.tail(4).index,inplace=True)
    indata = indata.T
    indata.columns = indata.iloc[0]
    indata = indata.drop(indata.index[0]).astype(float)
    n_seqs = indata.sum(axis=1).min()
    r_data = pandas2ri.py2rpy(indata)
    index_df = pd.DataFrame()
    for n in range(100):
        rare_data = vegan.rrarefy(r_data, n_seqs)  # rarefy to depth of min sum of all ASVS
        rare_data =pd.DataFrame(rare_data)
        index_list = []
        for i in range(len(rare_data)):
            pat = rare_data.iloc[i]
            index_list.append(index_type(pat.astype(dtype)))    
        index_df = pd.DataFrame(index_list)
        index_df = index_df.rename(columns={0: n})
        index_df[n] = index_df
    index_df.index = indata.index
    index_df = index_df.T
    index_df = pd.concat([index_df,clin], join='inner', ignore_index=True)
    index_df.rename({index_df.index[-1]: 'Cancer', index_df.index[-2]: 'Stage', 
                     index_df.index[-3]: 'Sex', index_df.index[-4]: 'Age'}, inplace=True)
    index_df = index_df.T.drop(columns=['Age', 'Sex', 'Stage'])
    c_group = index_df[index_df.Cancer == 'C']
    c_group = index_df[index_df.Cancer == 'C'].drop(columns = 'Cancer')
    n_group = index_df[index_df.Cancer == 'N'].drop(columns = 'Cancer')
    return c_group.astype(float), n_group.astype(float)


def perform_mwu(c_group, n_group):
    ''' Calculate MWU for cancer and normal diversity indices '''
    stat, p = mannwhitneyu(c_group[0].dropna(), n_group[0].dropna())
    return stat,p


def make_boxplot(c_group, n_group, mwu_p, out_prefix, index_type):
    ''' Make boxplot for cancer and normal diversity indices and show MWU p value '''
    c_group = c_group.rename(columns={0: 'index'})
    c_group['group'] = 'cancer'
    n_group = n_group.rename(columns={0: 'index'})
    n_group['group'] = 'normal'
    df = pd.concat([c_group, n_group])
    mwu_p = str(mwu_p.round(4))
    bx = sns.boxplot(x=df['group'], y=df['index'], palette='Blues')
    plt.text(1, 1.05 , 'p_value: ' + mwu_p, fontsize=10, horizontalalignment='right', verticalalignment='bottom', transform=plt.gca().transAxes)
    plt.title(index_type)
    plt.savefig(str(out_prefix) + index_type + '_' + date.today().strftime('%y%m%d') + '.png', dpi=200)
    plt.close()
    
    
def main():
    infile = Path(sys.argv[1])
    out_prefix = Path(sys.argv[2])
    # Shannon
    c_group_shannon, n_group_shannon = rare_diversity(infile, shannon, float)
    shannon_mwu_stat, shannon_mwu_p = perform_mwu(c_group_shannon, n_group_shannon)
    make_boxplot(c_group_shannon, n_group_shannon, shannon_mwu_p, out_prefix, 'shannon')
    # Simpson
    c_group_simpson, n_group_simpson = rare_diversity(infile, simpson, float)
    simpson_mwu_stat, simpson_mwu_p = perform_mwu(c_group_simpson, n_group_simpson)
    make_boxplot(c_group_simpson, n_group_simpson, simpson_mwu_p, out_prefix, 'simpson')
    # Chao1
    c_group_chao1, n_group_chao1 = rare_diversity(infile, chao1, float)
    chao1_mwu_stat, chao1_mwu_p = perform_mwu(c_group_chao1, n_group_chao1)
    make_boxplot(c_group_chao1, n_group_chao1, chao1_mwu_p, out_prefix, 'chao1')
    # Ace
    c_group_ace, n_group_ace = rare_diversity(infile, ace, int)
    ace_mwu_stat, ace_mwu_p = perform_mwu(c_group_ace, n_group_ace)
    make_boxplot(c_group_ace, n_group_ace, ace_mwu_p, out_prefix, 'ace')   
    
    

if __name__ == "__main__":
    main()   