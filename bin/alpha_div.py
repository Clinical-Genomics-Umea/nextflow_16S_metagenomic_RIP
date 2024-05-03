import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import date
from pathlib import Path
from skbio.diversity.alpha import shannon, simpson, ace, chao1
from scipy.stats import mannwhitneyu


'''  '''

def group_indata(infile):
    ''' Sort input into cancer and normal group '''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    indata = indata.T
    indata.columns = indata.iloc[0]
    indata = indata.drop(indata.index[0])  
    indata = indata.drop(columns=['Age', 'Sex', 'Stage'])
    c_group = indata[indata.Cancer == 'C'].drop(columns = 'Cancer')
    n_group = indata[indata.Cancer == 'N'].drop(columns = 'Cancer')
    return c_group, n_group


def calc_shannon(group):
    ''' Calculate shannon index for all samples in group '''
    group = group.astype(int)
    print(group)
    shannon_list = []
    for i in range(len(group)):
        pat = group.iloc[i].to_numpy()
        print(pat)
        print(shannon(pat))
        shannon_list.append(shannon(pat))
    shannon_df = pd.DataFrame(shannon_list)
    shannon_df.index = group.index
    return shannon_df


def calc_simpson(group):
    ''' Calculate simpson index for all samples in group '''
    group = group.astype(int)
    simpson_list = []
    for i in range(len(group)):
        pat = group.iloc[i].to_numpy()
        simpson_list.append(simpson(pat))
    simpson_df = pd.DataFrame(simpson_list)
    simpson_df.index = group.index
    return simpson_df


def calc_chao1(group):
    ''' Calculate chao1 richness estimator for all samples in group '''
    group = group.astype(int)
    chao_list = []
    for i in range(len(group)):
        pat = group.iloc[i].to_numpy()
        chao_list.append(chao1(pat))
    chao_df = pd.DataFrame(chao_list)
    chao_df.index = group.index
    return chao_df


def calc_ace(group):
    ''' Calculate ace metric for all samples in group '''
    group = group.astype(int)
    ace_list = []
    for i in range(len(group)):
        pat = group.iloc[i].to_numpy()
        ace_list.append(ace(pat))
    ace_df = pd.DataFrame(ace_list)
    ace_df.index = group.index
    return ace_df    


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
    plt.savefig(str(out_prefix) + index_type + date.today().strftime('%y%m%d') + '.png', dpi=200)
    plt.close()




def main():
    infile = Path(sys.argv[1])
    out_prefix = Path(sys.argv[2])
    c_group, n_group = group_indata(infile)
    # Shannon
    shannon_cancer_grp = calc_shannon(c_group)
    shannon_normal_grp = calc_shannon(n_group)
    shannon_mwu_stat, shannon_mwu_p = perform_mwu(shannon_cancer_grp, shannon_normal_grp)
    make_boxplot(shannon_cancer_grp, shannon_normal_grp, shannon_mwu_p, out_prefix, '_shannon_')
    # Simpson
    simpson_cancer_grp = calc_simpson(c_group)
    simpson_normal_grp = calc_simpson(n_group)
    simpson_mwu_stat, simpson_mwu_p = perform_mwu(simpson_cancer_grp, simpson_normal_grp)
    make_boxplot(simpson_cancer_grp, simpson_normal_grp, simpson_mwu_p, out_prefix, '_simpson_')
    # Chao1
    chao_cancer_grp = calc_chao1(c_group)
    chao_normal_grp = calc_chao1(n_group)
    chao_mwu_stat, chao_mwu_p = perform_mwu(chao_cancer_grp, chao_normal_grp)
    make_boxplot(chao_cancer_grp, chao_normal_grp, chao_mwu_p, out_prefix, '_chao1_')
    # Ace
    ace_cancer_grp = calc_ace(c_group)
    ace_normal_grp = calc_ace(n_group)
    ace_mwu_stat, ace_mwu_p = perform_mwu(ace_cancer_grp, ace_normal_grp)
    make_boxplot(ace_cancer_grp, ace_normal_grp, ace_mwu_p, out_prefix, '_ace_')
    
    

if __name__ == "__main__":
    main()