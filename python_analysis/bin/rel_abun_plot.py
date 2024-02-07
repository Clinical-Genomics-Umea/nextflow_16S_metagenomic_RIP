import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import colorcet as cc
from datetime import date
from pathlib import Path


''' 
Plot ASV abundance on phylum and genus level. Relative abundance is calculated per group (C /N)
'''


def group_phylum(infile):
    ''' Sum counts on phylum level for cancer and normal groups '''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    clin = indata.tail(4)
    indata.drop(indata.tail(4).index,inplace=True)
    phylum_names = indata['Unnamed: 0'].str.split('_').str[1]
    indata = indata.drop(columns='Unnamed: 0').astype(float)
    indata['phylum'] = phylum_names
    phyl_sums = indata.groupby('phylum').sum()
    phyl_sums = pd.concat([phyl_sums,clin], join='inner')
    phyl_sums.rename({phyl_sums.index[-1]: 'Cancer', phyl_sums.index[-2]: 'Stage', 
                      phyl_sums.index[-3]: 'Sex', phyl_sums.index[-4]: 'Age'}, inplace=True)              
    phyl_sums = phyl_sums.T
    phyl_sums = phyl_sums.loc[:,(phyl_sums.columns!='Age') & (phyl_sums.columns!='Stage') & (phyl_sums.columns!='Sex')].groupby('Cancer').sum()
    rel_phyl = phyl_sums.T.div(phyl_sums.sum(axis=1), axis=1)
    return rel_phyl.T


def group_genus(infile):
    ''' 
    Sum counts on genus level for cancer and normal groups
    Note that all NAs at genus level will be grouped  
    '''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    clin = indata.tail(4)
    indata.drop(indata.tail(4).index,inplace=True)
    genus_names = indata['Unnamed: 0'].str.split('_').str[5]
    indata = indata.drop(columns='Unnamed: 0').astype(float)
    indata['genus'] = genus_names
    gen_sums = indata.groupby('genus').sum()
    gen_sums = pd.concat([gen_sums,clin], join='inner')
    gen_sums.rename({gen_sums.index[-1]: 'Cancer', gen_sums.index[-2]: 'Stage', 
                      gen_sums.index[-3]: 'Sex', gen_sums.index[-4]: 'Age'}, inplace=True)
    gen_sums = gen_sums.T
    gen_sums = gen_sums.loc[:,(gen_sums.columns!='Age') & (gen_sums.columns!='Stage') & (gen_sums.columns!='Sex')].groupby('Cancer').sum()

    if 'NA' in gen_sums.columns:
        gen_sums = gen_sums.rename(columns={'NA':'Other'})
    else:
        pass
    kept_gen_sums = gen_sums[gen_sums.columns[gen_sums.max() > 114000]] # filter out the smallest genuses since they will not be visible in plot, hence remove from legend 
    remove_gen_sums = gen_sums[gen_sums.columns[~(gen_sums.max() > 114000)]]
    tmp = pd.DataFrame(remove_gen_sums.sum(axis=1))
    out_gen_sums = kept_gen_sums.copy()
    out_gen_sums['Other']  = kept_gen_sums['Other'] + tmp[0]
    tmp = out_gen_sums.pop('Other')
    out_gen_sums = pd.concat([out_gen_sums, tmp], 1)
    rel_gen = out_gen_sums.T.div(out_gen_sums.sum(axis=1), axis=1)   
    return rel_gen.T

    
def make_barplot(sums, type, out_prefix):
    ''' Make a stacked barplot '''
    sums.rename(index={'C': 'Cancer', 'N': 'Normal'}, inplace=True)
    my_colors = ListedColormap(sns.color_palette(cc.glasbey, n_colors=100))
    bars = sums.plot(kind='bar', stacked=True, title = 'By ' + type, colormap=my_colors)
    plt.xticks(rotation=0)
    bars.set(xlabel=None)
    plt.subplots_adjust(right=0.6) # make space to the right of figure to place legend
    bars.legend(loc='center', bbox_to_anchor=(1.3, 0.5), fontsize=5 , ncol=1)
    plt.savefig(str(out_prefix) + type + date.today().strftime('%y%m%d') + '.png', dpi=400)
    plt.show()
    plt.close()
    


def main():
    infile = Path(sys.argv[1])
    out_prefix = Path(sys.argv[2])
    # Phylum
    phyl_sums = group_phylum(infile)
    make_barplot(phyl_sums, '_phylum_', out_prefix)
    # Genus
    genus_sums = group_genus(infile)
    make_barplot(genus_sums, '_genus_', out_prefix)




if __name__ == "__main__":
    main()