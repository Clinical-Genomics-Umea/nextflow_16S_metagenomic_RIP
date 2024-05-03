import sys
import csv 
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import date

'''
Script to calculate and plot the mean and std of p.micra abundance for each community
'''

def get_counts(infile, comm_file):
    ''' Get p.micra counts for samples in communities  '''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    communities = pd.DataFrame(pd.read_csv(comm_file, sep=',', header=0))
    indata.set_index('Unnamed: 0', inplace=True)
    micra = indata.loc['Bacteria_Firmicutes_Clostridia_Peptostreptococcales-Tissierellales_Family-XI_Parvimonas_micra.1']
    asv_dict = dict(zip(micra.index, micra.values))
    tmp = [communities[col].str.split(':', expand=True) for col in communities.columns]
    tmp = pd.concat(tmp, axis=1).drop(columns=[0])
    tmp.columns = communities.columns
    tmp = [tmp[col].str.strip() for col in tmp.columns]
    tmp = pd.concat(tmp, axis=1)
    communities = tmp.replace(asv_dict)
    return communities


def calc_plot_stats(communities, out_prefix):
    ''' Calculate and plot the mean and std for each community '''
    communities = communities.astype(float)
    mean = communities.mean()
    std = communities.std()
    mean.plot(kind='bar', legend=False,  figsize=(10, 5), yerr=std)
    plt.ylabel('Mean rel abundance')
    plt.title('Mean rel abundance of p.micra for each community')
    plt.savefig(str(out_prefix) + '_rel_abun_p.micra_' + date.today().strftime('%y%m%d') + '.png', dpi=200)
    plt.close()   

def main():
    infile = Path(sys.argv[1])
    comm_file = Path(sys.argv[2])
    out_prefix = Path(sys.argv[3])
    counts = get_counts(infile, comm_file)
    calc_plot_stats(counts, out_prefix)


if __name__ == '__main__':
    main()
