import sys
import csv 
import pandas as pd
from pathlib import Path
from datetime import date
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
#import numpy as np

pandas2ri.activate()
fisher = robjects.r['fisher.test']

def clin_var(infile1, infile2):
    ''' Output clin data for samples included in cluster analysis '''
    metadata = pd.DataFrame(pd.read_csv(infile1, sep=',', header=0))
    memlist = pd.DataFrame(pd.read_csv(infile2, sep=',', header=0))
    memlist = memlist.drop(columns= 'Unnamed: 0')
    memlist = memlist.rename(columns={'ID': 'Novogene_ID'})
    clin_out = metadata.merge(memlist, on='Novogene_ID')
    print(clin_out)
    return clin_out


def comm_mem(infile3, clin_out, out_prefix):
    '''
    For each community member group calculate Fisher's exact test for all clinical variables. 
    Note that Fisher will round exptected values to nearest integer.
    '''

    clin = clin_out
    
    
    #for all clin variables:
    comm_out = []
    for var in range(1,len(clin.columns)):
        print('Running var ' + str(clin.columns[var]))
        tot = clin.iloc[:,var].value_counts().to_dict()
        if 'Nan' in tot.keys():
            del tot['Nan']
        sum_tot = sum(tot.values())   
        tot_frac = {key: value / sum_tot for key, value in tot.items()}

        indata = pd.DataFrame(pd.read_csv(infile3, sep=',', header=0))
        

        # for all communities:     
        for col in range(len(indata.columns)):   
            comm = pd.DataFrame(indata.iloc[:,col])
            comm.columns = ['sample']
            comm[['num', 'Novogene_ID']] = comm['sample'].str.split(':', expand=True)
            comm = comm.drop(columns=['sample', 'num'])
            comm = pd.DataFrame((comm['Novogene_ID'].str.strip()))
            comm = comm.merge(clin)
            if len(comm) == 1: # if only one member in community -> skip
                pass
            else:    
                comm = comm.iloc[:,var].value_counts().to_dict()
                if 'Nan' in comm.keys():
                    del comm['Nan']
                for key in tot.keys():
                    if key not in comm.keys():
                        comm[key] = 0
                comm_tot = sum(comm.values()) 
                comm_frac = {key: value * comm_tot for key, value in tot_frac.items()}

                for_fisher = pd.DataFrame({'real': list(comm.values()), 'prob': list(comm_frac.values())})
                test = fisher(for_fisher)
                fisher_p = test[0][0]
                comm_out.append(str(clin.columns[var]) + ' ' + str(col) + ' ' + str(list(comm.values())) + ' ' + str(list(comm_frac.values())) + ' ' + str(fisher_p))
                print('Community ' + str(col) + ': ' + str(fisher_p))
    comm_out = pd.DataFrame(comm_out)              
    comm_out.columns = ['clin_var']
    comm_out[['clin_var', 'community', 'real', 'prob', 'fisher_p']] = comm_out['clin_var'].str.split('\s(?![^\[\]]*\])', expand=True)
    comm_out.to_csv(str(out_prefix) + '_fisher_communities_' + date.today().strftime('%y%m%d') + '.csv', index=False)


def main():
    clin_file = Path(sys.argv[1])
    list_file = Path(sys.argv[2])
    comm_file = Path(sys.argv[3])
    out_prefix = Path(sys.argv[4])
    clin_info = clin_var(clin_file, list_file)
    comm_mem(comm_file, clin_info, out_prefix)


if __name__ == "__main__":
    main()
