import sys
import csv 
import pandas as pd
from pathlib import Path
from datetime import date
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import numpy as np

pandas2ri.activate()
fisher = robjects.r['fisher.test']

def clin_var(infile1, infile2):
    ''' Output clin data for samples included in cluster analysis '''
    metadata = pd.DataFrame(pd.read_csv(infile1, sep=',', header=0))
    memlist = pd.DataFrame(pd.read_csv(infile2, sep=',', header=0))
    memlist = memlist.drop(columns= 'Unnamed: 0')
    memlist = memlist.rename(columns={'ASV': 'Novogene_ID'})
    clin_out = metadata.merge(memlist, on='Novogene_ID')
    return clin_out


def comm_mem(infile3, clin_out):
    ''' For each community member group .... '''

    clin = clin_out
    
    # for all clin variables:
    # for var in range(1,len(clin.columns)):
    #     print('Running var ' + str(var))
    #     print(clin.columns[var])
    #     tot = clin.iloc[:,var].value_counts().to_dict()
    #     sum_tot = sum(tot.values())   
    #     tot_frac = {key: value / sum_tot for key, value in tot.items()}

    #     indata = pd.DataFrame(pd.read_csv(infile3, sep=',', header=0))
        

    #     # for all communities:
    #     for col in range(len(indata.columns)):   
    #         comm = pd.DataFrame(indata.iloc[:,col])
    #         comm.columns = ['sample']
    #         comm[['num', 'Novogene_ID']] = comm['sample'].str.split(':', expand=True)
    #         comm = comm.drop(columns=['sample', 'num'])
    #         comm = pd.DataFrame((comm['Novogene_ID'].str.strip()))
    #         comm = comm.merge(clin)
    #         if len(comm) == 1: # if only one member in community -> skip
    #             pass
    #         else:    
    #             comm = comm.iloc[:,1].value_counts().to_dict()
    #             for key in tot.keys():
    #                 if key not in comm.keys():
    #                     comm[key] = 0
    #             comm_tot = sum(comm.values()) 
    #             comm_frac = {key: value * comm_tot for key, value in tot_frac.items()}

    #             for_fisher = pd.DataFrame({'real': list(comm.values()), 'prob': list(comm_frac.values())})
    #             print(for_fisher)
    #             test = fisher(for_fisher)
    #             fisher_p = test[0][0]
    #             print('Community ' + str(col) + ': ' + str(fisher_p))        
     

    tot = clin['Local_3grp'].value_counts().to_dict()
    sum_tot = sum(tot.values())   
    tot_frac = {key: value / sum_tot for key, value in tot.items()}

    indata = pd.DataFrame(pd.read_csv(infile3, sep=',', header=0))
    

    # # for all communities:
    # for col in range(len(indata.columns)):   
    comm = pd.DataFrame(indata.iloc[:,0])
    comm.columns = ['sample']
    comm[['num', 'Novogene_ID']] = comm['sample'].str.split(':', expand=True)
    comm = comm.drop(columns=['sample', 'num'])
    comm = pd.DataFrame((comm['Novogene_ID'].str.strip()))
    comm = comm.merge(clin)
    if len(comm) == 1: # if only one member in community -> skip
        pass
    else:    
        comm = comm.iloc[:,1].value_counts().to_dict()
        for key in tot.keys():
            if key not in comm.keys():
                comm[key] = 0
        comm_tot = sum(comm.values()) 
        comm_frac = {key: value * comm_tot for key, value in tot_frac.items()}
        print(comm.values())
        print(comm_frac.values())


    #    for_fisher = pd.DataFrame({'real': list(comm.values()), 'prob': list(comm_frac.values())})
    #     print(for_fisher)
    #     test = fisher(for_fisher)
    #     fisher_p = test[0][0]
    #     print('Community ' + str(col) + ': ' + str(fisher_p))        
    

    









def main():
    clin_file = Path(sys.argv[1])
    list_file = Path(sys.argv[2])
    comm_file = Path(sys.argv[3])
    clin_info = clin_var(clin_file, list_file)
    comm_info = comm_mem(comm_file, clin_info)


if __name__ == "__main__":
    main() 