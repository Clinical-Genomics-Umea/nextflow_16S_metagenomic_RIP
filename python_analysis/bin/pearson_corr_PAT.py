import numpy as np
import sys
import pandas as pd
from pathlib import Path
from datetime import date

''' Calculates Person correlations for cancer and normal groups respectively based on PATIENTS. Put file prefix names as sys.argv[2] and [3] '''

def get_cancer_data(infile):
    ''' Read normalised table and subset to get cancer samples '''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    pat_table = indata.T
    pat_table.columns = pat_table.iloc[0]
    pat_table = pat_table.drop(pat_table.index[0])
    pat_table = pat_table[pat_table.Cancer == 'C']  
    pat_table = pat_table.drop(['Age', 'Sex', 'Stage', 'Cancer'], axis=1).astype('float')
    pat_table = pat_table.T
    return pat_table

def get_normal_data(infile):
    ''' Read normalised table and subset to get normal samples '''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    pat_table = indata.T
    pat_table.columns = pat_table.iloc[0]
    pat_table = pat_table.drop(pat_table.index[0])
    pat_table = pat_table[pat_table.Cancer == 'N']  
    pat_table = pat_table.drop(['Age', 'Sex', 'Stage', 'Cancer'], axis=1).astype('float')
    pat_table = pat_table.T
    return pat_table

def calc_corr_cancer(pat_table, out_prefix):
    ''' Takes a table of floats, calculates Pearson correlation and writes the upper triangle and lower+upper to csv's '''
    print(pat_table.shape)
    p_corr = pat_table.corr()
    p_corr.index.name = None
    p_corr.to_csv(str(out_prefix) + '_cancer_full_' + date.today().strftime('%y%m%d') + '.csv')
    p_corr = p_corr.where(np.triu(np.ones(p_corr.shape)).astype(bool))
    p_corr.to_csv(str(out_prefix) + '_cancer_upper_' + date.today().strftime('%y%m%d') + '.csv')
    return p_corr

def calc_corr_normal(pat_table, out_prefix):
    ''' Takes a table of floats, calculates Pearson correlation and writes the upper triangle and lower+upper to csv's '''
    print(pat_table.shape)
    p_corr = pat_table.corr()
    p_corr.index.name = None
    p_corr.to_csv(str(out_prefix) + '_normal_full_' + date.today().strftime('%y%m%d') + '.csv')
    p_corr = p_corr.where(np.triu(np.ones(p_corr.shape)).astype(bool))
    p_corr.to_csv(str(out_prefix) + '_normal_upper_' + date.today().strftime('%y%m%d') + '.csv')
    return p_corr    

def melt_corr_cancer(corr_table, out_prefix):
    ''' Melts a correlation table and writes that to a csv '''
    corr_table = corr_table.stack().reset_index()
    corr_table.columns = ['Row', 'Column', 'p_corr']
    corr_table.to_csv(str(out_prefix) + '_cancer_' + date.today().strftime('%y%m%d') + '.csv', index=False)

def melt_corr_normal(corr_table, out_prefix):
    ''' Melts a correlation table and writes that to a csv '''
    corr_table = corr_table.stack().reset_index()
    corr_table.columns = ['Row', 'Column', 'p_corr']
    corr_table.to_csv(str(out_prefix) + '_normal_' + date.today().strftime('%y%m%d') + '.csv', index=False)    


def main():
    infile = Path(sys.argv[1])
    pat_table_cancer = get_cancer_data(infile)
    pat_table_normal = get_normal_data(infile)
    corr_table_cancer = calc_corr_cancer(pat_table_cancer, Path(sys.argv[2]))
    corr_table_normal = calc_corr_normal(pat_table_normal, Path(sys.argv[2]))
    melt_corr_cancer(corr_table_cancer, Path(sys.argv[3]))
    melt_corr_normal(corr_table_normal, Path(sys.argv[3]))

if __name__ == "__main__":
    main()    