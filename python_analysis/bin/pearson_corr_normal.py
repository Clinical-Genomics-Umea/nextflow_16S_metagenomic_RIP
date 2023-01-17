import numpy as np
import sys
import pandas as pd
from pathlib import Path


def get_data(infile):
    ''' Read normalised table and subset to get cancer samples '''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    asv_table = indata.T 
    asv_table.columns = asv_table.iloc[0]
    asv_table = asv_table.drop(asv_table.index[0])
    asv_table = asv_table[asv_table.Cancer == 'N']  
    asv_table = asv_table.drop(['Age', 'Sex', 'Stage', 'Cancer'], axis=1).astype('float')   
    return asv_table

def calc_corr(asv_table, outfile):
    ''' Takes a table of floats, calculates Pearson correlation and writes the upper triangle to a csv '''
    print(asv_table.shape)
    p_corr = asv_table.corr()
    p_corr.index.name = None
    p_corr = p_corr.where(np.triu(np.ones(p_corr.shape)).astype(bool))
    p_corr.to_csv(outfile)
    return p_corr

def melt_corr(corr_table, outfile):
    ''' Melts a correlation table and writes that to a csv '''
    corr_table = corr_table.stack().reset_index()
    corr_table.columns = ['Row', 'Column', 'p_corr']
    corr_table.to_csv(outfile, index= False)


def main():
    infile = Path(sys.argv[1])
    asv_table = get_data(infile)
    corr_table = calc_corr(asv_table, Path(sys.argv[2]))
    melt_corr(corr_table, Path(sys.argv[3]))

if __name__ == "__main__":
    main()    