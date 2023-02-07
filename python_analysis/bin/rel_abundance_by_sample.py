import sys
import pandas as pd
from pathlib import Path

''' Normalize by calculating relative abundance - use output matrix in all donwstream analyses '''

def rel_abundance(infile, outfile):
    ''' Take a Lulu output file and calculate relative abundance by sample '''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    indata = indata.T
    indata.columns = indata.iloc[0]
    indata = indata.drop(indata.index[0])
    sample_ids = pd.DataFrame(indata.index)
    clin = indata[['Age', 'Sex', 'Stage', 'Cancer']]
    indata = indata.drop(columns=['Age', 'Sex', 'Stage', 'Cancer'])
    
    asv_names = indata.columns   
    rel_data = indata.T.astype(float)
    rel_data = rel_data.div(rel_data.sum(axis=0), axis=1)
    rel_data = pd.concat([rel_data, clin.T])
    rel_data.to_csv(outfile)

def main():
    infile = Path(sys.argv[1])
    outfile = Path(sys.argv[2])
    rel_abundance(infile, outfile)

if __name__ == "__main__":
    main()

