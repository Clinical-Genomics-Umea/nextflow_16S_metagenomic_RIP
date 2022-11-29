#!/usr/bin/env python

import sys
import pandas as pd
from pathlib import Path

pd.options.display.max_colwidth = 400
pd.options.display.max_rows = 50
pd.set_option('display.max_columns', 12)
pd.set_option('display.width', 0)

def filter(infile, outfile):
    """Filter Lulu output - remove ASVs with phylum "NA" and ASVs found in <10% of samples"""
    dada = pd.DataFrame(pd.read_csv(infile, sep='\t', header=0))
    
    # Check for ASVs with phylum 'NA' and remove them
    dada_filt = dada.loc[~((dada.index.str.contains("Bacteria_NA")) |  \
                           (dada.index.str.startswith("NA")))] 

    # Remove ASVs present in <10% of samples
    out_data = dada_filt[(dada_filt==0).sum(axis=1)/len(dada_filt.columns) <=  0.9] 
    out_data.to_csv(outfile, header=True)

def main():
    infile = Path(sys.argv[1])
    outfile = Path(sys.argv[2])
    filter(infile, outfile)

if __name__ == "__main__":
    main()