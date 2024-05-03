#!/usr/bin/env python

import sys
import pandas as pd
from pathlib import Path

pd.options.display.max_colwidth = 400
pd.options.display.max_rows = 50
pd.set_option('display.max_columns', 12)
pd.set_option('display.width', 0)

def add_clin(indata, clin_data, outfile):
    """Add clinical characteristics to filtered Lulu output"""
    lulu_data = pd.DataFrame(pd.read_csv(indata, sep=',', header=0))
    clin_char = pd.DataFrame(pd.read_csv(clin_data, sep= '\t', header=0, dtype={'Stage': pd.Int64Dtype()}))
        
    lulu_data = lulu_data.set_index("Unnamed: 0")
    clin_char = clin_char.T
    clin_char.columns =  clin_char.iloc[0]
    clin_char = clin_char.drop(clin_char.index[0])
   
    df_concat = pd.concat([lulu_data, clin_char], join="inner")
    df_concat.to_csv(outfile, header=True)

def main():
    indata = Path(sys.argv[1])
    clin_data = Path(sys.argv[2])
    outfile = Path(sys.argv[3])
    add_clin(indata, clin_data, outfile)

if __name__ == "__main__":
    main()