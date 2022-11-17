#!/usr/bin/env python

from pathlib import Path
import sys
import pandas as pd

def filter(infile, outfile):
    """
    Collapse duplicate ASVs by summing up their counts
    """
    table_in = pd.DataFrame(pd.read_csv(infile, sep='\t', header=0))
    filt = table_in.groupby("OTU_ID", sort=False).sum()
    filt.to_csv(outfile, header=True, sep='\t')

def main():
    infile = Path(sys.argv[1])
    outfile = Path(sys.argv[2])
    filter(infile, outfile)

if __name__ == "__main__":
    main()