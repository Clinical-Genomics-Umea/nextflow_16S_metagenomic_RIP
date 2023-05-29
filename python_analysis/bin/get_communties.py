import sys
import csv 
import pandas as pd
from pathlib import Path
from datetime import date

''' Take community file and list and output file with community members names '''

def make_df(infile):
    ''' Make df from communities txt file ''' 
    with open(infile) as inf:
        lines = inf.readlines()
        comm = {}
        for line in lines:
            if line.startswith('Community'):
                key = line.strip().split(':')[0]
                comm[key] = []
            else:
                value = line.strip()
                comm[key].append(value)
    df = pd.DataFrame({key: pd.Series(value) for key, value in comm.items()})
    return df   


def id_df(df, ids, out_prefix):
    ''' Put ID names in df using IDs list'''
    comm = df.copy()
    reader = csv.reader(open(ids, 'r'))
    ids = {}
    for key, value in reader:
        ids[key] = value
    df.replace(ids, inplace=True)
    out = comm.apply(lambda col: col.astype (str)) + " : " + df     
    out.to_csv(str(out_prefix) + '_community_mem_names_' + date.today().strftime('%y%m%d') + '.csv', index=False)


def main():
    infile = Path(sys.argv[1])
    ids = Path(sys.argv[2])
    out_prefix = Path(sys.argv[3])   
    df = make_df(infile)
    id_df(df, ids, out_prefix)


if __name__ == "__main__":
    main() 