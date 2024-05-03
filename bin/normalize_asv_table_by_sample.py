import numpy as np
import sys
import pandas as pd
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, MinMaxScaler


def norm_asvs(infile, outfile):
    ''' Take a Lulu output file and normalize by sample using MinMax scaling'''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    indata = indata.T
    indata.columns = indata.iloc[0]
    indata = indata.drop(indata.index[0])
    sample_ids = pd.DataFrame(indata.index)
    clin = indata[['Age', 'Sex', 'Stage', 'Cancer']]
    indata = indata.drop(columns=['Age', 'Sex', 'Stage', 'Cancer'])
    asv_names = indata.columns 
    indata = indata.to_numpy()
    data = indata.astype(float)

    preprocessor = Pipeline([("scaler", MinMaxScaler())])                
    pipe = Pipeline([("preprocessor", preprocessor)])
    pipe.fit(data)
    pcadf = pd.DataFrame(pipe["preprocessor"].fit_transform(data.T).T)
    
    out_df = sample_ids.merge(pcadf, left_index=True, right_index=True)
    out_df = out_df.set_index(['0_x'])
    out_df.columns = asv_names
    out_df.index.name = None
    out_df = out_df.merge(clin, left_index=True, right_index=True)
    out_df = out_df.T
    out_df.index.name = None
    out_df = out_df.reset_index()
    out_df.rename(columns={'index':'Unnamed: 0'}, inplace=True)
    out_df.to_csv(outfile, index=False)

def main():
    infile = Path(sys.argv[1])
    outfile = Path(sys.argv[2])
    norm_asvs(infile, outfile)

if __name__ == "__main__":
    main()

