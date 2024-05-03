import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, MinMaxScaler


def norm_asvs(infile, outfile):
    ''' Take a Lulu output file and normalize using MinMax scaling'''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    indata = indata.T
    indata.columns = indata.iloc[0]
    indata = indata.drop(indata.index[0])
    sample_ids = pd.DataFrame(indata.index)
    clin = indata[['Age', 'Sex', 'Stage', 'Cancer']]
    indata = indata.drop(columns=['Age', 'Sex', 'Stage', 'Cancer'])
    asv_names = indata.columns 
    indata = indata.to_numpy()
    data = indata.astype(np.float)

    preprocessor = Pipeline([("scaler", MinMaxScaler())])                
    pipe = Pipeline([("preprocessor", preprocessor)])
    pipe.fit(data)
    pcadf = pd.DataFrame(pipe["preprocessor"].transform(data))
    
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
    norm_asvs("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/lulu_out_filtered_clin_noS9_221201.csv",
               "/home/lindak/project/nextflow_16S_metagenomic_RIP/data/norm_lulu_filtered_221222.csv")

main()              