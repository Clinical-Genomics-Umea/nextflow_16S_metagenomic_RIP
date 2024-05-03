import pandas as pd

indata = pd.DataFrame(pd.read_csv("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/lulu_out_filtered_clin_221129.csv",sep=',', header=0))
indata = indata.T
indata.columns = indata.iloc[0]
indata = indata.drop(indata.index[0])
out = indata[indata.Stage != '9']
out = out.T
out.index.name = None
out.to_csv("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/lulu_out_filtered_clin_noS9_221201.csv")