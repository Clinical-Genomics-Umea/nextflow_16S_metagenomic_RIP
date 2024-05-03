import pandas as pd
import os


''' Filter asv table from mannWU where p>0.05  '''

mann_data = pd.DataFrame(pd.read_csv("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/mannWU_norm_filtered_221222.csv", sep=',', header=0))
#mann_data = pd.DataFrame(pd.read_csv("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/mann_test.csv", sep=',', header=0))

bac_data = pd.DataFrame(pd.read_csv("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/norm_lulu_filtered_221222.csv", sep=',', header=0))
#bac_data = pd.DataFrame(pd.read_csv("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/in_test.csv", sep=',', header=0))
bac_data = bac_data.T
bac_data.columns = bac_data.iloc[0]
bac_data = bac_data[1:]
bac_data = bac_data[bac_data.Cancer == 'N']
bac_data = bac_data.T
bac_data = bac_data.reset_index()

norm_data = mann_data[mann_data['p_value'] > 0.05]
norm_data = norm_data[['ASV']]

norm_data = norm_data.merge(bac_data, left_on='ASV', right_on='Unnamed: 0')
norm_data = norm_data.drop(['Unnamed: 0'], axis=1)
norm_data.set_index('ASV', inplace=True)
norm_data = norm_data.astype('float')

# Retain only ASVs where no normal sample = 0
#out  = norm_data[(norm_data != 0).all(axis=1)]
#print(out.index)


## Retain only ASVs where at most 5 normal samples of 95 to be 0 at some ASVs
mask = (norm_data == 0).sum(1) >= len(norm_data.columns) / 20 
out = norm_data[~mask]
print(out.index)

lulu_data = pd.DataFrame(pd.read_csv("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/norm_lulu_filtered_221222.csv", sep=',', header=0))
lulu_data.set_index('Unnamed: 0', inplace=True)

plot_data = lulu_data.drop(out.index)
plot_data = plot_data.reset_index()
plot_data.to_csv("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/norm_lulu_filtered_clin_noS9_NoNormASVs_221222.csv", index=False)


