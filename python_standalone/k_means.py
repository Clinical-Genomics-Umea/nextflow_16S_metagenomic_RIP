import tarfile
import urllib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, MinMaxScaler


indata = pd.DataFrame(pd.read_csv("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/lulu_out_filtered_clin_noU274.csv",sep=',', header=0))
indata = indata.T
indata.columns = indata.iloc[0]
indata = indata.drop(indata.index[0])

# Only cancer samples
########################
indata = indata[indata.Cancer == 'C']

# Get the truth labels
inlabels = indata['Stadium']
inlabels = pd.DataFrame(inlabels, index=indata.index, copy=True)
inlabels = inlabels['Stadium'].to_numpy()

####################################################

## Cancer and Normal samples
#########################

# # Get the truth labels
# inlabels = indata['Stadium']
# inlabels = pd.DataFrame(inlabels, index=indata.index, copy=True)
# inlabels = inlabels['Stadium'].to_numpy()
# inlabels = np.array(inlabels).astype(str)
# inlabels[inlabels == 'nan'] = 'N'
# print(inlabels)

####################################################


# Make numpy array
indata = indata.drop(columns=['Cancer', 'Sex', 'Age', 'Lokal', 'Stadium'])
indata = indata.to_numpy()
data = indata.astype(np.float)

# Make labels into ints
label_encoder = LabelEncoder()
true_labels = label_encoder.fit_transform(inlabels)
n_clusters = len(label_encoder.classes_)

# Scale and do dim reduction
preprocessor = Pipeline([("scaler", MinMaxScaler()), # MinmaxScaler when can not assume norm distribution
                ("pca", PCA(n_components=2, random_state=42))])

clusterer = Pipeline([("kmeans", KMeans(
                n_clusters=3,
                init="k-means++",
                n_init=50,
                max_iter=500,
                random_state=42))])

pipe = Pipeline([("preprocessor", preprocessor), ("clusterer", clusterer)])
pipe.fit(data)


preprocessed_data = pipe["preprocessor"].transform(data)
pred_labels = pipe["clusterer"]["kmeans"].labels_
#print("Silhouette score: " , silhouette_score(preprocessed_data, pred_labels))
#print("ARI score: " , adjusted_rand_score(true_labels, pred_labels))

## Plot
pcadf = pd.DataFrame(pipe["preprocessor"].transform(data), columns=["component_1", "component_2"])
pcadf["pred"] = pipe["clusterer"]["kmeans"].labels_
pcadf["true"] = label_encoder.inverse_transform(true_labels)


plt.style.use("fivethirtyeight")
plt.figure(figsize=(10, 10))

scat = sns.scatterplot(
    x="component_1",
    y="component_2",
    s=50,
    data=pcadf,
    hue="pred",
    style="true",
    palette="Set2",
)

scat.set_title("K-means clustering no_U274 n_clust=3")
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
plt.legend(loc= "upper right",bbox_to_anchor=(1, 1))
plt.savefig("/home/lindak/project/nextflow_16S_metagenomic_RIP/plots/k_means_clust_noU274_n_clust3.png", dpi=1000)
#plt.show()
