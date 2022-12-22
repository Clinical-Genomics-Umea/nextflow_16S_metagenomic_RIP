import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler


indata = pd.DataFrame(pd.read_csv("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/lulu_out_filtered_clin_noS9_221201.csv",sep=',', header=0))
indata = indata.T
indata.columns = indata.iloc[0]
indata = indata.drop(indata.index[0])
indata = indata[indata.Cancer == 'C']
# Make numpy array
indata = indata.drop(columns=['Cancer', 'Sex', 'Age', 'Stage'])
indata = indata.to_numpy()
data = indata.astype(np.float)


kmeans_kwargs = {
    "init": "k-means++",
    "n_init": 50 ,
    "max_iter": 500,
    "random_state": 42,
}

# A list holds the SSE values for each k
sse = []
for k in range(1, 11):
    kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
    kmeans.fit(indata)
    sse.append(kmeans.inertia_)


plt.style.use("fivethirtyeight")
plt.figure(figsize=(10, 10))
plt.plot(range(1, 11), sse)
plt.xticks(range(1, 11))
plt.xlabel("Number of Clusters")
plt.ylabel("SSE")
plt.savefig("/home/lindak/project/nextflow_16S_metagenomic_RIP/plots/elbow_noS9_221201.png", dpi=1000)
 