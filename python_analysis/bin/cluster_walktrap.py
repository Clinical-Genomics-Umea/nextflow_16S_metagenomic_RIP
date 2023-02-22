import numpy as np
import sys
import pandas as pd
from pathlib import Path
from datetime import date
import igraph as ig
import matplotlib.pyplot as plt


'''
Cluster adjacency matrix using the walktrap algorithm (Pons and Latapy, Computing Communities in Large Networks using random walks) 
Put Pearson corr cutoff as sys.argv[2] and output file/figure prefix name as sys.argv[3] 
'''

def make_adj_matrix(infile, cutoff, out_prefix):
    ''' Make an adjacency matrix from a pearson correlation matrix and output csv with ASV names of nodes '''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    indata = indata.set_index('Unnamed: 0')
    indata = indata.dropna(how='all', axis=0).dropna(how='all', axis=1)
    np.fill_diagonal(indata.values, 0)
    adj = pd.DataFrame((np.where(indata >= float(str(cutoff)), 1, indata)))
    adj = pd.DataFrame((np.where(adj <= float(str(cutoff)), 0, adj)))
    adj.index = indata.index
    bacs = pd.DataFrame(adj.index)
    bacs.columns = ['ASV']
    bacs.to_csv(str(out_prefix) + '_list_' + date.today().strftime('%y%m%d') + '.csv')
    adj.reset_index(inplace=True)
    adj = adj.drop(columns='Unnamed: 0')
    return adj


def cluster_adj_matrix(adj_matrix, out_prefix):
    ''' Cluster adjacency matrix and print .png'''
    graph = ig.Graph.Adjacency(adj_matrix, mode='undirected') # Make graph from adjacency matrix
    v_dendro = ig.Graph.community_walktrap(graph) # Create vertex dendrogram
    communities =  v_dendro.as_clustering() # Convert into VertexClustering for plotting

    #Print community members to file
    with open(str(out_prefix) + '_communities_' + date.today().strftime('%y%m%d') + '.txt', 'w') as outfile:
        for i, community in enumerate(communities):
            outfile.write(f"Community {i}:" + '\n')
            for v in community:
                outfile.write(f"\t{graph.vs[v]['name']}" + '\n')
    outfile.close()                
    
    graph.vs['label'] = graph.vs['name'] # Set labels from names

    # Set community colors
    num_communities = len(communities)
    palette1 = ig.RainbowPalette(n=num_communities)
    for i, community in enumerate(communities):
        graph.vs[community]["color"] = i
        community_edges = graph.es.select(_within=community)
        community_edges["color"] = i

    # Plot the communities
    fig1, ax1 = plt.subplots()
    ax1.set_title(str(out_prefix).split('/')[-1], fontsize=20)

    ig.plot(
        communities,
        target=ax1,
        mark_groups=True,
        palette=palette1,
        vertex_size=0.1,
        edge_width=0.5
    )
    fig1.set_size_inches(20, 20)
    fig1.savefig(str(out_prefix) + '_' + date.today().strftime('%y%m%d') + '.png', dpi=100)


def main():
    infile = Path(sys.argv[1])
    cutoff = Path(sys.argv[2])
    out_prefix = Path(sys.argv[3])  
    adj_matrix = make_adj_matrix(infile, cutoff, out_prefix)
    cluster_adj_matrix(adj_matrix, out_prefix)


if __name__ == "__main__":
    main()        