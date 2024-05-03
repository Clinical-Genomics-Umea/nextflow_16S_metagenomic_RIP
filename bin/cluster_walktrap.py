import numpy as np
import sys
import pandas as pd
from pathlib import Path
from datetime import date
import igraph as ig
import matplotlib.pyplot as plt
import random


'''
Cluster adjacency matrix using the walktrap algorithm (Pons and Latapy, Computing Communities in Large Networks using random walks) 
Sparsity cut off is 0.5% 
'''


def make_adj_matrix_spars(infile, out_prefix):
    ''' 
    Make an adjacency matrix from a pearson correlation matrix and output csv with ASV/sample names of nodes
    The generated adj. matrix contains 0.5% sparsity (the 0.5% highest correlations)      
    '''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    indata = indata.set_index('Unnamed: 0')
    indata = indata.dropna(how='all', axis=0).dropna(how='all', axis=1)
    np.fill_diagonal(indata.values, 0)
    out = indata.to_numpy()
    adj = out.copy()
    print(out.size - len(out))
    frac = (out.size - len(out)) * 0.005
    indx = np.argwhere(np.isin(out, np.sort(out, axis=None)[-(int(frac)):])) # get indices of top 0.5% values
    adj[indx[:, 0], indx[:, 1]] = 1
    inv_indx = np.argwhere(~(np.isin(out, np.sort(out, axis=None)[-(int(frac)):]))) # get all other indices
    adj[inv_indx[:, 0], inv_indx[:, 1]] = 0
    adj = pd.DataFrame(adj)
    adj.index = indata.index
    adj.columns = indata.columns
    print(adj.shape)
    adj = adj.loc[(adj!=0).any(axis=1)]
    adj = adj.loc[:,adj.any()]  
    adj.to_csv(str(out_prefix) + '_adj_matrix_' + date.today().strftime('%y%m%d') + '.csv')   
    bacs = pd.DataFrame(adj.index)
    bacs.columns = ['ID']
    bacs.to_csv(str(out_prefix) + '_list_' + date.today().strftime('%y%m%d') + '.csv')
    adj.reset_index(inplace=True)
    adj = adj.drop(columns='Unnamed: 0')
    adj.columns = adj.index
    print('Number of edges included: ' + str(int(adj.values.sum())))
    return adj
   


def cluster_adj_matrix(adj_matrix, out_prefix):
    ''' Cluster adjacency matrix and print .png'''
    random.seed(1111)
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
        vertex_label_size=15,
        edge_width=1
    )
    fig1.set_size_inches(20, 20)
    fig1.savefig(str(out_prefix) + '_' + date.today().strftime('%y%m%d') + '.png', dpi=100)


def main():
    infile = Path(sys.argv[1])
    out_prefix = Path(sys.argv[2])
    adj_matrix = make_adj_matrix_spars(infile, out_prefix)
    cluster_adj_matrix(adj_matrix, out_prefix)


if __name__ == "__main__":
    main()        