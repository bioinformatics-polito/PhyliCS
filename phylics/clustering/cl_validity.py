#!/usr/bin/env python
from math import nan
import os, sys
import argparse
import pandas as pd
import numpy as np
from hdbscan.validity import validity_index
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

def connectivity(points, labels):
    """
        Measures how often nearest neighbors are put in the same cluster.
        The connectivity has a value between zero and ∞ and should be minimized.

        C = sum_i sum_j (x_i,nn_ij)

        i = current point
        nnij = jth nearest neighbor of i
        x_i,nnij = 0, if they are in the same cluster,
                 = 1/j, if they are not (bigger values for nearest points)
        L = number of nearest neighbors
        
    """
    connectivity_score = 0.0
    N = len(points)
    L = round(N/2)
    neigh = NearestNeighbors(n_neighbors=L)
    nbrs = neigh.fit(points)
    distances, neighbors = nbrs.kneighbors(points) #indices stores a list of tuples of L nearest neighbors
    for neighbors_group in neighbors:
        connectivity_score_part = 0.0
        Ci = labels[neighbors_group[0]] #cluster of current cell
        for j in range(1, L): #nearest neighbor to the Lth nearest neighbor
            Cj = labels[neighbors_group[j]]
            if Ci == Cj:
                xij = 0
            else:
                xij = 1/j
            connectivity_score_part = connectivity_score_part + xij
        
        connectivity_score = connectivity_score + connectivity_score_part
    return connectivity_score
    

def delta_fast(ck, cl, distances):
    values = distances[np.where(ck)][:, np.where(cl)]
    values = values[np.nonzero(values)]

    return np.min(values)
    
def big_delta_fast(ci, distances):
    values = distances[np.where(ci)][:, np.where(ci)]
    #values = values[np.nonzero(values)]
            
    return np.max(values)

def dunn_fast(points, labels):
    #https://github.com/jqmviegas/jqm_cvi/blob/master/jqmcvi/base.py

    """ Dunn index - FAST (using sklearn pairwise euclidean_distance function)

    The Dunn Index is the ratio of the smallest distance between observations
    not in the same cluster to the largest intra-cluster distance.
    The Dunn Index has a value between zero and ∞, and should be
    maximized.
    
    Parameters
    ----------
    points : np.array
        np.array([N, p]) of all points
    labels: np.array
        np.array([N]) labels of all points
    """
    distances = euclidean_distances(points)
    ks = np.sort(np.unique(labels))
    
    deltas = np.ones([len(ks), len(ks)])*1000000
    big_deltas = np.zeros([len(ks), 1])
    
    l_range = list(range(0, len(ks)))
    
    for k in l_range:
        for l in (l_range[0:k]+l_range[k+1:]):
            deltas[k, l] = delta_fast((labels == ks[k]), (labels == ks[l]), distances)
        
        big_deltas[k] = big_delta_fast((labels == ks[k]), distances)

    di = np.min(deltas)/np.max(big_deltas)
    return di
    

def main():
    parser = argparse.ArgumentParser(description="Internal cluster validation.")

    parser.add_argument("-c", "--clusters_file", help="Clusters file", required=True)
    parser.add_argument("-f", "--file", help="Raw data file.", required=True, type=str)
    parser.add_argument("-s", "--stats_file", help="Stats file", required=True)
    parser.add_argument('-p', '--preproc', action='store_true', help="To be specified if the input file requires preprocessing")
    parser.add_argument("-o", "--out_prefix", help="Output prefix", required=True)

    args = parser.parse_args()

    clusters = pd.read_csv(args.clusters_file, sep="\t")
    clusters.columns = ['cell', 'cluster'] #rename columns, otherwise the cell column is named 'Unnamed: 0'

    if args.preproc:
        matrix = pd.read_csv(args.file, sep="\t", usecols = lambda column : column not in ['CHR', 'START', 'END']).transpose()
    else:
        matrix = pd.read_csv(args.file, index_col=0, sep="\t")

    stats = pd.read_csv(args.stats_file, index_col=0, sep="\t", header=None).transpose()
    #to fix an error
    stats = stats.loc[:,~stats.columns.duplicated()]

    labels_ = clusters['cluster'].values
    
    """
    Silhouette width

    The Silhouette value measures the degree of confidence in the clustering as-
    signment of a particular observation, with well-clustered observations having
    values near 1 and poorly clustered observations having values near −1.

    The Silhouette Width thus lies in the interval
    [−1, 1], and should be maximized.
    """
    try:
        S = silhouette_score(matrix.values, labels_)
    except:
        S = nan

    #Dunn Index

    try:
        DI = dunn_fast(matrix.values, labels_)
    except:
        DI=nan

    #Connectivity
    C = connectivity(matrix.values, labels_)

    # Fast Density Based Cluster Validity index
    # https://www.dbs.ifi.lmu.de/~zimek/publications/SDM2014/DBCV.pdf

    try:
        DBCV = validity_index(matrix.values, labels_)
    except:
        DBCV=nan
        """
        https://github.com/scikit-learn-contrib/hdbscan/issues/127
        
        Because I simply found that all the points of the given cluster are exactly 
        the same i.e. all pairwise distances of the points of this cluster are 0, and 
        then naturally 0**(-1/d) raises the error.
        """
    """
    print("silhouette_score\t" + str(S))
    print("dunn_index\t" + str(DI))
    print("connectivity\t" + str(C))
    print("density_based_cluster_validity\t" + str(DBCV))
    """
    with open(args.out_prefix + "_internal_validation_scores.csv", 'w+') as f:
        f.write("\n")
        f.write("silhouette_score\t" + str(S) + "\n")
        f.write("dunn_index\t" + str(DI) + "\n")
        f.write("connectivity\t" + str(C) + "\n")
        f.write("density_based_cluster_validity\t" + str(DBCV)+ "\n")
    
    
if __name__ == "__main__":
    sys.exit(main())