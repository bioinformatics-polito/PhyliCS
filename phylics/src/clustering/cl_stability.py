#!/usr/bin/env python
import os, sys
import random
import argparse
import pandas as pd
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import Birch
from sklearn.cluster import DBSCAN
import hdbscan
from sklearn.base import ClusterMixin
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from kneed import KneeLocator
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

# Wrapper function for Gaussian Mixture Modelling
class GMClusters(GaussianMixture, ClusterMixin):

    def __init__(self, n_clusters=2, **kwargs):
        kwargs["n_components"] = n_clusters
        super(GMClusters, self).__init__(**kwargs)

    def fit(self, X):
        super(GMClusters, self).fit(X)
        self.labels_ = self.predict(X)
        return self 

def APN(X, clusters, K, model, n_iter=None):
    # Average Proportion of Non-overlap (APN)
    """
        Measures the proportion of observations not placed in the same cluster by 
        clustering based on the full data and clustering based on the
        data with a single column removed.

        APN(K) = (1 / MN) sum_i( sum_j( 1 - n(C_il intersect C_i0)/n(Ci0) ) )
    
        C_i0 = obs_cluster containing observation i
        C_il = cluster_without_col_l containing observation i

        M = number of columns
        N = number of observations
        n_Ci0 = Ci0 size

        0 <= APN <= 1, with optimal value = 0 
    """
    if n_iter is not None:
        # Select randomly n_iter columns to be removed in the score computation
        # to limit the execution time
        columns = random.sample(range(max(X.columns.values)), n_iter)
        M = len(columns)
    else:
        columns = X.columns.values
        M = len(columns)

    cells = X.index.values
    N = len(cells)
    # for each dropped column

    apn = 0.0
    for l in columns:
        df = X.drop(l, axis=1)
        clusterer =  model.fit(df.values)
        # recompute clusters
        new_clusters = pd.DataFrame({"cell" : df.index.values, "cluster" : clusterer.labels_})
        # for each cell
        apn_l = 0.0
        for c in cells:
            #compute the proportion of non-overlapping between the original cluster (C0)
            # containing c and the new one (Cl)
            label_0 = clusters['cluster'][clusters['cell'] == c].values[0]
            label_l = new_clusters['cluster'][new_clusters['cell'] == c].values[0]

            C0 = clusters['cell'][clusters['cluster'] == label_0]
            Cl = new_clusters['cell'][new_clusters['cluster'] == label_l]
            overlap_proportion = len(set(C0).intersection(set(Cl))) / len(C0)

            apn_l = apn_l + (1 - overlap_proportion) 
        apn = apn + apn_l / N     
    return ( apn / M) 

def AD(X, clusters, K, clust_method, **kwargs):
    # Average Distance (AD)
    return ad 

def ADM(X, clusters, K, clust_method, **kwargs):
    # Average Distance between Means (ADM)
    return adm

def FOM(X, clusters, K, clust_method, **kwargs):
    # Figure of Merit (FOM)
    return fom 

def main():
    parser = argparse.ArgumentParser(description="Cluster Stability Estimator.")

    parser.add_argument("-c", "--clusters_file", help="Clusters file", required=True)
    parser.add_argument("-f", "--file", help="Raw data file.", required=True, type=str)
    parser.add_argument("-s", "--stats_file", help="Stats file", required=True)
    parser.add_argument("-m", "--method", 
            choices=["affinity", "agglomerative_average", "agglomerative_complete", "agglomerative_single", 
            "agglomerative_ward", "birch", "dbscan", "gaussian", "hdbscan", "kmeans"], 
            help="Clustering method", required=True)
    parser.add_argument('-p', '--preproc', action='store_true', help="To be specified if the input file requires preprocessing")
    parser.add_argument('-n', '--n_iter', help="Iteration number limit (meta-heuristic to limit the execution time in case of a high number of dimensions).", 
            type=int, default=None)
    #parser.add_argument("-o", "--out_prefix", help="Output prefix", required=True)

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

    #outprefix = args.out_prefix

    method = args.method

    k = int(stats["n_clusters_unfiltered"].values[0])

    # seed the rand generator for reproducibility of results
    random.seed(42)

    if "DBSCAN_eps" in stats.columns.values:
        eps = float(stats["DBSCAN_eps"].values[0])
    if "HDBSCAN_min_cluster_size" in stats.columns.values:
        min_cluster_size = int(stats["HDBSCAN_min_cluster_size"].values[0])
  
    if method == "affinity":
        model = AffinityPropagation()
        apn = APN(matrix, clusters, k, model, args.n_iter)
    elif method == "agglomerative_average":
        model = AgglomerativeClustering(n_clusters=k, affinity="l1", linkage="average")
        apn = APN(matrix, clusters, k, model, args.n_iter)
    elif method == "agglomerative_complete":
        model = AgglomerativeClustering(n_clusters=k, affinity="l1", linkage="complete")
        apn = APN(matrix, clusters, k, model, args.n_iter)
    elif method == "agglomerative_single":
        model = AgglomerativeClustering(n_clusters=k, affinity="l1", linkage="single")
        apn = APN(matrix, clusters, k, model, args.n_iter)
    elif method == "agglomerative_ward":
        model = AgglomerativeClustering(n_clusters=k, linkage="ward")
        apn = APN(matrix, clusters, k, model, args.n_iter)
    elif method == "birch":
        model = Birch(branching_factor=20, n_clusters=k, threshold=0.5, compute_labels=True)
        apn = APN(matrix, clusters, k, model, args.n_iter)
    elif method == "dbscan":
        model = DBSCAN(eps=eps)
        apn = APN(matrix, clusters, k, model, args.n_iter)
    elif method == "hdbscan":
        model = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=1, cluster_selection_epsilon=0.5)
        apn = APN(matrix, clusters, k, model, args.n_iter)
    elif method == "gaussian":
        model =  GMClusters(n_components=k)
        apn = APN(matrix, clusters, k, model, args.n_iter)
    elif method == "kmeans":
        model = KMeans(n_clusters=k)
        apn = APN(matrix, clusters, k, model, args.n_iter)

    
    with open(args.stats_file, 'a') as f:
        f.write("\n")
        f.write("APN\t" + str(apn) + "\n")
    
if __name__ == "__main__":
    sys.exit(main())