#!/usr/bin/env python

import hdbscan
import sys
import os 
import argparse
import operator
import pandas as pd 
import numpy as np
from itertools import combinations

from sklearn.metrics import silhouette_samples, silhouette_score, adjusted_rand_score

def dbscan_predict(model, X):

    nr_samples = X.shape[0]

    y_new = np.ones(shape=nr_samples, dtype=int) * -1

    for i in range(nr_samples):
        diff = model.components_ - X[i, :]  # NumPy broadcasting

        dist = np.linalg.norm(diff, axis=1)  # Euclidean distance

        shortest_dist_idx = np.argmin(dist)

        if dist[shortest_dist_idx] < model.eps:
            y_new[i] = model.labels_[model.core_sample_indices_[shortest_dist_idx]]

    return y_new

def check_nmin(nmin):
    class CheckValid(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            #print("args " + str(args))        
            if(len(values) < nmin):
                 msg='argument "{f}" requires at least {n} values'.format(f=self.dest, n=nmin)
                 raise argparse.ArgumentTypeError(msg)  
            setattr(args, self.dest, values)
    return CheckValid

def check_valid():
    class CheckValid(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            #print("args " + str(args))        
            if(len(values) != len(args.files)):
                 msg='argument "{f}" requires at least {n} values'.format(f=self.dest, n=len(args.files))
                 raise argparse.ArgumentTypeError(msg)  
            setattr(args, self.dest, values)
    return CheckValid

def ReadCnvFile(sample):
    cnvs = pd.read_csv(sample, sep="\t")
    cnvs = cnvs.drop(['CHR', 'START', 'END'], axis=1).transpose()
    return cnvs 

def ReadClusteringFile(clustering):
    clusters = pd.read_csv(clustering, sep="\t", index_col=0)
    return clusters

def cluster_size(matrix):
    model_dict = {}
    for min_size in range(3,50+1):
            model = hdbscan.HDBSCAN(min_cluster_size=min_size, min_samples=1, cluster_selection_epsilon=0.5).fit(matrix.values) 
            count = 0
            for p in model.probabilities_:
                if p < 0.05:
                    count = count +1
            model_dict[min_size] = count
    min_cluster_size=min(model_dict.items(), key=operator.itemgetter(1))[0]
    return min_cluster_size

def ClusterData(sample):
    #cluster_selection_method='eom', cluster_selection_epsilon=0.5, metric='l1'
    model = hdbscan.HDBSCAN(min_cluster_size=cluster_size(sample), min_samples=1, cluster_selection_epsilon=0.5, prediction_data=True).fit(sample.values)
    #clusters = pd.DataFrame(model.labels_, index=sample.index, columns=["cluster"])
    return model

def ClusteringOverlappingScore(sample1, sample2, clustering1, clustering2):
    #add sample2 cells to sample1 and compute clusters
    clusters, _ = hdbscan.approximate_predict(clustering1, sample2)
    ari1 = adjusted_rand_score(clustering2.labels_, clusters)
    print(np.unique(clusters))

    #add sample1 cells to sample2 and compute clusters
    clusters, _ = hdbscan.approximate_predict(clustering2, sample1)
    ari2 =  adjusted_rand_score(clustering1.labels_, clusters)
    print(np.unique(clusters))

    return (ari1 + ari2) / 2


def main():
    parser = argparse.ArgumentParser(description="Multi-sample clustering overlapping score.")

    parser.add_argument("files", help="List of input files",
                            type=argparse.FileType('r'), nargs='+', action=check_nmin(2))
    parser.add_argument("outprefix", help="Output path prefix", 
                            type=str)
    parser.add_argument("-c", "--clusterings", required=False, help="List of cluster files (in the same order of input files)",
                            type=argparse.FileType('r'), nargs='+', action=check_valid())
    parser.add_argument("-n", "--names", required=False, help="List of sample names (in the same order of input files)",
                            type=str, nargs='+', action=check_valid())
    #parser.add_argument("-ca", "--clustering_algorithm")

    #parser.add_argument("-p", "--preprocessing")

    

    args = parser.parse_args()

    # Read input data and store them in a data structure
    samples = {}

    if args.names is not None:
        names = args.names
    else:
        names = np.arange(len(args.files))

    for sample, name in zip(args.files, names):
        cnvs = ReadCnvFile(sample)
        samples[name] = cnvs

    #Step1: if not specified on the command line compute clusterings for each sample
    clusterings = {}
    if args.clusterings is not None:
        for clustering in args.clusterings:
            clusterings[name] = ReadClusteringFile(clustering)
    else:
        for name, sample in samples.items():
            clusterings[name] = ClusterData(sample)

    #Step2: for each pair of samples, compute the clustering overlapping score
    # 1. Add, one by one, each cell from sample2 to sample1 and recluster. Keep track of the cluster to which the cell has been assigned
    # 2. Compute silhouette score for the assignments of cells from sample2 to clusters of sample1
    # 3. Switch the samples and repeat points 1. and 2.
    # 4. Compute an average score: it reflects how similar the cluster structure of the two samples is
    clustering_overlapping_scores = {}
    for samples_pair in combinations(names, 2):
        clustering_overlapping_scores[samples_pair[0], samples_pair[1]] = ClusteringOverlappingScore(samples[samples_pair[0]], samples[samples_pair[1]], 
                                                        clusterings[samples_pair[0]], clusterings[samples_pair[1]])

    for names, score in clustering_overlapping_scores.items():
        print("{sample1} - {sample2}: {score}".format(sample1=names[0], sample2=names[1], score=score))

if __name__ == "__main__":
    sys.exit(main())