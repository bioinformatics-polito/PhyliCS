#!/usr/bin/envs python

import sys
import argparse
import numpy as np 
import pandas as pd 

def check_number_of_labels(n_labels, n_samples):
    """Check that number of labels are valid.
    Parameters
    ----------
    n_labels : int
        Number of labels
    n_samples : int
        Number of samples
    """
    if not 1 < n_labels < n_samples:
        raise ValueError("Number of labels is %d. Valid values are 2 "
                         "to n_samples - 1 (inclusive)" % n_labels)

def a(x, coph_dist_matrix, n, cells):
    # Mean intra-cluster distance
    # a(x) =  1 / (ni - 1) sum_y_belonging_Ci,y_different_from_x( coph_d(x, y) ) 
    d = 0.0
    for y in cells:
        if y != x:
            coph_d = coph_dist_matrix.loc[x,y]
            d = d + coph_d
    return (d / (n - 1))

def b(x, coph_dist_matrix, Ci, labels, clusters):
    # Min mean inter-cluster distance, with respect to current cluster Ci
    # b(x) = min_j_different_from_i [ (1 / nj) sum_j_belonging_Cj( coph_d(x, y) ) ]  
    D = []
    for Cj in labels:
        if Cj != Ci:
            cells = clusters['cell'][clusters['cluster'] == Cj]
            n = len(cells)
            d = 0.0 #mean distance between x and all points in cluster Cj
            for y in cells:
                coph_d = coph_dist_matrix.loc[x,y]
                d = d + coph_d
            D.append(d / n)
    return min(D)

def cophenetic_affinity_score(clusters, labels, coph_dist_matrix):
    # Parameters: 
    #   clusters: clusters file
    #   coph_dist_matrix: pairwise cophenetic distance matrix, computed on a dendrogram
    # Return values:
    #    cophenetic affinity score
    # 
    # This score measures the level of affinity between the cluster assignment computed by 
    # a clustering algorithm and the dispersion of the same data within the evolutionary tree
    # which generated them. 
    #
    # Mathematical formulation:
    # cophenetic_affinity_score (C.A.S) = 1/NC * sum_i( 1/n_i * sum_x_belonging_C_i( (b(x) - a(x)) / max[ b(x), a(x) ] ) )
    #
    # a(x)  = The mean cophenetic distance between a sample and all other points in the same class, 
    #           computed on the dendrogram
    # b(x)  = The mean cophenetic distance between a sample and all other points in the next nearest cluster,
    #           computed on the dendrogram
    # NC    = Clusters number 
    # C_i   = Each cluster
    # 
    # Compute average Cluster Cophenetic Affinity score
    """
    cluster_cardinality = clusters.groupby('cluster').count() #idx = cluster, column = number of cells
    singletons = cluster_cardinality[cluster_cardinality['cell'] == 1].index.values
    n_singletons = len(singletons)
    #print("singletons_n: " + str(singletons_n))
    
    clusters = clusters[~clusters['cluster'].isin(singletons)]
    labels = clusters['cluster'].unique()
    """
    NC = len(labels)
    n_samples = len(clusters)
    check_number_of_labels(NC, n_samples)
    S = 0.0
    for C in labels:
        # compute single cluster coph aff score
        cells = clusters['cell'][clusters['cluster'] == C] 
        s = 0.0
        n = len(cells)
        for x in cells:
            a_x = a(x, coph_dist_matrix, n, cells)
            b_x = b(x, coph_dist_matrix, C, labels, clusters)
            s = s + ((b_x - a_x)/max(b_x, a_x))
        s = s  / n
        #print("Cluster {} - score: {} (cluster size = {})".format(C, s, n))
        S = S + s*n/n_samples
    #S = S/NC #cophenetic affinity score
    return S  

def compute_pvalue(clusters, coph_dist_matrix, obsv_score):
    n_perm = 1000
    perm_scores = []
    for i in range(n_perm):
        clusters_perm = pd.DataFrame(columns=['cell', 'cluster'])
        clusters_perm['cell'] = clusters['cell']
        clusters_perm['cluster'] = np.random.permutation(clusters['cluster'].values)
        labels_perm = clusters_perm['cluster'].unique()
        s = cophenetic_affinity_score(clusters_perm, labels_perm, coph_dist_matrix)
        perm_scores.append(s)
    p_val = 0.0 
    for s in perm_scores:
        if s > obsv_score:
            p_val = p_val + 1
    return (p_val / n_perm)
          
def main():
    parser = argparse.ArgumentParser(description="Cophenetic Affinity Score.")

    parser.add_argument("-f", "--file", help="Clusters file.", required=True, type=str)
    #parser.add_argument("-s", "--stats", help="Stats file.", required=True, type=str)
    parser.add_argument("-m", "--matrix", help="Cophenetic distance file.", required=True, type=str)
    parser.add_argument("-o", "--outprefix", help="Output prefix path.", required=True, type=str)

    args = parser.parse_args()

    clusters = pd.read_csv(args.file, sep="\t")
    matrix = pd.read_csv(args.matrix, index_col=0, sep="\t")
    #outprefix = args.outprefix

    clusters.columns = ['cell', 'cluster'] #rename columns, otherwise the cell column is named 'Unnamed: 0'

    # filter out cells labeled as outliers (dbscan, hdbscan) and singleton clusters

    clusters = clusters[clusters['cluster'] != -1]
    n_outliers = len(clusters[clusters['cluster'] == -1])   
    #print("outliers_n: " + str(outliers_n))

    labels = clusters['cluster'].unique()
    n_labels_unfilt = len(labels)
    #print("unfiltered_clusters_n: " + str(len(labels)))

    
    cluster_cardinality = clusters.groupby('cluster').count() #idx = cluster, column = number of cells
    singletons = cluster_cardinality[cluster_cardinality['cell'] == 1].index.values
    n_singletons = len(singletons)
    #print("singletons_n: " + str(singletons_n))
    
    clusters = clusters[~clusters['cluster'].isin(singletons)]
    labels = clusters['cluster'].unique() 
    n_labels_filt = len(labels)
    #print("filtered_clusters_n: " + str(len(labels)))

    try:
        S = cophenetic_affinity_score(clusters, labels, matrix)
        #pvalue = compute_pvalue(clusters, matrix, S)
        #print("cophenetic_affinity_score: " + str(S))
    except ValueError:
        S = 'NaN'
        #pvalue = 0.0
       #print("Too fiew clusters.")
    
    with open(args.outprefix + "_coph_aff_score.csv", 'w+') as f:
        """
        f.write("\n")
        f.write("n_outliers\t" + str(n_outliers) + "\n")
        f.write("n_clusters_unfiltered\t" + str(n_labels_unfilt)+ "\n")
        f.write("n_singletons\t" + str(n_singletons)+ "\n")
        f.write("n_clusters_filtered\t" + str(n_labels_filt)+ "\n")
        """
        f.write("cophenetic_affinity_score\t" + str(S)+ "\n")
    
    #print("Cophenetic_affinity_score: " + str(S))
if __name__ == "__main__":
        sys.exit(main())
  
