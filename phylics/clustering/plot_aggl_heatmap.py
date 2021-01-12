#!/usr/bin/envs python

"""
=========================================
Plot Hierarachical Clustering Dendrogram 
=========================================
This example plots the corresponding dendrogram of a hierarchical clustering
using AgglomerativeClustering and the dendrogram method available in scipy.
"""

import sys, os
import argparse

import numpy as np
import pandas as pd

import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering

def plot_dendrogram(X, model, chr_limits, chr_boundaries, chrN_list, pos_list, out_prefix):

    # Children of hierarchical clustering
    children = model.children_

    # Distances between each pair of children
    # Since we don't have this information, we can use a uniform one for plotting
    distance = np.arange(children.shape[0])

    # The number of observations contained in each cluster level
    no_of_observations = np.arange(2, children.shape[0]+2)

    # Create linkage matrix and then plot the dendrogram
    linkage_matrix = np.column_stack([children, distance, no_of_observations]).astype(float)

    # Plot the corresponding dendrogram
    #dendrogram(linkage_matrix, labels=labels)
    #plt.savefig("/home/bioeda/bio/spatial/snakess/dataset/synth/dendrogram.png")

    # Cluster colors
    color_palette = sns.color_palette("hls", len(np.unique(model.labels_)))
    cluster_lut = dict(zip(np.unique(model.labels_), color_palette))
    cluster_colors=pd.Series(model.labels_, index=X.index).map(cluster_lut)

    #Plot seaborn clustermap
    #cbar_kws={"ticks":np.arange(0,7,1)}
    h = sns.clustermap(X, 
                    col_cluster=False, 
                    row_linkage=linkage_matrix,
                    yticklabels = False,
                    row_colors=cluster_colors,
                    cmap='RdBu_r',
                    vmin=0, vmax=6,
                    center = 2,
                    #norm=divnorm,
                    #cbar_kws=cbar_kws
                    )
    #h.cax.set_position([0.05, .2, .03, .45])
    ax = h.ax_heatmap
    #place vertical lines to identify chromosomes
    for pos in chr_limits:
            ax.axvline(x=pos, color='black')

    #place chromosome ticks at the right position
    ax.xaxis.set_major_locator(ticker.FixedLocator(pos_list))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter((chrN_list)))
    ax.tick_params(axis='x', rotation=0, labelsize=14, which='major')

    ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_boundaries))
    ax.tick_params(axis='x', length=20, which='minor')

    ax.set_xlabel("chromosomes", fontsize=14, fontweight='bold')
    ax.set_ylabel("cells", fontsize=14, fontweight='bold')
       
    plt.gcf().suptitle('Cnv profiles', fontsize=16, fontweight='bold' )
    plt.gcf().set_size_inches(37, 21)
    plt.savefig(out_prefix + '_clusters_heatmap.png')
    plt.clf()


def main():
    parser = argparse.ArgumentParser(description="Plot agglomerative clustering heatmap and dendrogram.")

    parser.add_argument("-s", "--stats", required=True, help="Stats file")
    parser.add_argument("-f", "--file", required=True, help="CNV file")
    parser.add_argument("-o", "--out_prefix", required=True, help="Output prefix")

    parser.add_argument("-l", "--linkage", required=True, choices=["average", "complete", "single", "ward"])

    args = parser.parse_args()

    stats = pd.read_csv(args.stats, index_col=0, sep="\t", header=None).transpose()
    stats = stats.loc[:,~stats.columns.duplicated()]
    k = int(stats['n_clusters_unfiltered'].values[0])
    
    df = pd.read_csv(args.file, sep="\t")
    cnvs = df.drop(['CHR', 'START', 'END'], axis=1).transpose()
    boundaries = df[['CHR', 'START', 'END']].copy()
    chr_limits = boundaries.index[boundaries['END'].isin(boundaries.groupby('CHR', sort=False)['END'].max().values)].tolist()
    chr_boundaries = np.append(0, chr_limits)
    chr_list = boundaries['CHR'].unique().tolist()
    chrN_list = []

    for x in chr_list:
        x = x[3:] #remove 'chr' for readability
        chrN_list.append(x)

    #compute the position where chromosome labels will be placed on the plots
    start = 0
    pos_list = []
    for end in chr_limits:
        pos_list.append((start+end)/2)
        start = end+1

    if args.linkage == "ward":
        affinity = 'euclidean'
    else:
        affinity = 'l1'

    model = AgglomerativeClustering(n_clusters=k, affinity=affinity, linkage=args.linkage)

    model = model.fit(cnvs.values)
    plot_dendrogram(cnvs, model, chr_limits, chr_boundaries, chrN_list, pos_list, args.out_prefix)

if __name__ == "__main__":
    sys.exit(main())