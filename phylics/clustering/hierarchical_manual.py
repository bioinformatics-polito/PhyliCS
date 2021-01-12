#!/usr/bin/envs python 
import sys, os
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import average, complete, single, ward, fcluster
from cophenetic_affinity import cophenetic_affinity_score
from sklearn.metrics import silhouette_score
#scipy.cluster.hierarchy.linkage(y, method='single', metric='euclidean', optimal_ordering=False)

def plot_dendrogram(linkage_matrix, X, labels, chr_limits, chr_boundaries, chrN_list, pos_list, out_prefix):

    # Cluster colors
    color_palette = sns.color_palette("hls", len(np.unique(labels)))
    cluster_lut = dict(zip(np.unique(labels), color_palette))
    cluster_colors=pd.Series(labels, index=X.index).map(cluster_lut)

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
    plt.savefig(out_prefix + '_clusters_heatmap_pca.png')
    plt.clf()

def main():
    parser = argparse.ArgumentParser(description="Hierarchical clustering")
    parser.add_argument("-f", "--file", required=True, help="CNV file")
    parser.add_argument("-o", "--out_prefix", required=True, help="Output prefix")
    parser.add_argument("-m", "--matrix", help="Cophenetic distance file.", required=True, type=str)

    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, index_col=0, sep="\t")

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

    k = 10 #from the observation of the original heatmap

    # Low variance feature filtering 
    selector = VarianceThreshold() #filters out features with 0 variance from cell to cell (not informative for clustering)
    cnvs_filtered = pd.DataFrame(data=selector.fit_transform(cnvs), index=cnvs.index)


    # PCA (keep only significant features)
    #n_samples = len(df)
    #n_features = len(df.columns)
    pca = PCA(n_components=31)
    pca_result = pd.DataFrame(data=pca.fit_transform(cnvs_filtered.values), index=cnvs.index)

    Z_average = (average(pdist(pca_result, metric='cityblock')))
    labels_avg = fcluster(Z_average, k, 'maxclust')
    Z_complete = (complete(pdist(pca_result, metric='cityblock')))
    labels_complete = fcluster(Z_complete, k, 'maxclust')
    Z_single = (single(pdist(pca_result, metric='cityblock')))
    labels_single = fcluster(Z_single, k, 'maxclust')
    Z_ward = (ward(pdist(pca_result)))
    labels_ward = fcluster(Z_ward, k, 'maxclust')

    #plot_dendrogram(Z_average, cnvs, labels_avg, chr_limits, chr_boundaries, chrN_list, pos_list, args.out_prefix + "_average")
    #plot_dendrogram(Z_complete, cnvs, labels_complete, chr_limits, chr_boundaries, chrN_list, pos_list, args.out_prefix + "_complete")
    #plot_dendrogram(Z_single, cnvs, labels_single, chr_limits, chr_boundaries, chrN_list, pos_list, args.out_prefix + "_single")
    #plot_dendrogram(Z_ward, cnvs, labels_ward, chr_limits, chr_boundaries, chrN_list, pos_list, args.out_prefix + "_ward")

    print("Average linkage")
    S_avg = cophenetic_affinity_score(pd.DataFrame({'cell':cnvs.index, 'cluster':labels_avg}), np.unique(labels_avg), matrix)
    print("Complete linkage")
    S_complete = cophenetic_affinity_score(pd.DataFrame({'cell':cnvs.index, 'cluster':labels_complete}), np.unique(labels_complete), matrix)
    print("Single linkage")
    S_single = cophenetic_affinity_score(pd.DataFrame({'cell':cnvs.index, 'cluster':labels_single}), np.unique(labels_single), matrix)
    print("Ward linkage")
    S_ward = cophenetic_affinity_score(pd.DataFrame({'cell':cnvs.index, 'cluster':labels_ward}), np.unique(labels_ward), matrix)

    print("Coph. aff. score - agglomerative (average): {} ".format(S_avg))
    print("Coph. aff. score - agglomerative (complete): {} ".format(S_complete))
    print("Coph. aff. score - agglomerative (single): {} ".format(S_single))
    print("Coph. aff. score - agglomerative (ward): {} ".format(S_ward))

    """
    sil_avg =  silhouette_score(cnvs.values, labels_avg)
    sil_complete =  silhouette_score(cnvs.values, labels_complete)
    sil_single =  silhouette_score(cnvs.values, labels_single)
    sil_ward =  silhouette_score(cnvs.values, labels_ward)
    print("\n")
    print("Sil. score - agglomerative (average): {} ".format(sil_avg))
    print("Sil. score - agglomerative (complete): {} ".format(sil_complete))
    print("Sil. score - agglomerative (single): {} ".format(sil_single))
    print("Sil. score - agglomerative (ward): {} ".format(sil_ward))
    """
if __name__ == "__main__":
    sys.exit(main())