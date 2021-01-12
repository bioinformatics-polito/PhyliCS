#!/usr/bin/env python
import os
import sys
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from matplotlib import cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def ReadData(segcopy):
    df = pd.read_csv(segcopy, sep="\t")
    cnvs = df.drop(['CHR', 'START', 'END'], axis=1).transpose()

    cells = cnvs.index.tolist()

    #remove "cell" prefix
    index_ = []
    for cell in cells:
        cellN = cell[4:]
        index_.append(int(cellN))

    cnvs.index = pd.Index(data=index_, dtype=int)

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
    return cnvs, chr_boundaries, chrN_list, chr_limits, pos_list

def ReadClusters(clusters_file):
    clusters = pd.read_csv(clusters_file, sep="\t", index_col=0)
    return clusters

def AssociateClusterData(cnvs, clusters):
    cnvs['cluster'] = clusters['ClusterNumber']
    return cnvs

def ClusterHeatmap(cnvs, chr_boundaries, chrN_list, chr_limits, pos_list, outprefix):
    cnvs = cnvs.sort_values(by='cluster')
    
    sorted_labels = cnvs['cluster'].values

    color_palette = sns.color_palette("hls", len(np.unique(sorted_labels)))
    cluster_colors = [color_palette[x-1] if x >= 0
                      else (0.5, 0.5, 0.5)
                      for x in sorted_labels]
    #print(row_colors)
    
    cnvs = cnvs.drop(['cluster'], axis=1)
    cbar_kws={"ticks":np.arange(0,7,1)}
    
    sns.set(font_scale=2)
    h = sns.clustermap(cnvs,
            row_cluster=False,
            col_cluster=False,
            yticklabels = False,
            row_colors=cluster_colors,
            cmap='RdBu_r',
            vmin=0, vmax=6,
            center = 2,
            #norm=divnorm,
            cbar_kws=cbar_kws)

    h.cax.set_position([0.2, .2, .03, .45])
    for label in np.unique(sorted_labels-1):
        h.ax_col_dendrogram.bar(0, 0, color=color_palette[label],label=label, linewidth=0)
    h.ax_col_dendrogram.legend(loc="center",  ncol=len(sorted_labels))

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

    ax.set_xlabel("chromosomes", fontsize=16, fontweight='bold')
    ax.set_ylabel("cells", fontsize=16, fontweight='bold')
       
    #plt.gcf().suptitle('Cnv profiles', fontsize=16, fontweight='bold' )
    plt.gcf().set_size_inches(37, 21)
    plt.savefig(os.path.join(outprefix, 'newick_clusters_heatmap.png'), bbox_inches = 'tight')
    plt.clf()

def main():
    parser = argparse.ArgumentParser(description="Draw newick cluster heatmap")
    parser.add_argument("-f", "--file", type=str, required=True, help="CNV data file.")
    parser.add_argument("-c", "--clustersfile", type=str, required=True, help="Clusters file.")
    parser.add_argument("-o", "--outprefix", type=str, required=True, help="Output prefix.")

    args = parser.parse_args()

    
    cnvs, chr_boundaries, chrN_list, chr_limits, pos_list = ReadData(args.file)
    clusters = ReadClusters(args.clustersfile)
    cnvs = AssociateClusterData(cnvs, clusters)
    ClusterHeatmap(cnvs, chr_boundaries, chrN_list, chr_limits, pos_list, args.outprefix)


if __name__ == "__main__":
    sys.exit(main())
