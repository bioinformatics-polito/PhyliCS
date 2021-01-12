import sys, os
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def main():
    parser = argparse.ArgumentParser(description="Plot cluster heatmap")

    parser.add_argument("-c", "--clusters", required=True, help="Clusters file")
    parser.add_argument("-f", "--file", required=True, help="CNV file")
    parser.add_argument("-o", "--out_prefix", required=True, help="Output prefix")

    args = parser.parse_args()

    clusters = pd.read_csv(args.clusters, sep="\t")
    clusters.columns = ['cell', 'cluster'] 

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
    
    cnvs['cluster'] = clusters['cluster'].values
    cnvs = cnvs.sort_values(by='cluster')
    sorted_labels = cnvs['cluster'].values

    color_palette = sns.color_palette("hls", len(np.unique(cnvs['cluster'].values)))
    cluster_colors = [color_palette[x] if x >= 0
                      else (0.5, 0.5, 0.5)
                      for x in sorted_labels]
   
    cnvs = cnvs.drop(['cluster'], axis=1)
    cbar_kws={"ticks":np.arange(0,7,1)}
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
    h.cax.set_position([0.05, .2, .03, .45])
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
    plt.savefig(args.out_prefix + '_clusters_heatmap.png')
    plt.clf()


if __name__ == "__main__":
    sys.exit(main())