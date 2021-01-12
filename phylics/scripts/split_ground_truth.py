#!/usr/bin/env python
import sys
import os
import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as ticker
import matplotlib.colors as colors
import numpy as np
import argparse
import seaborn as sns


def Heatmap(segcopy, outprefix):
    cnvs = segcopy.drop(['CHR', 'START', 'END'], axis=1).transpose()
    boundaries = segcopy[['CHR', 'START', 'END']].copy()
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
   
    cbar_kws={"ticks":np.arange(0,7,1)}
    
    sns.set(font_scale=2)
    h = sns.clustermap(cnvs,
            row_cluster=False,
            col_cluster=False,
            yticklabels = False,
            cmap='RdBu_r',
            vmin=0, vmax=6,
            center = 2,
            #norm=divnorm,
            cbar_kws=cbar_kws)

    h.cax.set_position([0.2, .2, .03, .45])

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
    plt.savefig(outprefix + '_heatmap.png', bbox_inches = 'tight')
    plt.clf()
    

def SplitSamplesHom(cnvs, clusters, outprefix):

    #labels = clusters['ClusterNumber'].unique()

    boundaries = pd.DataFrame(cnvs[['CHR', 'START', 'END']])
    cnvs = cnvs.drop(['CHR', 'START', 'END'], axis=1).transpose()

    cells = cnvs.index.tolist()
    index_ = []
    for cell in cells:
        cellN = cell[4:]
        index_.append(int(cellN))
    cnvs.index = pd.Index(data=index_, dtype=int)
    cnvs['subsample'] = clusters['ClusterNumber']


    #median subsample size and std_dev
    sizes = cnvs.groupby('subsample').size()
    sizes.describe().to_csv(os.path.join(outprefix, "newick_clusters_sizes_stats.csv"))

    first_q = sizes.quantile(0.25)
    median = sizes.median()
    third_q = sizes.quantile(0.75)

    sns.distplot(sizes)
    plt.gca().axvline(first_q, color='k', linestyle='--')
    plt.gca().axvline(median, color='r', linestyle='--')
    plt.gca().axvline(third_q, color='b', linestyle='--')

    plt.legend({'25%':first_q, '50%':median, '75%':third_q})
    plt.savefig(os.path.join(outprefix, "newick_clusters_size_distribution.png"))
    plt.clf()

    #filter out subsamples with a number of cells smaller than the first quartile
    cnvs = cnvs.groupby('subsample').filter(lambda x : x['subsample'].count() > first_q)
    subsamples = dict(iter(cnvs.groupby('subsample')))

    #downsample very large clusters (> 3rd quartile)
    new_sizes = {}
    for key, sample in subsamples.items():
        if sample.shape[0] > third_q:
            sample = sample.sample(n=int(round(third_q)))
            subsamples[key] = sample
        new_sizes[key] = sample.shape[0]

    for k in subsamples.keys():
        subsample = subsamples[k].drop('subsample', axis=1).transpose()
        idx = 0
        for cname in boundaries.columns:
            subsample.insert(loc=idx, column=cname, value=boundaries[cname])
            idx = idx + 1
        subsamples[k] = subsample

    return subsamples, new_sizes

def SplitSamplesHet(cnvs, new_sizes):
    boundaries = pd.DataFrame(cnvs[['CHR', 'START', 'END']])
    cnvs = cnvs.drop(['CHR', 'START', 'END'], axis=1).transpose()

    cells = cnvs.index.tolist()
    index_ = []
    for cell in cells:
        cellN = cell[4:]
        index_.append(int(cellN))
    cnvs.index = pd.Index(data=index_, dtype=int)

    subsamples = {}
    for idx, size in zip(np.arange(len(new_sizes)),new_sizes):
        subsamples[idx] = cnvs.sample(n=size)
        #remove from cnvs rows which have been already assigned
        cnvs = cnvs.loc[~cnvs.index.isin(subsamples[idx].index)]
        
    
    for k in subsamples.keys():
        subsample = subsamples[k].transpose()
        idx = 0
        for cname in boundaries.columns:
            subsample.insert(loc=idx, column=cname, value=boundaries[cname])
            idx = idx + 1
        subsamples[k] = subsample
    return subsamples


def main():
    parser = argparse.ArgumentParser(description="Splits simulations in subsamples.")
    parser.add_argument("-f", "--file", required=True, type=str, help="SegCopy")
    parser.add_argument("-n", "--npy", required=True,type=str, help="Subtrees npy")
    parser.add_argument("-o", "--outprefix", required=True,type=str, help="Output path prefix")
    
    args = parser.parse_args()

    segcopy = pd.read_csv(args.file, sep="\t")
    #clusters = pd.read_csv(args.clusters, sep="\t", index_col=0)
    clones = np.load(args.npy)
    sizes = []

    for v in np.unique(clones):
        sizes.append(np.count_nonzero(clones == v))

    # 1. heterogeneous samples
    #outdir = os.path.join(args.outprefix, "het")
    if os.path.exists(args.outprefix) == False:
        os.mkdir(args.outprefix)
    subsamples_het = SplitSamplesHet(segcopy, sizes)
    for label in subsamples_het.keys():
        Heatmap(subsamples_het[label], os.path.join(args.outprefix, "SegCopy_"+str(label)))
        subsamples_het[label].to_csv(os.path.join(args.outprefix, "SegCopy_" + str(label)), sep="\t", index=False)

    # 2. homogeneous samples
    """
    outdir = os.path.join(args.outprefix, "hom")
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    subsamples_hom = SplitSamplesHom(segcopy, new_sizes.values(), outdir)
    for label in subsamples_hom.keys():
        Heatmap(subsamples_hom[label], os.path.join(outdir, "sample"+str(label)))
        subsamples_hom[label].to_csv(os.path.join(outdir, "sample"+str(label)+"_SegCopy"), sep="\t", index=False)
    """

if __name__ == "__main__":
    sys.exit(main())
