#!/usr/bin/env python
import sys
import os
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import argparse
import seaborn as sns


def DownSample():
    return NotImplemented

def SplitSamples(cnvs, clusters, outprefix):

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

    #filter out subsamples with a number of cells smaller than the first quartile
    cnvs = cnvs.groupby('subsample').filter(lambda x : x['subsample'].count() > first_q)
    subsamples = dict(iter(cnvs.groupby('subsample')))

    #downsample very large clusters (> 3rd quartile)
    
    
    for k in subsamples.keys():
        subsample = subsamples[k].drop('subsample', axis=1).transpose()
        idx = 0
        for cname in boundaries.columns:
            subsample.insert(loc=idx, column=cname, value=boundaries[cname])
            idx = idx + 1
        subsamples[k] = subsample

    return subsamples


def main():
    parser = argparse.ArgumentParser(description="Splits simulations in subsamples.")
    parser.add_argument("-f", "--file", required=True, type=str, help="SegCopy")
    parser.add_argument("-c", "--clusters", required=True,type=str, help="Newick clusters")
    parser.add_argument("-o", "--outprefix", required=True,type=str, help="Output path prefix")
    
    args = parser.parse_args()

    segcopy = pd.read_csv(args.file, sep="\t")
    clusters = pd.read_csv(args.clusters, sep="\t", index_col=0)

    subsamples = SplitSamples(segcopy, clusters, args.outprefix)

    '''
    for label in subsamples.keys():
        subsamples[label].to_csv(args.outprefix+"sample"+str(label)+"_SegCopy", sep="\t", index=False)
    '''


if __name__ == "__main__":
    sys.exit(main())
