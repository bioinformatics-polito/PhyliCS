#!/usr/bin/envs python

import sys, os
import argparse

import pandas as pd 
import numpy as np 

import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from ete3 import Tree, NodeStyle, TreeStyle

def main():
    parser = argparse.ArgumentParser(description="Plot clusters on newick tree.")

    parser.add_argument("-t", "--tree", required=True, 
                        help="Newick tree file", type=str)
    parser.add_argument("-c", "--clusters", required=True, 
                        help="Clusters file", type=str)
    parser.add_argument("-o", "--outpath", required=True, 
                        help="Output folder path", type=str)
    
    args = parser.parse_args()

    clusters = pd.read_csv(args.clusters, sep="\t")
    labels = clusters["ClusterNumber"]
    
    color_palette = sns.color_palette("hls", len(np.unique(labels))).as_hex()
    cluster_lut = dict(zip(np.unique(labels), color_palette))
    cluster_colors=pd.Series(labels).map(cluster_lut)
    cluster_colors_map = dict(zip(clusters["SequenceName"], cluster_colors))

    t = Tree(args.tree)

    # Basic tree style
    ts = TreeStyle()
    ts.mode = 'c'

    # Creates an independent node style for each leaf, which is
    # initialized with the corresponding cluster color.
    for n in t.traverse():
        if n.is_leaf():      
            nstyle = NodeStyle()  
            nstyle["fgcolor"] = cluster_colors_map[int(n.name)]
            nstyle["size"] = 15
            n.set_style(nstyle)

    t.render(os.path.join(args.outpath, "tree_clusters.png"),tree_style=ts)

if __name__ == "__main__":
    sys.exit(main())