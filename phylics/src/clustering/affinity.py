#!/usr/bin/envs python

import time
import umap
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import AffinityPropagation
from yellowbrick.cluster import KElbowVisualizer

import argparse, os, sys

"""
def locate_elbow(matrix, outprefix):
    model = AffinityPropagation()
    visualizer = KElbowVisualizer(model, k=(3,15), metric='calinski_harabasz')
    visualizer.fit(matrix.values)
    visualizer.finalize()
    plt.savefig(outprefix+"_affinity_elbow.png")
    plt.clf()
    return visualizer.elbow_value_
"""

def cluster_and_output(matrix, outprefix, preproc):
    start_time = time.time()
    model = AffinityPropagation().fit(matrix.values)
    end_time = time.time()
    clusters = pd.DataFrame(model.labels_, index=matrix.index, columns=["cluster"])
    clusters.to_csv(outprefix+"_affinity_clusters.csv", sep="\t")
    with open(outprefix+"_affinity_performance.csv", 'w+') as f:
        f.write("computation_time\t"+str(end_time - start_time))
    if preproc:
        standard_embedding = umap.UMAP(random_state=42).fit_transform(matrix.values)
        plt.scatter(standard_embedding[:,0],standard_embedding[:,1], c=model.labels_, cmap='rainbow')
    else:
        plt.scatter(matrix.values[:,0],matrix.values[:,1], c=model.labels_, cmap='rainbow')
    plt.suptitle("Affinity propagation clustering result")
    plt.savefig(outprefix+"_affinity_scatter.png")
    plt.clf()
    

def main():
    parser = argparse.ArgumentParser(description='Affinity propagation clustering')
    parser.add_argument('-f', '--file', type=str, required=True,
                    help='Principal components matrix.')
    parser.add_argument('-p', '--preproc', action='store_true', 
                    help="To be specified if the input file requires preprocessing")
    parser.add_argument('-o', '--outprefix', type=str, required=True,
                    help='Output prefix')

    args = parser.parse_args()

    outprefix = args.outprefix
    if args.preproc:
        matrix_input = pd.read_csv(args.file, sep="\t", usecols = lambda column : column not in ['CHR', 'START', 'END']).transpose()
    else:
        matrix_input = pd.read_csv(args.file, index_col=0, sep='\t')

    #k = locate_elbow(matrix_input, outprefix)
    cluster_and_output(matrix_input, outprefix, args.preproc)

if __name__ == "__main__":
    sys.exit(main())



