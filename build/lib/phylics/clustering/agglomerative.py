#!/usr/bin/envs python
import operator
import time
import pandas as pd
import umap
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from yellowbrick.cluster import KElbowVisualizer
from sklearn.metrics import silhouette_score

import argparse, os, sys

def silhouette_score_(matrix, outprefix, linkage):
    s_dict = {}
   
    if linkage == "ward":
        affinity="euclidean"
    else:
        affinity="l1"   
    for k in range(2,21):
        model = AgglomerativeClustering(n_clusters=k, affinity=affinity, linkage=linkage).fit(matrix.values)
        s_dict[k] = silhouette_score(matrix.values, model.labels_)
        
    optK = max(s_dict.items(), key=operator.itemgetter(1))[0]
    lists = sorted(s_dict.items())
    x, y = zip(*lists) 
    plt.title("Agglomerative model (linkage= " + linkage + ") silhouette score \nfor determining number of clusters\n",fontsize=16)
    plt.scatter(x=x,y=y,s=150,edgecolor='k')
    plt.axvline(x=optK, color='k', linestyle='--', label="optK")
    plt.grid(True)
    plt.legend(fontsize=18)
    plt.xlabel("Number of clusters",fontsize=14)
    plt.ylabel("Silhouette score",fontsize=15)
    plt.xticks([i for i in range(2,21)],fontsize=14)
    plt.yticks(fontsize=15)
    plt.savefig(outprefix+"_agglomerative" + linkage +"_silhouette_score.png")
    plt.clf()
    return optK

def locate_elbow(matrix, outprefix, linkage):
    if linkage == "ward":
        affinity="euclidean"
    else:
        affinity="l1"
    model = AgglomerativeClustering(affinity=affinity, linkage=linkage)
    visualizer = KElbowVisualizer(model, k=(2,50), metric='silhouette')
    visualizer.fit(matrix.values)
    visualizer.finalize()
    plt.savefig(outprefix+"_agglomerative_" + linkage + "_elbow.png")
    plt.clf()
    return visualizer.elbow_value_


def cluster_and_output(k, matrix, outprefix, linkage, preproc):
    start_time = time.time()
    if linkage == "ward":
        affinity="euclidean"
    else:
        affinity="l1"
    model = AgglomerativeClustering(n_clusters=k, affinity=affinity, linkage=linkage).fit(matrix.values)
    end_time = time.time()
    clusters = pd.DataFrame(model.labels_, index=matrix.index, columns=["cluster"]) 
    clusters.to_csv(outprefix+"_agglomerative_" + linkage + "_clusters.csv", sep="\t")
    with open(outprefix+"_agglomerative_" + linkage + "_performance.csv", 'w+') as f:
        f.write("computation_time\t"+str(end_time - start_time))
    if preproc:
        standard_embedding = umap.UMAP(random_state=42).fit_transform(matrix.values)
        plt.scatter(standard_embedding[:,0],standard_embedding[:,1], c=model.labels_, cmap='rainbow')
    else:
        plt.scatter(matrix.values[:,0],matrix.values[:,1], c=model.labels_, cmap='rainbow')
    plt.suptitle("Agglomerative clustering result (linkage = {})".format(linkage))
    plt.savefig(outprefix+"_agglomerative_" + linkage + "_scatter.png")
    plt.clf()
    

def main():
    parser = argparse.ArgumentParser(description='Agglomerative clustering')
    parser.add_argument('-f', '--file', type=str, required=True,
                    help='Principal components matrix.')
    parser.add_argument('-p', '--preproc', action='store_true', 
                    help="To be specified if the input file requires preprocessing")
    parser.add_argument('-o', '--outprefix', type=str, required=True,
                    help='Output prefix')
    parser.add_argument('-l', '--linkage', type=str, required=True,
                    help='Linkage method')

    args = parser.parse_args()

    if args.preproc:
        matrix_input = pd.read_csv(args.file, sep="\t", usecols = lambda column : column not in ['CHR', 'START', 'END']).transpose()
    else:
        matrix_input = pd.read_csv(args.file, index_col=0, sep='\t')    

    outprefix = args.outprefix
    linkage = args.linkage

    k = silhouette_score_(matrix_input, outprefix, linkage)
    cluster_and_output(k, matrix_input, outprefix, linkage, args.preproc)

if __name__ == "__main__":
    sys.exit(main())