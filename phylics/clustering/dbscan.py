import time
import umap
import pandas as pd
import seaborn as sns
import numpy as np
from kneed import KneeLocator
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN


import argparse, os, sys

def find_opt_eps(matrix, outprefix):
    """
    In layman’s terms, we find a suitable value for epsilon 
    by calculating the distance to the nearest n points for each point, 
    sorting and plotting the results. Then we look to see where the 
    change is most pronounced (think of the angle between your arm 
    and forearm) and select that as epsilon.

    https://iopscience.iop.org/article/10.1088/1755-1315/31/1/012012/pdf
    """
    
    X = matrix.values
    neigh = NearestNeighbors(n_neighbors=2)
    nbrs = neigh.fit(X)
    distances, indices = nbrs.kneighbors(X)
    distances = np.sort(distances, axis=0)
    y = distances[:,1]
    x = range(0, len(y))
    kn = KneeLocator(x, y, curve='convex', direction='increasing')
    eps = y[kn.knee]
    """
    plt.xlabel('samples')
    plt.ylabel('eps')
    plt.plot(x, y, 'bx-')
    plt.axhline(y=eps, color='k', linestyle='--', label="optEps")
    plt.gca().grid(True)
    plt.legend(fontsize=18)
    plt.savefig(outprefix+"_dbscan_eps.png")
    plt.clf()
    """
    return eps

def cluster_and_output(matrix, outprefix, eps, preproc):
    start_time = time.time()
    model = DBSCAN(eps=eps).fit(matrix.values)
    end_time = time.time()
    clusters = pd.DataFrame(model.labels_, index=matrix.index, columns=["cluster"])
    clusters.to_csv(outprefix+"_dbscan_clusters.csv", sep="\t")
    with open(outprefix+"_dbscan_performance.csv", 'w+') as f:
        f.write("computation_time\t"+str(end_time - start_time))
        f.write("DBSCAN_eps\t"+str(eps))
    color_palette = sns.color_palette("rainbow", len(np.unique(model.labels_)))
    cluster_colors = [color_palette[x] if x >= 0
                      else (0.5, 0.5, 0.5)
                      for x in model.labels_]
    if preproc:
        standard_embedding = umap.UMAP(random_state=42).fit_transform(matrix.values)
        plt.scatter(standard_embedding[:,0],standard_embedding[:,1], c=model.labels_, cmap='rainbow')
    else:
        plt.scatter(matrix.values[:,0],matrix.values[:,1], c=model.labels_, cmap='rainbow')
    plt.suptitle("DBSCAN clustering result")
    plt.savefig(outprefix+"_dbscan_scatter.png")
    plt.clf()
    

def main():
    parser = argparse.ArgumentParser(description='DBSCAN clustering')
    parser.add_argument('-f', '--file', type=str, required=True,
                    help='Principal components matrix.')
    parser.add_argument('-p', '--preproc', action='store_true', 
                    help="To be specified if the input file requires preprocessing")
    parser.add_argument('-o', '--outprefix', type=str, required=True,
                    help='Output prefix')


    args = parser.parse_args()

    if args.preproc:
        matrix_input = pd.read_csv(args.file, sep="\t", usecols = lambda column : column not in ['CHR', 'START', 'END']).transpose()
    else:
        matrix_input = pd.read_csv(args.file, index_col=0, sep='\t')
    outprefix = args.outprefix

    eps = find_opt_eps(matrix_input, outprefix)
    with open(outprefix+"_dbscan_performance.csv", 'a') as f:
        f.write("DBSCAN_eps\t"+str(eps))
    #cluster_and_output(matrix_input, outprefix, eps, args.preproc)

if __name__ == "__main__":
    sys.exit(main())


