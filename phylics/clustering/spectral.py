import time
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import SpectralClustering
from yellowbrick.cluster import KElbowVisualizer

import argparse, os, sys

def locate_elbow(matrix, outprefix):
    model = SpectralClustering(assign_labels="discretize")
    visualizer = KElbowVisualizer(model, k=(2,20))
    visualizer.fit(matrix.values)
    visualizer.finalize()
    plt.savefig(outprefix+"_spectral_elbow.png")
    plt.clf()
    return visualizer.elbow_value_


def cluster_and_output(k, matrix, outprefix):
    start_time = time.time()
    model = SpectralClustering(n_clusters=k,assign_labels="discretize").fit(matrix.values)
    end_time = time.time()
    clusters = pd.DataFrame(model.labels_, index=matrix.index, columns=["cluster"])  
    clusters.to_csv(outprefix+"_spectral_clusters.csv", sep="\t")
    with open(outprefix+"_spectral_performance.csv", 'w+') as f:
        f.write("computation_time\t"+str(end_time - start_time))
    plt.scatter(matrix.values[:,0],matrix.values[:,1], c=model.labels_, cmap='rainbow')
    plt.suptitle("Spectral clustering result")
    plt.savefig(outprefix+"_spectral_scatter.png")
    plt.clf()
    

def main():
    parser = argparse.ArgumentParser(description='Spectral clustering')
    parser.add_argument('-f', '--file', type=str, required=True,
                    help='Principal components matrix.')
    parser.add_argument('-o', '--outprefix', type=str, required=True,
                    help='Output prefix')

    args = parser.parse_args()
    matrix_input = pd.read_csv(args.file, index_col=0, sep='\t')
    outprefix = args.outprefix

    k = locate_elbow(matrix_input, outprefix)
    cluster_and_output(k, matrix_input, outprefix)

if __name__ == "__main__":
    sys.exit(main())


