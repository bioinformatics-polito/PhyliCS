import time
import umap
import numpy as np
import operator
import seaborn as sns
import pandas as pd
from kneed import KneeLocator
from sklearn.neighbors import NearestNeighbors
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
from sklearn.base import ClusterMixin
from sklearn.mixture import GaussianMixture
from yellowbrick.cluster import KElbowVisualizer


#gmm = mixture.GaussianMixture(n_components=4, covariance_type="full").fit(pca_result)
#labels = gmm.predict(pca_result)

import argparse, os, sys

class GMClusters(GaussianMixture, ClusterMixin):
    def __init__(self, n_clusters=2, **kwargs):
        kwargs["n_components"] = n_clusters
        super(GMClusters, self).__init__(**kwargs)
    def fit(self, X):
        super(GMClusters, self).fit(X)
        self.labels_ = self.predict(X)
        return self 

class GMM_res:
    def __init__(self, k, bic):
        self.k = k
        self.bic = bic

def GMM(X, n_components, **kwargs):
    gm = GaussianMixture(n_components=n_components, **kwargs).fit(X)
    return GMM_res(n_components, -gm.bic(X))

def find_opt_K(matrix, outprefix, n_jobs):
    gm_bic= {}
    """
    for i in range(2,21):
        
        gm = GaussianMixture(n_components=i, warm_start=True, random_state=0).fit(matrix.values)
        #print("BIC for number of cluster(s) {}: {}".format(i,gm.bic(X_scaled)))
        #print("Log-likelihood score for number of cluster(s) {}: {}".format(i,gm.score(X_scaled)))
        #print("-"*100)
        gm_bic[i] = -gm.bic(matrix.values)
    """
    k_bic_list = Parallel(n_jobs=n_jobs, prefer="threads")(delayed(GMM)(X=matrix.values, n_components=i, warm_start=True, random_state=0) for i in range(2, 20))
    print(k_bic_list)
    for k_bic in k_bic_list:
        gm_bic[k_bic.k] = k_bic.bic
    optK = max(gm_bic.items(), key=operator.itemgetter(1))[0]
    plt.title("The Gaussian Mixture model BIC \nfor determining number of clusters\n",fontsize=16)
    plt.scatter(x=list(gm_bic.keys()),y=np.log(list(gm_bic.values())),s=150,edgecolor='k')
    plt.axvline(x=optK, color='k', linestyle='--', label="optK")
    plt.grid(True)
    plt.legend(fontsize=18)
    plt.xlabel("Number of clusters",fontsize=14)
    plt.ylabel("Log of Gaussian mixture BIC score",fontsize=15)
    plt.xticks([i for i in range(2,21)],fontsize=14)
    plt.yticks(fontsize=15)
    plt.savefig(outprefix+"_gaussian_bic_score.png")
    plt.clf()
    return optK


def cluster_and_output(k, matrix, outprefix, preproc):
    start_time = time.time()
    model = GaussianMixture(n_components=k, random_state=0).fit(matrix.values)
    end_time = time.time()
    labels_ = model.predict(matrix.values)
    clusters = pd.DataFrame(labels_, index=matrix.index, columns=["cluster"])
    clusters.to_csv(outprefix+"_gaussian_clusters.csv", sep="\t")
    with open(outprefix+"_gaussian_performance.csv", 'w+') as f:
        f.write("computation_time\t"+str(end_time - start_time))
    if preproc:
        standard_embedding = umap.UMAP(random_state=42).fit_transform(matrix.values)
        plt.scatter(standard_embedding[:,0],standard_embedding[:,1], c=labels_, cmap='rainbow')
    else:
        plt.scatter(matrix.values[:,0],matrix.values[:,1], c=labels_, cmap='rainbow')
    plt.suptitle("Gaussian Mixture Modelling result")
    plt.savefig(outprefix+"_gaussian_scatter.png")
    plt.clf()
    

def main():
    parser = argparse.ArgumentParser(description='Gaussian Mixture Modelling clustering')
    parser.add_argument('-f', '--file', type=str, required=True,
                    help='Principal components matrix.')
    parser.add_argument('-p', '--preproc', action='store_true', 
                    help="To be specified if the input file requires preprocessing")
    parser.add_argument('-j', '--jobs', type=int, 
                    help="Number of jobs to parallelize optimal K research.")
    parser.add_argument('-o', '--outprefix', type=str, required=True,
                    help='Output prefix')


    args = parser.parse_args()

    n_jobs = args.jobs

    if args.preproc:
        matrix_input = pd.read_csv(args.file, sep="\t", usecols = lambda column : column not in ['CHR', 'START', 'END']).transpose()
    else:
        matrix_input = pd.read_csv(args.file, index_col=0, sep='\t')
    outprefix = args.outprefix

    k = find_opt_K(matrix_input, outprefix, n_jobs)
    cluster_and_output(k, matrix_input, outprefix, args.preproc)

if __name__ == "__main__":
    sys.exit(main())