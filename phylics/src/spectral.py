import sys
import operator
import argparse
import scipy
from scipy.sparse import csgraph
# from scipy.sparse.linalg import eigsh
from numpy import linalg as LA
import numpy as np
from math import sqrt
import pandas as pd
import seaborn as sns
import matplotlib 
from sklearn import metrics
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import SpectralClustering 
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

def getAffinityMatrix(coordinates, k = 7):
    """
    Calculate affinity matrix based on input coordinates matrix and the numeber
    of nearest neighbours.
    
    Apply local scaling based on the k nearest neighbour
        References:
    https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering.pdf
    """
    # calculate euclidian distance matrix
    dists = squareform(pdist(coordinates)) 
    
    # for each row, sort the distances ascendingly and take the index of the 
    #k-th position (nearest neighbour)
    knn_distances = np.sort(dists, axis=0)[k]
    knn_distances = knn_distances[np.newaxis].T
    
    # calculate sigma_i * sigma_j
    local_scale = knn_distances.dot(knn_distances.T)

    affinity_matrix = dists * dists
    affinity_matrix = -affinity_matrix / local_scale
    # divide square distance matrix by local scale
    affinity_matrix[np.where(np.isnan(affinity_matrix))] = 0.0
    # apply exponential
    affinity_matrix = np.exp(affinity_matrix)
    np.fill_diagonal(affinity_matrix, 0)
    return affinity_matrix

def eigenDecomposition(A, plot = True, topK = 5):
    """
    :param A: Affinity matrix
    :param plot: plots the sorted eigen values for visual inspection
    :return A tuple containing:
    - the optimal number of clusters by eigengap heuristic
    - all eigen values
    - all eigen vectors

    This method performs the eigen decomposition on a given affinity matrix,
    following the steps recommended in the paper:
    1. Construct the normalized affinity matrix: L = D−1/2ADˆ −1/2.
    2. Find the eigenvalues and their associated eigen vectors
    3. Identify the maximum gap which corresponds to the number of clusters
    by eigengap heuristic

    References:
    https://papers.nips.cc/paper/2619-self-tuning-spectral-clustering.pdf
    http://www.kyb.mpg.de/fileadmin/user_upload/files/publications/attachments/Luxburg07_tutorial_4488%5b0%5d.pdf
    """
    L = csgraph.laplacian(A, normed=True)
    n_components = A.shape[0]

    # LM parameter : Eigenvalues with largest magnitude (eigs, eigsh), that is, largest eigenvalues in
    # the euclidean norm of complex numbers.
#     eigenvalues, eigenvectors = eigsh(L, k=n_components, which="LM", sigma=1.0, maxiter=5000)
    eigenvalues, eigenvectors = LA.eig(L)

    if plot:
        plt.title('Largest eigen values of input matrix')
        plt.scatter(np.arange(len(eigenvalues)), eigenvalues)
        plt.grid()

    # Identify the optimal number of clusters as the index corresponding
    # to the larger gap between eigen values
    index_largest_gap = np.argsort(np.diff(eigenvalues))[::-1][:topK]
    nb_clusters = index_largest_gap + 1

    return nb_clusters, eigenvalues, eigenvectors


parser = argparse.ArgumentParser(description="hdbscan")


parser.add_argument("input", metavar='SegCopy', action='store',
        help='cnvs file.', type=str)
#parser.add_argument("outdir", metavar='outdir', action='store',
#        help='Output directory path.', type=str)

parser.add_argument("n_pcs", metavar='n_pcs', action='store',
        help='Number of principal components.', type=int)


args=parser.parse_args()

input_file = args.input
#outdir = args.outdir
npcs = args.n_pcs


df = pd.read_csv(input_file, sep="\t")

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

#preliminary dimensionality reduction
pca = PCA(n_components=npcs)
pca_results = pca.fit_transform(cnvs)

#tsne
N=len(cnvs)
optPerp=round(sqrt(N))
tsne = TSNE(n_components=2, perplexity=optPerp, n_iter=10000)
tsne_results = tsne.fit_transform(pca_results)

plt.scatter(*tsne_results.T, s=50, linewidth=0, alpha=0.5)
plt.show()

#spectral clustering
affinity_matrix = getAffinityMatrix(cnvs, k = 5)
k, _,  _ = eigenDecomposition(affinity_matrix)
print(f'Optimal number of clusters {k}')
plt.show()

nb_clusters = k[3]

clustering = SpectralClustering(n_clusters=nb_clusters, assign_labels="discretize", random_state=0).fit(pca_results)
y_pred = clustering.labels_

#plt.figure(figsize=(14,6))
#plt.subplot(121)
plt.title(f'Spectral clustering results ')
plt.scatter(*tsne_results.T, s=50, c = y_pred,  linewidth=0, alpha=0.5)
plt.show()


cnvs['cluster'] = y_pred
cnvs = cnvs.sort_values(by='cluster')
labels = cnvs['cluster'].values
color_palette = sns.color_palette("hls", len(np.unique(y_pred)))
cluster_colors = [color_palette[x] if x >= 0
                  else (0.5, 0.5, 0.5)
                  for x in labels]

cbar_kws={"ticks":np.arange(0,7,1)}
h = sns.clustermap(cnvs,
	row_cluster=False,
        col_cluster=False,
        yticklabels = False,
        row_colors=cluster_colors,
        cmap='RdBu_r',
        vmin=0, vmax=6,
        center=2,
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
plt.show()

