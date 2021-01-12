import sys
import operator
import argparse
import numpy as np
from math import sqrt
import pandas as pd
import seaborn as sns
import matplotlib 
from sklearn import metrics
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

parser = argparse.ArgumentParser(description="hdbscan")


parser.add_argument("input", metavar='SegCopy', action='store',
        help='cnvs file.', type=str)
parser.add_argument("outdir", metavar='outdir', action='store',
        help='Output directory path.', type=str)

parser.add_argument("n_pcs", metavar='n_pcs', action='store',
        help='Number of principal components.', type=int)


args=parser.parse_args()

input_file = args.input
outdir = args.outdir
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
tsne_results = tsne.fit_transform(cnvs.values)

#dbscan 
db = DBSCAN(eps=0.1, min_samples=10).fit(pca_results)
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)

print('Estimated number of clusters: %d' % n_clusters_)
print('Estimated number of noise points: %d' % n_noise_)
#print("Silhouette Coefficient: %0.3f"
      #% metrics.silhouette_score(pca_results, labels))

color_palette = sns.color_palette("hls", len(np.unique(labels)))
cluster_colors = [color_palette[x] if x >= 0
                  else (0.5, 0.5, 0.5)
                  for x in labels]

plt.scatter(*tsne_results.T, s=90, linewidth=0, c=cluster_colors, alpha=0.5)
plt.gcf().suptitle('dbscan clusters', fontsize=16, fontweight='bold' )
plt.gcf().set_size_inches(12,12)
plt.savefig(outdir+'/dbscan_clusters.png')
plt.clf()

