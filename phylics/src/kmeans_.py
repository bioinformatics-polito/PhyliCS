from sklearn.cluster import KMeans
import numpy as np
import sys
import umap
import operator
import argparse
from math import sqrt
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

parser = argparse.ArgumentParser(description="hdbscan")


parser.add_argument("input", metavar='SegCopy', action='store',
        help='cnvs file.', type=str)
parser.add_argument("outdir", metavar='outdir', action='store',
        help='Output directory path.', type=str)

parser.add_argument("n_pcs", metavar='n_pcs', action='store',
        help='Number of principal components.', type=int)

parser.add_argument("k", metavar='K', action='store',
        help='k-means centroids.', type=int)


args=parser.parse_args()

input_file = args.input
outdir = args.outdir
npcs = args.n_pcs
k = args.k


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
pca_result = pca.fit_transform(cnvs)

#umap
N=len(cnvs)
standard_embedding = umap.UMAP(
                                #n_neighbors=15,
                                #min_dist=0.0,
                                n_components=2,
                                random_state=42

                            ).fit_transform(pca_result)


#kmeans
kmeans = KMeans(n_clusters=k, random_state=0).fit(pca_result)

color_palette = sns.color_palette("hls", len(np.unique(kmeans.labels_)))
cluster_colors = [color_palette[x] if x >= 0
                  else (0.5, 0.5, 0.5)
                  for x in kmeans.labels_]

plt.scatter(standard_embedding[: , 0],
            standard_embedding[: , 1],
            c=cluster_colors,
            #s=0.9,
            alpha=0.5)

plt.gcf().suptitle('hdbscan clusters pca', fontsize=16, fontweight='bold' )
plt.gcf().set_size_inches(12,12)
#plt.show()
plt.savefig(outdir+'/hdbscan_umap.png')
plt.clf()


#kmeans+cnv profiles
cnvs['cluster'] = kmeans.labels_
cnvs = cnvs.sort_values(by='cluster')
labels = cnvs['cluster'].values
#color_palette = sns.color_palette("Spectral", len(np.unique(cnvs['cluster'][cnvs['cluster'] >= 0].values)))
cluster_colors = [color_palette[x] if x >= 0
                  else (0.5, 0.5, 0.5)
                  for x in labels]
#print(row_colors)
cbar_kws={"ticks":np.arange(0,7,1)}
h = sns.clustermap(cnvs,
        row_cluster=False,
        col_cluster=False,
        yticklabels = False,
        row_colors=cluster_colors,
        cmap='RdBu_r',
        vmin=0, vmax=6,
        center = 2,
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
#plt.show()
plt.savefig(outdir+'/clusters_heatmap.png')
plt.clf()


