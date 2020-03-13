import sys
import operator
import argparse
import numpy as np
from math import sqrt
import pandas as pd
import seaborn as sns
import matplotlib 
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import hdbscan
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

parser = argparse.ArgumentParser(description="hdbscan")


parser.add_argument("input", metavar='SegCopy', action='store',
        help='cnvs file.', type=str)
parser.add_argument("outdir", metavar='outdir', action='store',
        help='Output directory path.', type=str)

parser.add_argument("n_pcs", metavar='n_pcs', action='store',
        help='Number of principal components.', type=int)

parser.add_argument("--algorithm", metavar='algorithm', action='store', choices=['best', 'generic', 'prims_kdtree', 'prims_balltree', 'boruvka_kdtree', 'boruvka_balltree'], default='best',
        help='Number of principal components.', type=str)


args=parser.parse_args()

input_file = args.input
outdir = args.outdir
npcs = args.n_pcs
algorithm = args.algorithm


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


#hdbscan
clusterer_dict = {}
for min_size in range(3,15+1):
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_size, min_samples=1,  cluster_selection_method='eom', algorithm=algorithm).fit(pca_results)
  
    count = 0
    for p in clusterer.probabilities_:
        if p < 0.05:
            count = count +1
    clusterer_dict[min_size] = count

optSize=min(clusterer_dict.items(), key=operator.itemgetter(1))[0]
print("optimal min cluster size: " + str(optSize))


clusterer = hdbscan.HDBSCAN(min_cluster_size=optSize, min_samples=1, cluster_selection_epsilon=0.5, cluster_selection_method='eom', algorithm=algorithm, prediction_data=True).fit(pca_results)

#print("n probabilties: " + str(len(clusterer.probabilities_)))
#print(clusterer.probabilities_)

#print(clusterer.labels_)
color_palette = sns.color_palette("hls", len(np.unique(clusterer.labels_)))
cluster_colors = [color_palette[x] if x >= 0
                  else (0.5, 0.5, 0.5)
                  for x in clusterer.labels_]

plt.scatter(*tsne_results.T, s=90, linewidth=0, c=cluster_colors, alpha=0.5)
plt.gcf().suptitle('hdbscan clusters', fontsize=16, fontweight='bold' )
plt.gcf().set_size_inches(12,12)
plt.savefig(outdir+'/hdbscan_clusters.png')
plt.clf()

clusterer.condensed_tree_.plot(select_clusters=True,selection_palette=color_palette)
plt.gcf().suptitle('hdbscan condensed tree', fontsize=16, fontweight='bold' )
plt.gcf().set_size_inches(37,21)
plt.savefig(outdir+'/hdbscan_condensed_tree.png')
plt.clf()


#hdbscab+cnv profiles
cnvs['cluster'] = clusterer.labels_
cnvs = cnvs.sort_values(by='cluster')
labels = cnvs['cluster'].values
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
plt.savefig(outdir+'/clusters_heatmap.png')
plt.clf()

#cnvs[cnvs["cluster"] == 0].to_csv(outdir+"/to_be_filtered.tsv", sep="\t")

"""
#soft clustering
soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
color_palette = sns.color_palette('hls', 20)
cluster_colors = [color_palette[np.argmax(x)]
                  for x in soft_clusters]
plt.scatter(*tsne_results.T, s=90, linewidth=0, c=cluster_colors, alpha=0.50)
plt.gcf().suptitle('hdbscan soft clusters', fontsize=16, fontweight='bold' )
plt.gcf().set_size_inches(12,12)
plt.savefig(outdir+'/hdbscan_soft_clusters.png')
plt.clf()

soft_cluster_labels = []
for c in soft_clusters:
    soft_cluster_labels.append(np.argmax(c))

#soft-hdbscan+cnvs
#hdbscab+cnv profiles
cnvs['cluster'] = soft_cluster_labels
cnvs = cnvs.sort_values(by='cluster')
labels = cnvs['cluster'].values

cluster_colors = [color_palette[x] if x >= 0
                  else (0.5, 0.5, 0.5)
                  for x in labels]

cnvs = cnvs.drop(['cluster'], axis=1)
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
plt.savefig(outdir+'/soft_clusters_heatmap.png')
plt.clf()
"""





