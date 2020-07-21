import sys
import operator
import argparse
import numpy as np
from math import sqrt
import pandas as pd
import seaborn as sns
from matplotlib import cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import umap
import hdbscan
from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA


parser = argparse.ArgumentParser(description="hdbscan")


parser.add_argument("input", metavar='SegCopy', action='store',
        help='cnvs file.', type=str)

parser.add_argument("input2", metavar='SegStats', action='store',
        help='cell info file.', type=str)

parser.add_argument("outdir", metavar='outdir', action='store',
        help='Output directory path.', type=str)

parser.add_argument("n_pcs", metavar='n_pcs', action='store',
        help='Number of principal components.', type=int)

parser.add_argument("--algorithm", metavar='algorithm', action='store', choices=['best', 'generic', 'prims_kdtree', 'prims_balltree', 'boruvka_kdtree', 'boruvka_balltree'], default='best',
        help='Number of principal components.', type=str)

parser.add_argument("--min_cluster_size", metavar='N', action='store',
        help='Min cluster size (hdbscan parameter).', type=str)


args=parser.parse_args()

input_file = args.input
input2 = args.input2
outdir = args.outdir
npcs = args.n_pcs
algorithm = args.algorithm
min_cluster_size = None

if args.min_cluster_size:
    min_cluster_size = min_cluster_size

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

# cell id dictionary
cell_id_dict = dict(zip(list(range(len(cnvs))), cnvs.index ))
cell_ids = cnvs.index 

# retrieve the list of 'semi-noisy' cells (mad > 3 quartile)
df_stats = pd.read_csv(input2, sep="\t")
sns.violinplot(x=df_stats['Disp'], scale='count')
plt.savefig(outdir + "/cells_disp.png")
plt.clf()


trd_quartile = df_stats['Disp'].quantile([0.5]).values[0]
noisy_cells = df_stats[df_stats['Disp'] >= trd_quartile].index.values

#df_metrics = pd.read_csv(input2, usecols = ['barcode', 'semi_noisy'])
#noisy_cells = df_metrics['barcode'][df_metrics['semi_noisy'] == 1].str[:-2].values

#preliminary dimensionality reduction

pca = PCA(n_components=npcs)
pca_result = pca.fit_transform(cnvs.values)

#umap
N=len(cnvs)
standard_embedding = umap.UMAP( 
                                #n_neighbors=15,
                                #min_dist=0.0,
                                n_components=2,
                                random_state=42
        
                            ).fit_transform(pca_result)
#plt.scatter(standard_embedding[:, 0], standard_embedding[:, 1], s=0.1, cmap='Spectral');

#standard_embedding = umap.UMAP(random_state=42).fit_transform(pca_result)

#hdbscan
if min_cluster_size == None:
    clusterer_dict = {}
    for min_size in range(3,15+1):
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_size, min_samples=1,  cluster_selection_method='eom', algorithm=algorithm).fit(pca_result) 
        count = 0
        for p in clusterer.probabilities_:
            if p < 0.05:
                count = count +1
        clusterer_dict[min_size] = count

    min_cluster_size=min(clusterer_dict.items(), key=operator.itemgetter(1))[0]
    print("optimal min cluster size: " + str(min_cluster_size))


clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=1, cluster_selection_epsilon=0.5, cluster_selection_method='eom', algorithm=algorithm, prediction_data=True).fit(pca_result)


clustered = (clusterer.labels_ >= 0)
unclustered = (clusterer.labels_ < 0)

clustered_cells = cell_ids[clustered]
unclustered_cells = cell_ids[unclustered]


# print distribution of MAD of unclustered data

"""
clustered_df = df_stats.loc[clustered_cells, : ]
sns.violinplot(x=clustered_df['Disp'], scale='count')
plt.savefig(outdir + "/clustered_cells_disp.png")
plt.clf()
"""

clusters_df = pd.DataFrame(columns=['Disp', 'Clustering'])

clustered_df = pd.DataFrame( data=df_stats.loc[clustered_cells, ['Disp']].values, index=clustered_cells, columns=['Disp'])
clustered_df['Clustering'] = 'Clustered'

unclustered_df = pd.DataFrame( data=df_stats.loc[unclustered_cells, ['Disp']].values, index=unclustered_cells, columns=['Disp'])
unclustered_df['Clustering'] = 'Unclustered'


clusters_df = clusters_df.append(clustered_df)
clusters_df = clusters_df.append(unclustered_df)

if len(unclustered_cells) > 0:
    #unclustered_df = df_stats.loc[unclustered_cells, : ]
    print("Number of unclustered cells: " + str(unclustered_df.shape[0]))

    intersection = []
    for x in unclustered_cells:
        if x in noisy_cells:
            intersection.append(x)

    perc_unclustered_noisy = len(intersection)/len(unclustered_cells)
    print("Percentage of unclustered data with high bin-to-bin dispertion: " + str(perc_unclustered_noisy))

    #mannwhitneyu test

    w, p = mannwhitneyu(clustered_df['Disp'], unclustered_df['Disp'])
    print("mannwhitneyu test p = " + str(p))


    #clustered_df = df_stats.loc[clustered_cells, : ]

    sns.set(style="whitegrid", font_scale=1.5)
    sns.violinplot(x='Clustering', y='Disp', data=clusters_df, scale='count')
    #sns.violinplot(x=unclustered_df['Disp'], scale='count',  c='green', label='unclustered', ax = axes[1])
    

    plt.gca().set_title('Statistical significance (mannwhitneyu test) = {}'.format(p))
    plt.gcf().set_size_inches(12,12)
    plt.savefig(outdir + "/cells_disp_violin_plot.png")
    plt.clf()
    


# hdbscan results

color_palette = sns.color_palette("hls", len(np.unique(clusterer.labels_)))
cluster_colors = [color_palette[x] if x >= 0
                  else (0.5, 0.5, 0.5)
                  for x in clusterer.labels_]

plt.scatter(standard_embedding[: , 0],
            standard_embedding[: , 1],
            c=cluster_colors,
            #s=0.9,
            alpha=0.5)

"""
plt.scatter(standard_embedding[clustered, 0],
            standard_embedding[clustered, 1],
            c=clusterer.labels_[clustered],
            s=0.9,
            cmap='Spectral');
"""
#plt.scatter(*tsne_results.T, s=90, linewidth=0, c=cluster_colors, alpha=0.5)
plt.gcf().suptitle('hdbscan clusters pca', fontsize=16, fontweight='bold' )
plt.gcf().set_size_inches(12,12)
#plt.show()
plt.savefig(outdir+'/hdbscan_umap.png')
plt.clf()

clustered_perc = np.sum(clustered) / cnvs.values.shape[0]
print("Percentage of clustered data: " + str(clustered_perc))

"""
clusterer.condensed_tree_.plot(select_clusters=True,selection_palette=color_palette)
plt.gcf().suptitle('hdbscan condensed tree', fontsize=16, fontweight='bold' )
plt.gcf().set_size_inches(37,21)
plt.savefig(outdir+'/hdbscan_condensed_tree.png')
plt.clf()
"""

#hdbscab+cnv profiles
cnvs['cluster'] = clusterer.labels_
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

