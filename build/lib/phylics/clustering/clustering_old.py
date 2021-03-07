import time
import pandas as pd
import operator
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.feature_selection import VarianceThreshold
from sklearn import metrics
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from sklearn.cluster import Birch
from yellowbrick.cluster import KElbowVisualizer

import argparse, os, sys

def preprocessing(matrix, outprefix):
    # Low variance feature filtering
    selector = VarianceThreshold() #filters out features with 0 variance from cell to cell (not informative for clustering)
    df_filtered = pd.DataFrame(data=selector.fit_transform(matrix), index=matrix.index)

    # PCA (keep only significant features)
    #n_samples = len(df)
    #n_features = len(df.columns)
    PCs = 100
    pca = PCA(n_components=PCs)
    pca_result = pd.DataFrame(data=pca.fit_transform(df_filtered.values), index=matrix.index)

    expl_var_ratio = pca.explained_variance_ratio_
    expl_var_ratio_df = pd.DataFrame(data=pca.explained_variance_ratio_.reshape(-1, len(pca.explained_variance_ratio_)))

    expl_var_ratio_perm = pd.DataFrame() #columns = components, rows = permutations
    n_perm = 50
    for k in range(n_perm):
        df_perm = df_filtered.apply(np.random.permutation, 0)
        pca_perm = pca.fit_transform(df_perm.values)
        expl_var_ratio_single_df = pd.DataFrame(data=pca.explained_variance_ratio_.reshape(-1, len(pca.explained_variance_ratio_)))
        expl_var_ratio_perm = expl_var_ratio_perm.append(expl_var_ratio_single_df)
    mean_expl_var_ratio_perm = expl_var_ratio_perm.mean().values


    plt.figure()
    plt.plot(expl_var_ratio, color='green', label='Explained by PCs')
    plt.plot(mean_expl_var_ratio_perm, color='red', label='Explained by chance')

    plt.gca().set_xlabel('number of components',  fontsize=18)
    plt.gca().set_ylabel('cumulative explained variance',  fontsize=18)

    plt.xticks(np.arange(0, PCs, 1))

    plt.gca().legend(fontsize=18)
    plt.gca().grid(True)


    plt.suptitle("Explained variance ratio",  fontsize=22)
    plt.gca().tick_params(axis="x", labelsize=12)
    plt.gca().tick_params(axis="y", labelsize=12)
    plt.gcf().set_size_inches(30, 12)
    plt.savefig(outprefix+"_expl_var_ratio.png")
    plt.gcf().clf()


    # compute a pvalue  of how the observed variance is different from the permuted variance for each PC

    for c in expl_var_ratio_df.columns:
        expl_var_ratio_perm.loc[expl_var_ratio_perm[c] > expl_var_ratio_df.loc[0, c], c] = 1
        expl_var_ratio_perm.loc[expl_var_ratio_perm[c] <= expl_var_ratio_df.loc[0, c], c] = 0

    pvals = expl_var_ratio_perm.mean()

    optPCs = pvals[pvals.gt(0.05)].index[0]

    plt.figure()
    plt.plot(pvals, color='red')

    plt.axvline(x=optPCs, color='blue', linestyle='--', label="optPCs")

    plt.gca().set_xlabel('principal components', fontsize=18)
    plt.gca().set_ylabel('pvalue', fontsize=18)

    plt.xticks(np.arange(0, PCs, 1))
    plt.gca().grid(True)
    plt.legend(fontsize=18)
    plt.suptitle("Significance of principal components", fontsize=22)
    plt.gca().tick_params(axis="x", labelsize=12)
    plt.gca().tick_params(axis="y", labelsize=12)
    plt.gcf().set_size_inches(30, 12)
    plt.savefig(outprefix+"_pcs_pvals.png")
    plt.gcf().clf()

    #cumulative variance
    plt.figure()

    plt.plot(np.cumsum(expl_var_ratio))
    plt.gcf().set_size_inches(30, 12)
    plt.xticks(np.arange(0, PCs, 1))
    plt.gca().grid(True)
    plt.gca().tick_params(axis="x", labelsize=12)
    plt.gca().tick_params(axis="y", labelsize=12)
    plt.xlabel('Number of Components', fontsize=18)
    plt.ylabel('Variance (%)', fontsize=18) #for each component
    plt.title('Cumulative explained Variance ratio', fontsize=22)
    plt.savefig(outprefix+"_cum_expl_var_ratio.png")
    plt.gcf().clf()

    pca_result[pca_result.columns[0:optPCs]].to_csv(outprefix+"_principal_components.csv", sep="\t")
    return pca_result[pca_result.columns[0:optPCs]]

def calinski_harabasz(matrix, outprefix):
    d = {}
    for k in range(2, 21):
        model = Birch(n_clusters=k, branching_factor=20, threshold=0.5, compute_labels=True).fit(matrix.values)
        labels = model.labels_
        d[k] = metrics.calinski_harabasz_score(matrix.values, labels)
    k = max(d.items(), key=operator.itemgetter(1))[0]
    lists = sorted(d.items())
    x, y = zip(*lists) 
    plt.plot(x, y, 'bx-')
    plt.xlabel('k')
    plt.ylabel('calinski_harabaz_score')
    plt.axvline(x=k, color='k', linestyle='--', label="optK")
    plt.gca().grid(True)
    plt.legend(fontsize=18)
    plt.savefig(outprefix+"_birch_calinski_harabasz.png")
    plt.clf()
    return k


def locate_elbow(matrix, outprefix):
    model = Birch(branching_factor=20, threshold=0.5, compute_labels=True)
    visualizer = KElbowVisualizer(model, k=(2,20))
    visualizer.fit(matrix.values)
    visualizer.finalize()
    plt.savefig(outprefix+"_birch_elbow.png")
    plt.clf()
    return visualizer.elbow_value_


def cluster_and_output(k, matrix, outprefix):
    
    model = Birch(branching_factor=20, n_clusters=int(k), threshold=0.5, compute_labels=True).fit(matrix.values)
    clusters = pd.DataFrame(model.labels_, index=matrix.index, columns=["cluster"])
    clusters.to_csv(outprefix+"_birch_clusters.csv", sep="\t")
    color_palette = sns.color_palette("hls", len(np.unique(model.labels_)))
    cluster_colors = [color_palette[x]
                      for x in model.labels_]
    plt.scatter(matrix.values[:,0],matrix.values[:,1], c=cluster_colors)
    plt.suptitle("Birch clustering result")
    plt.savefig(outprefix+"_birch_scatter.png")
    plt.clf()
    return clusters

def heatmap(cnvs, clusters, pos_list, chrN_list, chr_limits, chr_boundaries, outprefix):
    #clusters+cnv profiles
    cnvs.loc[:,'cluster'] = clusters.cluster

    cnvs = cnvs.sort_values(by='cluster')

    sorted_labels = cnvs['cluster'].values
    color_palette = sns.color_palette("hls", len(np.unique(sorted_labels)))
    cluster_colors = [color_palette[x] 
                      for x in sorted_labels]
    #print(row_colors)

    cnvs = cnvs.drop(['cluster'], axis=1)
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
    plt.savefig(outprefix+'_clusters.heatmap.png')
    plt.clf()


def filter_noisy(cnvs, stats, outprefix):    
    threshold = stats['Disp'].quantile([0.9]).values[0]
    noisy_cells = stats[stats['Disp'] >= threshold].index.values
    #noisy = df_stats[df_stats['Disp'] > threshold].index.values

    sns.distplot(stats['Disp'], rug=True)
    plt.axvline(x=threshold, color='r', linestyle='--', label="90th percentile")
    plt.legend()
    plt.savefig(outprefix + "_noisy.distr.png")
    plt.clf()

    print("N cells: {} - N noisy cells {}".format(len(cnvs), len(noisy_cells)))
    cnvs = cnvs.loc[~cnvs.index.isin(noisy_cells)]
    return noisy_cells, cnvs

    # cell id dictionary
    cell_id_dict = dict(zip(list(range(len(cnvs))), cnvs.index ))
    cell_ids = cnvs.index
    

def main():
    parser = argparse.ArgumentParser(description='Birch clustering')
    parser.add_argument('-f', '--file', type=str, required=True,
                    help='CNV file.')
    parser.add_argument('-s', '--stats', type=str, required=True,
                        help="Stats file.")
    parser.add_argument('-o', '--outprefix', type=str, required=True,
                    help='Output prefix')


    args = parser.parse_args()

    df = pd.read_csv(args.file, sep='\t')
    stats = pd.read_csv(args.stats, sep="\t")
    outprefix = args.outprefix
    
    

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
    
    noisy_cells, filtered_cnvs = filter_noisy(cnvs, stats, outprefix)

    # cell id dictionary
    cell_id_dict = dict(zip(list(range(len(cnvs))), cnvs.index ))
    cell_ids = cnvs.index

    # Low variance feature filtering and PCA
    pca_result = preprocessing(filtered_cnvs, outprefix)
    
    k = locate_elbow(pca_result, outprefix)
    clusters = cluster_and_output(k, pca_result, outprefix)
    heatmap(filtered_cnvs, clusters, pos_list, chrN_list, chr_limits, chr_boundaries, outprefix)


if __name__ == "__main__":
    sys.exit(main())


