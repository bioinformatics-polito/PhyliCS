import sys
import random
import ntpath
import argparse
import operator
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib 
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.transforms as transforms
from scipy.cluster import hierarchy
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_samples, silhouette_score
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cophenet
from scipy.spatial.distance import pdist, cdist, squareform


def multi_sample_heatmap(cnvs, boundaries, samples_dict, method, metric, outdir, verbose):
    #compute the ending position of each chromosome to draw vertical lines on the plot
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
    
    merged_cnvs = cnvs.copy()
    merged_cnvs['sample'] = ''

    for idx, row in merged_cnvs.iterrows():
        merged_cnvs.loc[idx, 'sample'] = idx.split(":")[0]
    #print(merged_cnvs['sample'])
    #multilevel_index
    merged_cnvs.set_index([merged_cnvs.index, 'sample'], inplace=True)
    samples = samples_dict.keys()
    # Create a categorical palette to identify the samples
    sample_pal = sns.color_palette("bright", len(samples_dict.keys()))
    #sample_pal = sns.husl_palette(len(samples), h=.5)
    sample_lut = dict(zip(map(str, samples), sample_pal))
    # Convert the palette to vectors that will be drawn on the side of the matrix
    allsamples = merged_cnvs.index.get_level_values('sample')

    merged_cnvs = merged_cnvs.reset_index(level=1, drop=True)
    #print(merged_cnvs)
    sample_colors=pd.Series(allsamples, index=merged_cnvs.index).map(sample_lut)
    #print(sample_colors)
    cbar_kws={"ticks":np.arange(0,7,1)}


    h = sns.clustermap(merged_cnvs, method=method, metric=metric,
        col_cluster=False,
        yticklabels = False,
        row_colors=sample_colors,
        cmap='RdBu_r',
        vmin=0, vmax=6,
        center=2,
        #norm=divnorm,
        cbar_kws=cbar_kws)

    # Draw the legend bar for the classes                 
    for label in allsamples.unique():
        h.ax_row_dendrogram.bar(0, 0, color=sample_lut[label],
                            label=label, linewidth=0)
    print(h.ax_row_dendrogram)
    h.ax_row_dendrogram.legend(title='Samples', loc='center', bbox_to_anchor=(0.5,0.9), ncol=len(samples), fontsize=14, title_fontsize='x-large', bbox_transform=plt.gcf().transFigure)
    h.cax.set_position([0.05, .2, .03, .45])
    
    Z = h.dendrogram_row.linkage
    c, coph_dist = cophenet(Z,pdist(merged_cnvs, metric))
    textstr = '(cophenet coeff= ' + str(c)+')'

    print_msg("Cophenet coefficient: " + str(c), 1, verbose)

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
    
    plt.gcf().suptitle('Multi-sample heatmap\n' + textstr, fontsize=16, fontweight='bold' )
    plt.gcf().set_size_inches(37, 21)
    plt.savefig(outdir+'/multi_sample_heatmap.png')
    plt.clf()

    return Z

def heatmap(cnvs, boundaries, method, metric, outdir, clusters, verbose, sample=None):
    divnorm = colors.DivergingNorm(vmin=0, vcenter=2, vmax=12)
    
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

    if clusters:
        yticklabels = True
    else:
        yticklabels = False

    cbar_kws={"ticks":np.arange(0,13,1)}
    h = sns.clustermap(cnvs, method=method, metric=metric, col_cluster=False, yticklabels = yticklabels,  cmap='RdBu_r', vmin=0, vmax=12,norm=divnorm, cbar_kws=cbar_kws)
    Z = h.dendrogram_row.linkage
    
    if clusters == False:
        h.cax.set_position([0.05, .2, .03, .45])
        c, coph_dist = cophenet(Z,pdist(cnvs, metric))
        textstr = '(cophenet coeff= ' + str(c)+')'
        print_msg("cophenet coefficient: " + str(c), 1, verbose)


    ax = h.ax_heatmap
    #place vertical lines to identify chromosomes
    for pos in chr_limits:
        ax.axvline(x=pos, color='black')

    #place chromosome ticks at the right position
    ax.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter((chrN_list)))
    ax.tick_params(axis='x', rotation=0, labelsize=14)
    
    ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_boundaries))
    ax.tick_params(axis='x', length=20, which='minor')

    ax.set_xlabel("Chromosomes", fontweight='bold', fontsize=14)
    if clusters:
        ax.set_ylabel("Clusters", fontweight='bold', fontsize=14)
    else:
        ax.set_ylabel("cells", fontsize=14, fontweight='bold')

    plt.gcf().set_size_inches(37, 21)
    
    if clusters:
        plt.gcf().suptitle("Clusters mean CNV heatmap  " + sample, fontsize=16, fontweight='bold')
        plt.savefig(outdir+"/clusters_heatmap.png")
    else:
        plt.gcf().suptitle("CNV heatmap - Sample = " + sample + '\n' + textstr, fontsize=16, fontweight='bold')
        plt.savefig(outdir+"/heatmap.png")
    plt.clf()
    return Z


def tsne_pca(cnvs, tsne_iterations, metadata, outdir, verbose, perplexity=None):
    #preliminary dimensionality reduction with pca
    #min(n_samples, n_features)
    if perplexity == None:
            #suggested perplexity values are in the range [5,50] and
        #must be < number of points
        perplexity = min(len(cnvs)/2, 40)
    metadata['tsne_perplexity'] = perplexity
    print_msg("Computing tSNE (perplexity = %d, n_iterations = %d)"%(perplexity, tsne_iterations), 0, verbose)
    pca = PCA(n_components=min(len(cnvs), len(cnvs.columns)))
    pca_result = pca.fit_transform(cnvs)
    np.savetxt(outdir+"/varexp.tsv", pca.explained_variance_ratio_, delimiter="\t")

    #tsne

    tsne = TSNE(n_components=2, verbose=0, perplexity=perplexity, n_iter=tsne_iterations)
    tsne_results = tsne.fit_transform(pca_result) # con le pc
    dpca = pd.DataFrame({'tsne1':tsne_results[:,0],'tsne2':tsne_results[:,1]})
    plt.figure(figsize=(16,10))
    sns.scatterplot(x="tsne1", y="tsne2",data=dpca, legend="full")
    plt.suptitle("tSNE analysis")
    #plt.switch_backend('agg') # trick che si era reso necessario per salvare png su hactar (mi pare)
    plt.savefig(outdir+"/tsne.png")
    plt.gcf().clf()

def my_silhouette_score(Z, cnvs, outdir, verbose, metadata):
    print_msg("Computing the optimal number of clusters", 0, verbose)
    range_n_clusters = range(2,18)
    silhouette_avgs = {}
    i = 0
    j = 0
    NR = 4
    NC = 4

    fig, ax = plt.subplots(NR, NC)
    fig.set_size_inches(46, 20)
    silhouette_df = pd.DataFrame(columns=['n_clusters','avg silhouette score'])
    for n_clusters in range_n_clusters:
        # The silhouette coefficient can range from -1, 1
        ax[i][j].set_xlim([-1, 1])
        # The (n_clusters+1)*10 is for inserting blank space between silhouette
        # plots of individual clusters, to demarcate them clearly.
        ax[i][j].set_ylim([0, len(cnvs) + (n_clusters + 1) * 10])

        #Z = h.dendrogram_row.linkage
        cluster_labels = fcluster(Z, t = n_clusters, criterion = 'maxclust')

        # The silhouette_score gives the average value for all the samples.
        # This gives a perspective into the density and separation of the formed
        # clusters
        silhouette_avg = silhouette_score(cnvs, cluster_labels)
        silhouette_avgs[n_clusters] = silhouette_avg
        print_msg("n_clusters = " + str(n_clusters) + " - The average silhouette_score is :" + str(silhouette_avg), 1, verbose)
        silhouette_df = silhouette_df.append(pd.DataFrame([[n_clusters, silhouette_avg]], columns=['n_clusters','avg silhouette score']), ignore_index=True)
        # Compute the silhouette scores for each sample
        sample_silhouette_values = silhouette_samples(cnvs, cluster_labels)

        y_lower = 10
        for k in range(1, n_clusters+1):
            # Aggregate the silhouette scores for samples belonging to
            # cluster k, and sort them
            #print('k = ' + str(n_clusters) + ' \n\t cluster = ' + str(i) + '\n\t n_items = ' + str(len(sample_silhouette_values[cluster_labels == i])))
            kth_cluster_silhouette_values = \
                sample_silhouette_values[cluster_labels == k]

            kth_cluster_silhouette_values.sort()

            size_cluster_k = kth_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_k
            color = cm.nipy_spectral(float(k) / n_clusters)
            ax[i][j].fill_betweenx(np.arange(y_lower, y_upper),
                          0, kth_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

            # Label the silhouette plots with their cluster numbers at the middle
            ax[i][j].text(-0.05, y_lower + 0.5 * size_cluster_k, str(k))

            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples


            ax[i][j].set_xlabel("The silhouette coefficient values")
            ax[i][j].set_ylabel("Cluster label")

            # The vertical line for average silhouette score of all the values
            ax[i][j].axvline(x=silhouette_avg, color="red", linestyle="--")

            ax[i][j].set_yticks([])  # Clear the yaxis labels / ticks
            ax[i][j].set_xticks([-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
            ax[i][j].set_title(("n_clusters = %d" % n_clusters),fontweight='bold')

        if j < NC-1:
            j = j+1
        else:
            #columns iterator reached the last column:
            #update row iterator and reset column one
            j = 0
            i = i + 1


    plt.suptitle("Silhouette analysis for hierarchical clustering", fontsize=14, fontweight='bold')
    fig.savefig(outdir+"/silhouette_results.png")
    fig.clf()


    silhouette_df.to_csv(outdir+'/per_k_silhouette_scores.csv', index=False, sep='\t')

    opt_n_clusters, max_silhouette_avg =  max(silhouette_avgs.items(),key=operator.itemgetter(1))
    worst_n_clusters, min_silhouette_avg =  min(silhouette_avgs.items(),key=operator.itemgetter(1))
    print_msg("The maximum average silhouette_score is: " + str(max_silhouette_avg) + "- The suggested number of clusters is = " + str(opt_n_clusters), 1, verbose)
    fig2, ax1 = plt.subplots()
    ax1.set_ylim([min_silhouette_avg-0.2, max_silhouette_avg+0.2])

    ax1.plot(range_n_clusters, list(silhouette_avgs.values()), 'o')
    ax1.axvline(x=opt_n_clusters, color="red", linestyle="--")
    ax1.grid(True)

    ax1.set_xlabel("Number of clusters")
    ax1.set_ylabel("Average silouette score")
    fig2.suptitle("Silhouette analysis summary for clustering", fontsize=14, fontweight='bold')

    fig2.savefig(outdir+"/silhouette_summary.png")
    fig2.clf()

    n_clusters = opt_n_clusters
    metadata['n_clusters'] = n_clusters

    return n_clusters

def extract_clusters(Z, boundaries, n_clusters, cnvs, reclust, outdir, verbose):
    print_msg("Computing mean cnv profiles for clusters", 0, verbose)
    #Z = h.dendrogram_row.linkage
    
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

    
    cluster_labels = fcluster(Z, t=n_clusters, criterion='maxclust')

    if reclust:
         silhouette_avg = silhouette_score(cnvs, cluster_labels)

    cnvs['CLUSTER'] = cluster_labels
    cnvs.index.name = 'cellid'
    #cnvs['CLUSTER'].to_csv(outdir+"/clusters.csv", sep='\t', header=True)
    clusters = cnvs['CLUSTER'].unique()

    '''
    Create a new dataframe containing the mean cnvs for each cluster
    '''
    df_mean = pd.DataFrame(columns=clusters)
    df_results = pd.DataFrame(columns=['cluster', 'mean_ploidy', 'cells'])
    fig, axes = plt.subplots(n_clusters, 1)
    fig.set_size_inches(36, 10*n_clusters)
    i = 0
    for cluster in clusters:
        df_mean[cluster] = cnvs.loc[cnvs['CLUSTER'] == cluster].drop(['CLUSTER'], axis=1).mean().values
        cells = cnvs.loc[cnvs['CLUSTER'] == cluster].index.values
        mean_ploidy = df_mean[cluster].mean()
        df_clone = pd.DataFrame([[cluster, mean_ploidy, ','.join( cells )]], columns=['cluster', 'mean_ploidy', 'cells'])
        df_results = df_results.append(df_clone, ignore_index=True)

        print_msg("cluster %d, mean ploidy = %f"%(cluster,mean_ploidy), 1, verbose)

        df_mean.plot(y=cluster,  legend=False, ax=axes[i])
        axes[i].set_title('Cluster ' + str(cluster) + ', mean ploidy =' + str(mean_ploidy))
        axes[i].set_xlabel("Chromosomes", fontweight='bold', fontsize=14)
        axes[i].set_ylabel("Mean copy number", fontweight='bold', fontsize=14)

        #place vertical lines to identify chromosomes
        for pos in chr_limits:
            axes[i].axvline(x=pos, color='black')

        #place chromosome ticks at the right position
        axes[i].xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
        axes[i].xaxis.set_major_formatter(ticker.FixedFormatter((chrN_list)))
        axes[i].tick_params(axis='x', labelsize=14)

        axes[i].xaxis.set_minor_locator(ticker.FixedLocator(chr_boundaries))
        axes[i].tick_params(axis='x', length=20, which='minor')

        i = i+1


    fig.suptitle("Mean cnv profiles of clusters", fontsize=16, fontweight='bold')
    fig.savefig(outdir+"/mean_cnv_clusters.png", )
    fig.clf()

    #df_mean.to_csv(outdir+'/multi_sample_SegCopy_clusters.csv', sep='\t')
    df_results.to_csv(outdir+'/clusters.tsv', sep='\t', index=False)

    return df_mean

"""
def check_valid_sample(nmin):
    class CheckValid(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin<=len(values):
                msg='argument "{f}" requires at least {nmin} arguments'.format(
                        f=self.dest,nmin=nmin)
                raise argparse.ArgumentTypeError(msg)
            for value in values:
                arguments = value.split(":")
                if not len(arguments) == 2:
                    msg='argument "{f}" must be formatted as follows: sample:SegCopy'.format(
                            f=self.dest)
                    raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return CheckValid

def check_valid_N_clust():
    class CheckValid(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if values[0] <= 0:
                msg='argument "{f}" requires a number of clusters greater than 0'.format(
                        		f=self.dest)
               	raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return CheckValid

def check_valid():
    class CheckValid(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if values[0] <= 0:
                msg='argument "{f}" requires a number of clusters greater than 0'.format(
                                        f=self.dest)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return CheckValid
"""
def findsubsets(S,n):
    return set(itertools.combinations(S, n))

def partition (list_in, n):
    random.shuffle(list_in)
    return [list_in[i::n] for i in range(n)]

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)
    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def print_msg(msg, level, verbose):
	if verbose:    
            print("[" + path_leaf(sys.argv[0]) +"] " + ("-"*2*level) + " " + msg)

def print_line(verbose):
    if verbose:
        print("-"*80)
