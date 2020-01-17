#!/usr/bin/env python


from my_funcs import *
from check_funcs import *

from dicttoxml import dicttoxml
from xml.dom.minidom import parseString
import json
import random
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import matplotlib.colors as colors

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Single-sample analysis.")
    parser.add_argument("sample", metavar='sample_name', action='store',
        help='Sample name.', nargs=1, type=str)
    parser.add_argument("cnvs", metavar='SegCopy', action='store',
        help='Path to cnvs file.', nargs=1, type=str)

    parser.add_argument("results", metavar='results.txt', action='store',
        help='Path to stats file.', nargs=1, type=str)

    parser.add_argument("method", metavar='clust_method', action='store',
        help='Clustering method',
        nargs=1, type=str)

    parser.add_argument("metric", metavar='distance_metric', action='store',
        help='Distance metric',
        nargs=1, type=str)

    parser.add_argument("meta_format", choices=['json', 'xml'], action='store',
        help='Metadata file format.',
        nargs=1, type=str)

    parser.add_argument("outdir", metavar='outdir', action='store',
        help='Path to the desired output directory where the merged files have to be stored',
        nargs=1, type=str)

    parser.add_argument("--seed", metavar='N', action='store',
        help='Seed to initialize the pseudo-random generator used to perform the permutation test.',
        nargs=1, type=int)

    parser.add_argument("--n_permutations", metavar='N', action='store',
        help='Number of permutations to execute the permutation test for sample coesion score.',
        nargs=1, type=int)

    #parser.add_argument("--tsne_iterations", metavar='N', action='store', help='Number of iterations for tSNE computation.', nargs=1, type=int)

    #parser.add_argument("--tsne_perplexity", metavar='N', action='store', help='Perplexity value for tSNE computation.', nargs=1, type=int)

    parser.add_argument("--reclust", metavar='n_clusters', action=check_valid(),
        help='If this option is specified, only the clustering part is executed with the specified number of clusters, unless --reinit option is specified (see below).',
        nargs=1, type=int)

    parser.add_argument("--reinit", action='store_true',
        help='This option has effect only if combined with the --clustering option. It allows to recompute the entire analysis and then recluster with the specified number of clusters.')

    parser.add_argument("--verbose", action='store_true', help='Verbose execution.')


    args=parser.parse_args()


    sample = args.sample[0]
    cnvsf = args.cnvs[0]
    results = args.results[0]
    method = args.method[0]
    metric = args.metric[0]
    meta_format = args.meta_format[0]

    outdir = args.outdir[0]

    # default params
    n_permutations = 10
    #tsne_iterations = 5000

    #tnse_perplexity is computed afterwards if not specified as parameter

    reclust = False
    reinit = False
    verbose = False

    if args.n_permutations:
        n_permutations = args.n_permutations[0]

    if args.seed:
        random.seed(args.seed[0])
    '''
    if args.tsne_iterations:
        tsne_iterations = args.tsne_iterations[0]
    
    if args.tsne_perplexity:
        perplexity = args.tsne_perplexity[0]
    '''
    if args.reclust:
        reclust = True
        n_clusters = args.reclust[0]
        if args.reinit:
            reinit = True
    if args.verbose:
        verbose = True

    metadata = {}
    metadata['sample'] = sample
    if not reclust:
        metadata['analysis'] = 'complete'
    else:
        metadata['analysis'] = 'reclust'

    metadata['clustering_method'] = method
    metadata['distance_metric'] = metric
    metadata['n_permutations'] = n_permutations
    #metadata['tsne_iterations'] = tsne_iterations
    #if args.tsne_perplexity:
    #    metadata['tsne_perplexity'] = perplexity
    #else:
    #    metadata['tsne_perplexity'] = ''
    
    if reclust:
        metadata['n_clusters'] = n_clusters
    else:
        metadata['n_clusters'] = ''

    if args.seed:
        metadata['seed'] = seed
    else:
         metadata['seed'] = ''



    df = pd.read_csv(cnvsf, sep="\t", usecols = lambda column : column not in ['Unnamed: 103', 'Unnamed: 113', 'Unnamed: 32'])

    cnvs = df.drop(['CHR', 'START', 'END'], axis=1).transpose()
    boundaries = df[['CHR', 'START', 'END']].copy()

    if reclust and not reinit:
        '''
            Recompute  only Z linkage matrix

        '''
        print_msg("Reclustering cells", 0, verbose)
        Z = linkage(cnvs, method=method, metric=metric)
    else:

        '''
            heatmap
        '''
        print_msg("Complete analysis", 0, verbose)
        print_line(verbose)
        print_msg("Computing heatmap and phylogenetic tree (method = %s, metric = %s)"%(method, metric), 0, verbose)
        Z = heatmap(cnvs, boundaries, method, metric, outdir, False, verbose, sample=sample)

        '''
            tSNE
        
        print_line(verbose)
        if not args.tsne_perplexity:
            tsne_pca(cnvs = cnvs, tsne_iterations=tsne_iterations, metadata=metadata, outdir=outdir, verbose=verbose)
        else:
            tsne_pca(cnvs = cnvs, perplexity=perplexity, tsne_iterations=tsne_iterations, metadata=metadata, outdir=outdir, verbose=verbose)
        '''
        ''' 
            Plot distribution of mean ploidies over all cells
        '''
        print_line(verbose)
        print_msg("Plotting mean ploidy distribution", 0, verbose)
        
        keyword = 'Copy_Number'
        res_df = pd.read_csv(results, sep="\t").dropna()
        min_cn = np.min(res_df[keyword].values)
        max_cn = np.max(res_df[keyword].values)
        

        res_df[keyword].plot.kde(bw_method=0.05, grid=True)
        
        plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        plt.xticks(np.arange(0, max_cn+1, 0.1))
        plt.xlim(0, max_cn+1)
        #plt.ylim(-0.2, 7)
        plt.autoscale(enable=True, axis='y', tight=True)
        
        plt.gca().set_xlabel("Mean ploidies")
        plt.gca().set_ylabel("Ploidies distribution")

        plt.gcf().suptitle("Distribution plot of mean ploidies of all cells", fontsize=14, fontweight='bold')
        plt.gcf().set_size_inches(36,14)
        plt.savefig(outdir+"/mean_plody_distribution.svg")
        plt.gcf().clf()

        '''
            Compute and plot mean CNV profile

        '''
        print_line(verbose)
        print_msg("Computing mean CNV profile", 0, verbose)
        mean_df = cnvs.mean(axis=0)
        mean_ploidy = mean_df.mean()

        print_msg("mean ploidy = %f"%(mean_ploidy), 1, verbose)

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



        fig, ax = plt.subplots()
        fig.set_size_inches(36,10)
        mean_df.plot(legend=False, ax=ax)
        ax.set_title('Sample ' + sample + ', mean ploidy =' + str(mean_ploidy))
        ax.set_xlabel("Chromosomes", fontweight='bold', fontsize=14)
        ax.set_ylabel("Mean copy number", fontweight='bold', fontsize=14)

        #place vertical lines to identify chromosomes
        for pos in chr_limits:
            ax.axvline(x=pos, color='black')

        #place chromosome ticks at the right position
        ax.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter((chrN_list)))
        ax.tick_params(axis='x', labelsize=14)

        ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_boundaries))
        ax.tick_params(axis='x', length=20, which='minor')
        
        fig.suptitle("Mean cnv profiles of sample = " + sample, fontsize=16, fontweight='bold')
        fig.savefig(outdir+"/mean_cnv.png")
        fig.clf()


        '''
            Silhouette coefficient: optimal number of clusters
        '''
        print_line(verbose)
        n_clusters = my_silhouette_score(Z, cnvs, outdir, verbose, metadata)
    '''
        Extract clusters and compute mean cnv profiles

    '''
    print_line(verbose)

    df_mean = extract_clusters(Z, boundaries, n_clusters, cnvs, reclust, outdir, verbose)
    heatmap(df_mean.transpose(), boundaries, method, metric,outdir, True, verbose, sample=sample)

    if meta_format == 'xml':
        xml=dicttoxml(metadata,  custom_root='analysis', attr_type=False)
        dom = parseString(xml)
        with open(outdir+'/metadata.xml', 'w+') as f:
            f.write(dom.toprettyxml())
    else:
        with open(outdir+'/metadata.json', 'w+') as f:
            json.dump(metadata, f, indent=2)
