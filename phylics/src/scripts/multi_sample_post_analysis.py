#!/usr/bin/env python

# ==========================================================================
#                                  PhyliCS
# ==========================================================================
# This file is part of PhyliCS.
#
# TOOL is Free Software: you can redistribute it and/or modify it
# under the terms found in the LICENSE.rst file distributed
# together with this file.
#
# PhyliCS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# ==========================================================================
# Author: Marilisa Montemurro <marilisa.montemurro@polito.it>
# ==========================================================================
# multi_sample_post_analysis.py: Multi-sample analysis module
# ==========================================================================


from funcs import *
from check_funcs import *
#from dicttoxml import dicttoxml
#from xml.dom.minidom import parseString
import time
import json
import random
import warnings
import operator
import argparse
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
from joblib import Parallel, delayed
from sklearn.metrics import silhouette_samples, silhouette_score
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cophenet
from scipy.spatial.distance import pdist, cdist, squareform

def silhouette_(df, labels, i):
    if i%100 == 0:
        print_msg("iteration: %d"%i, 2, verbose)
    random.shuffle(labels)
    return silhouette_score(df, labels)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Multi-sample analysis.")
    parser.add_argument("samples", metavar='sample_name:SegCopy', 
            help='Sample name and cnvs filepath separated by ":". At least two samples must be provided.',
            nargs='+', type=str)

    parser.add_argument("method", metavar='clust_method', action='store',
            help='Clustering method.',
            nargs=1, type=str)

    parser.add_argument("metric", metavar='distance_metric', action='store',
            help='Distance metric.',
            nargs=1, type=str)

    """
    parser.add_argument("meta_format", choices=['json', 'xml'], action='store',
            help='Metadata file format.',
            nargs=1, type=str)

    """
    parser.add_argument("outdir", metavar='out_dir', action='store',
            help='Path to the output directory (must exist).',
            nargs=1, type=str)

    parser.add_argument("--n_permutations", metavar='N', action='store',
            help='Number of permutations to execute the permutation test for sample coesion score.',
            nargs=1, type=int)

    parser.add_argument("--seed", metavar='N', action='store',
            help='Seed to initialize the pseudo-random generator used to perform the permutation test.',
            nargs=1, type=int)

    '''
    parser.add_argument("--tsne_iterations", metavar='N', action='store',
            help='Number of iterations for tSNE computation.',
            nargs=1, type=int)

    parser.add_argument("--tsne_perplexity", metavar='N', action='store',
            help='Perplexity value for tSNE computation.',
            nargs=1, type=int)
    '''
    parser.add_argument("--reclust", metavar='N', 
            help='If this option is specified, only the clustering part is executed with the specified number of clusters, unless --reinit option is specified (see below).',
            nargs=1, type=int)

    parser.add_argument("--reinit", action='store_true',
            help='This option has effect only if combined with the --clustering option. It allows to recompute the entire analysis and then recluster with the specified number of clusters.')
    parser.add_argument("--verbose", action='store_true', help='Verbose execution.')

    parser.add_argument("--n_jobs", metavar='NJ', action='store', type=int, help='Number of parallel jobs to speed up pvalue computation')

    args=parser.parse_args()

    samples = args.samples
    method = args.method[0]
    metric = args.metric[0]
    #meta_format = args.meta_format[0]
    outdir = args.outdir[0]

    # default params
    n_jobs = 1
    n_permutations = 10
    #tsne_iterations = 5000

    #tnse_perplexity is computed afterwards if not specified as parameter

    reclust = False
    reinit = False
    verbose = False

    if args.n_jobs:
        n_jobs = args.n_jobs
    
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

    """
    metadata = {}
    if not reclust:
        metadata['analysis'] = 'complete'
    else:
        metadata['analysis'] = 'reclust'

    metadata['samples'] = samples
    metadata['clustering_method'] = method
    metadata['distance_metric'] = metric
    metadata['n_permutations'] = n_permutations
    #metadata['tsne_iterations'] = tsne_iterations
    '''
    if args.tsne_perplexity:
        metadata['tsne_perplexity'] = perplexity
    else:
        metadata['tsne_perplexity'] = ''
    '''
    if reclust:
        metadata['n_clusters'] = n_clusters
    else:
         metadata['n_clusters'] = ''

    if args.seed:
        metadata['seed'] = seed
    else:
         metadata['seed'] = ''
    """

    cnvs = pd.DataFrame()
    samples_dict = {}
    cell_ids = {}
    sample_count = 0

    '''
        CNV calls merging 
    '''

    print_msg( "CNV calls merging",0, verbose)
    for sample in samples:
        sample_name = sample.split(":")[0]
        segcopy = sample.split(":")[1]
        
        samples_dict[sample_name] = {}

        df = pd.read_csv(segcopy, sep="\t",  usecols = lambda column : column not in ['Unnamed: 103', 'Unnamed: 113', 'Unnamed: 32'])

        #if this is the first sample copy the header columns
        #to build  the cnvs merged dataframe
        if sample_count == 0:
            cnvs[['CHR', 'START', 'END']] = df[['CHR', 'START', 'END']]
        
        df = df.drop(['CHR', 'START', 'END'], axis=1)
        columns = df.columns
        new_columns =[]
        
        #change the column names to reidentify them
        #when merged
        for col in columns:
            new_columns.append(sample_name + ":" + col)
        cnvs[new_columns] = df[columns]

    cnvs.to_csv(outdir+'/SegCopy_merged', sep='\t', index=False)

    merged_cnvs = cnvs.drop(['CHR', 'START', 'END'], axis=1).transpose()
    boundaries = cnvs[['CHR', 'START', 'END']]

    print_line(verbose)

    if reclust and not reinit:
        '''
             Recompute  only Z linkage matrix
        '''
        print_msg("Reclustering cells", 0, verbose)
        Z = linkage(merged_cnvs, method=method, metric=metric)
    else:    
        '''
        Compute the average silhouette score by assigning cells to fixed clusters (their samples)
        to evaluate sample coesion with respect to samIple separation
        '''
        print_msg("Complete analysis", 0, verbose)
        
        print_line(verbose)
        print_msg("Heterogeneity score computation", 0, verbose)
        
        coesion_df = pd.DataFrame(columns=['samples', 'het_score', 'pvalue', 'n_permutations'])
        
        coesion_score  = {}

        sample_labels = []
        allcells = merged_cnvs.index.tolist()
        for cell in allcells:
            label = cell.split(':')[0]
            sample_labels.append(label)

        silhouette_avg = silhouette_score(merged_cnvs, sample_labels)
        coesion_score[str(list(set(sample_labels)))] = silhouette_avg
        '''
            permutation test: how many times we obtain a 'more extreme' (<=,=>) silhouette if we randomly shuffle the sample labels.
        '''
        #n_samples = sample_count
        silhouettes = []

        #print_msg("Permutation test (n_permutations = %d)"%n_permutations, 1, verbose)
       
        t0 = time.time()
        
        #for i in range(n_permutations):
           #print_msg("iteration {} --- shuffling".format(i), 2, verbose)
            #random.shuffle(sample_labels)
            #silhouettes.append(silhouette_score(merged_cnvs, sample_labels))
        silhouettes = Parallel(n_jobs=n_jobs, prefer="threads")(
                delayed(silhouette_)(merged_cnvs, sample_labels, i) for i in range(n_permutations))
        t = time.time() - t0
        print_msg("Permutation test done ({}s)".format(t), 2, verbose)

        
        #bilateral test
        pvalue = 0

        if silhouette_avg >= 0:
            for s in silhouettes:
                if s >= silhouette_avg:
                    pvalue += 1
        else:
            for s in silhouettes:
                if s <= silhouette_avg:
                    pvalue += 1

        pvalue = pvalue/len(silhouettes)
        coesion_df = coesion_df.append(pd.DataFrame([[str(list(set(sample_labels))), silhouette_avg, pvalue, n_permutations]], columns=['samples', 'het_score', 'pvalue', 'n_permutations']), ignore_index=True)
        print(coesion_df)
        print_msg("Samples: " + str(list(set(sample_labels))), 1, verbose)
        print_msg("The average heterogenity score is: %f (pvalue=%.4f)"%(silhouette_avg, pvalue), 1, verbose)

        '''
            Create a boxplot showing the silhouette score distribution
        '''

        warnings.filterwarnings("ignore")
        fig, ax = plt.subplots()
        #fig.set_size_inches(12, 5)
        s_df = pd.DataFrame(silhouettes, columns=['s'])
        sns.violinplot(y='s', data=s_df)
        ax.set_ylabel("heterogeneity score", fontsize=12, fontweight='bold')
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        
        l = ax.axhline(y=silhouette_avg, color='red', ls='--')
        l.set_label('Observed score')
        trans = transforms.blended_transform_factory(
            ax.get_yticklabels()[0].get_transform(), ax.transData)
        ax.text(0,silhouette_avg, "{:.4f}".format(silhouette_avg), color="red", transform=trans,
            ha="right", va="center")

        plt.legend()

        fig.suptitle("Average heterogeneity score: permutation test distribution.\nSamples = {}".format(str(list(set(sample_labels)))), fontsize=14, fontweight='bold')
        fig.savefig(outdir+"/het_violin_boxplot.png", bbox_inches='tight')
        fig.clf()

        '''
            if n_samples > 2, compute sample coesion score, for each sample, with respect to a super-cluster
            made of the cells of all other cells. 

            If the score improves with respect to the global score, cells are better assigned to these clusters,I
            which means that this sample is highly separate from the rest of the samples, which should be assigned
            to a single cluster, rather than to different ones (the other samples are very similar with respect to the 
            other one). On the contrary, if it enworsens, this means that cells are better assigned to previous ones.
        '''

        sample_names = list(set(sample_labels))
        n_samples = len(sample_names)
        counter = 1
        if n_samples > 2:
            print_line(verbose)
            print_msg("Sample-aggregate heterogeneity score",0, verbose)
            max_size = n_samples-1
            if n_samples%2 == 0: #even
                min_size = n_samples//2 -1
            else: #odd
                min_size = n_samples//2
            for N in range(max_size, min_size, -1):
                #generate all possible sample aggregation to 
                #evaluate if any of them improves the silhouette score
                superclusters = findsubsets(sample_names, N)
                #e.g. samples= ['a', 'b', 'c', 'c'], superclusters of N = 3 {('a', 'b', 'd'), ('a', 'c', 'd'), ('b', 'c', 'd'), ('a', 'b', 'c')}
                for supercluster in superclusters:
                    #iterate on all possible aggregations
                    #e.g. ('a', 'b', 'd')
                    label1 = ""
                    label2 = ""
                    new_labels = []
                    for c in supercluster:
                        if label1 == "":
                            label1 = c
                        else:
                            label1 += "+"+c #e.g. label1 = 'a+b+d"

                    for c in sample_names:
                        if c not in supercluster:
                            if label2 == "":
                                label2 = c
                            else:
                                label2 += "+"+c #e.g. label2 = 'c'

                    for cell in allcells:
                        sample = cell.split(':')[0]
                        if sample in supercluster: 
                            new_labels.append(label1) #assign label1 to the samples belonging to the supercluster
                        else:
                            new_labels.append(label2) #assign label2 to the samples belonging to the other supercluster
                
                    print_msg("Heterogeneity computation for sample aggregates", 1, verbose)
                    print_msg("Samples = %s,%s"%(label1, label2), 1, verbose)
                    silhouette_avg = silhouette_score(merged_cnvs, new_labels)
                    coesion_score[str(list(set(new_labels)))] = silhouette_avg
                    #permutation test
                    
                    silhouettes = []
                    #n_iterations = 1000
                    
                    
                    print_msg("Permutation test (n_permutations = %d)"%n_permutations, 1, verbose)
                    t0 = time.time()
                    ''''
                    for i in range(n_permutations):
                        random.shuffle(new_labels)
                        silhouettes.append(silhouette_score(merged_cnvs, new_labels))
                        if i%100 == 0:
                            print_msg("iteration: %d"%i, 2, verbose)
                    print_msg("iteration: %d"%i, 2, verbose)
                    '''
                    silhouettes = Parallel(n_jobs=n_jobs, prefer="threads")(
                                            delayed(silhouette_)(merged_cnvs, new_labels, i) for i in range(n_permutations))

                    t = time.time() -t0
                    print_msg("Permutation test done ({}s)".format(t), 2, verbose)
                    
                    #bilateral test
                    pvalue = 0
                    
                    if silhouette_avg >= 0:
                        for s in silhouettes:
                            if s >= silhouette_avg:
                                pvalue += 1
                    else:
                        for s in silhouettes:
                            if s <= silhouette_avg:
                                pvalue += 1
                    
                    pvalue = pvalue/len(silhouettes)
                    coesion_df = coesion_df.append(pd.DataFrame([[str(list(set(new_labels))), silhouette_avg, pvalue, n_permutations]], columns=['samples', 'het_score', 'pvalue', 'n_permutations']), ignore_index=True)

                    print_msg("The average heterogeneity score is: %f (pvalue=%.4f)"%(silhouette_avg, pvalue), 1, verbose)

                    warnings.filterwarnings("ignore")
                    fig, ax = plt.subplots()
                    #fig.set_size_inches(12, 5)
                    s_df = pd.DataFrame(silhouettes, columns=['s'])
                    sns.violinplot(y='s', data=s_df)
                    ax.set_ylabel("heterogeneity score", fontsize=12, fontweight='bold')
                    ax.get_xaxis().tick_bottom()
                    ax.get_yaxis().tick_left()

                    l = ax.axhline(y=silhouette_avg, color='red', ls='--')
                    l.set_label('Observed score')
                    trans = transforms.blended_transform_factory(
                        ax.get_yticklabels()[0].get_transform(), ax.transData)
                    ax.text(0,silhouette_avg, "{:.4f}".format(silhouette_avg), color="red", transform=trans,
                        ha="right", va="center")

                    plt.legend()

                    fig.suptitle("Average heterogeneity score: permutation test distribution.\nSamples = {}, {}".format(label1, label2), fontsize=14, fontweight='bold')
                    fig.savefig(outdir+"/het_score_boxplot_{}.png".format(counter), bbox_inches='tight')
                    fig.clf()
                    
                    counter = counter + 1

            opt_division, max_coesion_score =  max(coesion_score.items(),key=operator.itemgetter(1))
            
            print_msg("Maximum heterogeneity score: %f, clusters: %s"%(max_coesion_score, opt_division),1, verbose)
            coesion_df.to_csv(outdir+'/heterogeneity_scores.csv', index=False, sep='\t')
            #print(coesion_score.keys())

            exit(0) 

            fig2, ax1 = plt.subplots()
            fig2.set_size_inches(18, 7)
            
            cmap = cm.get_cmap("Rainbow")
            colors = cmap(np.linspace(0, 1, len(coesion_score.keys())))
            for i in range(len(coesion_score.keys())):
                ax1.plot(list(coesion_score.keys())[i], list(coesion_score.values())[i], 'o', color=colors[i], label=list(coesion_score.keys())[i])


            for i in range(len(coesion_score.keys())):
                txt =  "{:.4f}".format(list(coesion_score.values())[i])
                ax1.annotate(txt, xy= (i+0.003, list(coesion_score.values())[i]+0.003), fontsize=14)
                
                #ax1.axhline(xmax=i, y=list(coesion_score.values())[i], color='blue', ls='--')
            

            
            ax1.grid(True)
            ax1.set_xlabel("Samples")
            ax1.set_ylabel("Average heterogeneity score")
            fig2.suptitle("Heterogeneity scores for different sample aggragation", fontsize=14, fontweight='bold')        
            fig2.savefig(outdir+"/multi_aggregation_het_score.png")
            fig2.clf()
        coesion_df.to_csv(outdir+'/heterogeneity_scores.csv', index=False, sep='\t')
        '''
            Heatmap
        '''
        print_line(verbose)
        print_msg("Computing heatmap and phylogenetic tree (method = %s, metric = %s)"%(method, metric), 0, verbose)
        #sample label to visualize different colors for different samples in the
        # merged heatmap
        Z = multi_sample_heatmap(merged_cnvs, boundaries, samples_dict, method, metric, outdir, verbose)
        #print(merged_cnvs)
        
        '''
            tSNE
        '''
        '''
        #def tsne_pca(cnvs, perplexity=None, tsne_iterations, metadata, outdir, verbose)
        print_line(verbose)
        if not args.tsne_perplexity:
            tsne_pca(cnvs = merged_cnvs, tsne_iterations=tsne_iterations, metadata=metadata, outdir=outdir, verbose=verbose)
        else:
            tsne_pca(cnvs = merged_cnvs, perplexity=perplexity, tsne_iterations=tsne_iterations, metadata=metadata, outdir=outdir, verbose=verbose)
        '''
        '''
            Silhouette coefficient: optimal number of clusters
        '''
        print_line(verbose)
        n_clusters = my_silhouette_score(Z, merged_cnvs, outdir, verbose)

    '''
        Extract clusters and compute mean cnv profiles

    '''
    print_line(verbose)
    if n_clusters > 1:
        df_mean = extract_clusters(Z, boundaries, n_clusters, merged_cnvs, reclust, outdir, verbose)
        heatmap(df_mean.transpose(), boundaries, method, metric,outdir, verbose, sample)

    """
    if meta_format == 'xml':
        xml=dicttoxml(metadata,  custom_root='analysis', attr_type=False)
        dom = parseString(xml)
        with open(outdir+'/metadata.xml', 'w+') as f:
            f.write(dom.toprettyxml())
    else:
        with open(outdir+'/metadata.json', 'w+') as f:
            json.dump(metadata, f, indent=2)
    """
