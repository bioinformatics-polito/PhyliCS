#!/usr/bin/env python

# ==========================================================================
#                                  PhyliCS
# ==========================================================================
# This file is part of PhyliCS.
#
# PhyliCS is Free Software: you can redistribute it and/or modify it
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
# custom_types.py: this module implements phylics worker types
# ==========================================================================
from .preprocessing._highly_variant_features import _highly_variable_features
from .types import CnvData
from .preprocessing._pca import _pca
from .preprocessing._umap import _umap
from .clustering.metrics import silhouette_, davies_bouldin_, calinski_harabasz_, cluster_accuracy_
from .clustering.utils import ClusterConfig
from .clustering.cluster import cluster_
import numpy as np
import pandas as pd
from typing import Union, List, Optional, Sequence, Dict
from  sklearn.decomposition import PCA
from .utils import AnyRandom, _InitPos, clustering_func, partition
from .constants import K_METHODS, LINKAGE, EMBEDDINGS
from . import logging as logg
from .multisample.spatiality import het_score_
from joblib import Parallel, delayed

METRICS = ["euclidean", "l1", "l2", "manhattan", "cosine"]

def variable_features(X:CnvData, min_disp: Optional[float] = None, max_disp: Optional[float] = None,
    min_mean: Optional[float] = None, max_mean: Optional[float] = None, n_top_features: Optional[int] = None,
    n_bins: int = 20 ):
    df = _highly_variable_features(X = X, min_disp=min_disp, max_disp=max_disp, min_mean=min_mean,
            max_mean = max_mean, n_top_features = n_top_features, n_bins=n_bins)

    return df['highly_variable'], df['means'], df['dispersions'], df['dispersions_norm']
       

def pca(data: Union[CnvData, np.ndarray], n_comps: Optional[int] = None, jackstraw_perms: Optional[int] = None, svd_solver: str = 'arpack', random_state: AnyRandom = 0, 
            use_highly_variable: Optional[bool] = False):
    return _pca(data, n_comps, jackstraw_perms, svd_solver, random_state, use_highly_variable)

def umap(data: Union[CnvData, np.ndarray], n_neighbors: int = 15, n_components: int = 2, metric: str = 'euclidean', metric_kwds: Dict = None,
                min_dist: float = 0.5, spread: float = 1.0, maxiter: Optional[int] = None, alpha: float = 1.0, gamma: float = 1.0,
                fast: Union[bool, None] = False, negative_sample_rate: int = 5, local_connectivity: Union[int, None] = 1,
                init_pos: Union[_InitPos, np.ndarray, None] = 'spectral', random_state: AnyRandom = 0, a: Optional[float] = None,
                b: Optional[float] = None, use_highly_variable: Optional[bool] = False):
        return _umap(data, n_neighbors, n_components, metric, metric_kwds, min_dist, spread, maxiter, alpha, gamma,
                                fast, negative_sample_rate, local_connectivity, init_pos, random_state, a, b, use_highly_variable)

def nk_clust(data: Union[CnvData, np.ndarray], method:str, metric:Optional[Union[str, None]]="euclidean", 
            linkage:Optional[Union[str, None]]=None, embeddings:Optional[Union[str, None]]=None, 
            min_k:Optional[int]=2, max_k:Optional[int]=15, index:Optional[Union[str, List[str]]]="all",
            n_jobs:Optional[int]=1):
        """
        Calculate clustering statistics to propose to the user the best k for kmeans, agglomerative clustering and birch,
        according to three internal validation indices (silhoette, calinski-harabasz, davies-bouldin).

        It allows also to select the distance metric and the linkage method, for agglomerative clustering.
        """
        if method not in K_METHODS.values():
                raise ValueError("Accepted values for method are " + ', '.join(K_METHODS.values()) +"")
        if (linkage is not None) and (linkage not in LINKAGE.values()):
                raise ValueError("Accepted values for linkage are " + ', '.join(LINKAGE.values()) +"")
        if (metric is not None) and (metric not in METRICS):
                raise ValueError("Accepted values for metric are " + ', '.join(METRICS) +"")
        """
        if embeddings is not None:
            if embeddings == EMBEDDINGS.PCA:
                if "pca" not in cnvdata.uns.keys():
                    raise ValueError(
                        'Did not find cnv.uns[\'pca\']. '
                        'Consider running `Sample.pca()` first.'
                    )
                else:
                    data_comp = (cnvdata[:, cnvdata.uns['pca']["X"] ])
            elif embeddings == EMBEDDINGS.UMAP:
                if "umap" not in cnvdata.uns.keys():
                    raise ValueError(
                        'Did not find cnv.uns[\'umap\']. '
                        'Consider running `Sample.umap()` first.'
                    )
                else:
                    data_comp = (cnvdata[:, cnvdata.uns['umap']["X"] ])
            else:
                raise ValueError("Accepted values for embeggings are 'umap', 'pca'")
        else: 
            data_comp = (cnvdata.X)
        """
        if min_k < 2:
                raise ValueError("min_k must be at least equal to 2")

        if max_k < min_k:
                raise ValueError("max_k must be greater or equal to min_k")

        if method == K_METHODS.AGGLOMERATIVE and linkage == LINKAGE.WARD and metric != "euclidean":
                logg.warning("agglomerative clustering with ward linkage can accept only euclidean distance.\nFixed by default.")
                metric = "euclidean"
        silhouette = False
        db = False
        ch = False
        if isinstance(index, str):
                if index == "silhouette" or index == "all":
                        silhouette = True 
                if index == "db" or index == "all":
                        db = True
                if index == "ch" or index == "all":
                        ch = True 
        elif isinstance(index, List):
                if "silhouette" in index:
                        silhouette = True 
                if "db" in index:
                        db = True 
                if "ch" in index:
                        ch = True
        indices = {}
        if silhouette == True:
                indices["silhouette"] = silhouette_(data, method, metric, linkage, min_k, max_k, n_jobs)
        if db == True:
                indices["db"] = davies_bouldin_(data, method, metric, linkage,  min_k, max_k, n_jobs)
        if ch == True:
                indices["ch"] = calinski_harabasz_(data, method, metric, linkage, min_k, max_k, n_jobs)

        return indices
"""
def consensus_clustering(data: Union[CnvData, np.ndarray], method:str, cluster_configs:List[ClusterConfig],
                n_iter:Optional[int]=5, embeddings:Optional[Union[str, None]]=None, n_jobs:Optional[int]=1):
        cnvdata = data if isinstance(data, CnvData) else CnvData(data)

        if method not in clustering_func.keys():
                raise ValueError("Accepted values for method are [" + ', '.join(clustering_func.keys()) +"]")
        
        if embeddings is not None:
            if embeddings == EMBEDDINGS.PCA:
                if "pca" not in cnvdata.uns.keys():
                    raise ValueError(
                        'Did not find cnv.uns[\'pca\']. '
                        'Consider running `Sample.pca()` first.'
                    )
                else:
                    data_comp = (cnvdata[:, cnvdata.uns['pca']["X"] ])
            elif embeddings == EMBEDDINGS.UMAP:
                if "umap" not in cnvdata.uns.keys():
                    raise ValueError(
                        'Did not find cnv.uns[\'umap\']. '
                        'Consider running `Sample.umap()` first.'
                    )
                else:
                    data_comp = (cnvdata[:, cnvdata.uns['umap']["X"] ])
            else:
                raise ValueError("Accepted values for embeggings are ['umap', 'pca']")
        else: 
            data_comp = (cnvdata.X)

        return ConsensusCluster(clustering_func, configurations, n_iter).fit(data_comp)
"""
def cluster_accuracy(clusterer_, data, labels):
        return cluster_accuracy_(clusterer_, data, labels)

def clust_accuracy_ranking(data:Union[CnvData, np.ndarray], cluster_configurations:ClusterConfig, 
                        labels:Union[list, np.array], embeddings:Optional[Union[str, None]]=None, n_jobs:Optional[int]=1):
        cnvdata = data if isinstance(data, CnvData) else CnvData(data)

        if embeddings is not None:
            if embeddings == EMBEDDINGS.PCA:
                if "pca" not in cnvdata.uns.keys():
                    raise ValueError(
                        'Did not find cnv.uns[\'pca\']. '
                        'Consider running `Sample.pca()` first.'
                    )
                else:
                    data_comp = (cnvdata[:, cnvdata.uns['pca']["X"] ])
            elif embeddings == EMBEDDINGS.UMAP:
                if "umap" not in cnvdata.uns.keys():
                    raise ValueError(
                        'Did not find cnv.uns[\'umap\']. '
                        'Consider running `Sample.umap()` first.'
                    )
                else:
                    data_comp = (cnvdata[:, cnvdata.uns['umap']["X"] ])
            else:
                raise ValueError("Accepted values for embeggings are 'umap', 'pca'")
        else: 
            data_comp = (cnvdata.X)
        results = Parallel(n_jobs=n_jobs)(delayed(cluster_accuracy)(clusterer_, data_comp, labels) for clusterer_ in cluster_configurations.clusterer_)
        ranking = pd.concat(results, ignore_index=True)
        #rank results according to the average ranking of the four indices
        ranking['rank'] = ranking[['ari', 'ami', 'fmi', 'vm']].rank().mean(axis=1)
        ranking = ranking.sort_values(by='rank', ascending=False)
        ranking = ranking.drop('rank', axis=1)
        ranking = ranking.set_index('method')
        return ranking

def cluster(data:Union[CnvData, np.ndarray], method:str, embeddings:Optional[Union[str, None]]=None, **kwargs):
        cnvdata = data if isinstance(data, CnvData) else CnvData(data)
        if embeddings is not None:
            if embeddings == EMBEDDINGS.PCA:
                if "pca" not in cnvdata.uns.keys():
                    raise ValueError(
                        'Did not find cnv.uns[\'pca\']. '
                        'Consider running `Sample.pca()` first.'
                    )
                else:
                    data_comp = (cnvdata[:, cnvdata.uns['pca']["X"] ])
            elif embeddings == EMBEDDINGS.UMAP:
                if "umap" not in cnvdata.uns.keys():
                    raise ValueError(
                        'Did not find cnv.uns[\'umap\']. '
                        'Consider running `Sample.umap()` first.'
                    )
                else:
                    data_comp = (cnvdata[:, cnvdata.uns['umap']["X"] ])
            else:
                raise ValueError("Accepted values for embeggings are 'umap', 'pca'")
        else: 
            data_comp = (cnvdata.X)
        
        return cluster_(data_comp, method, **kwargs)

def segregation_score(data:CnvData, n_jobs:int=1, verbose:bool=False):
    labels = data.obs['sample'].values
    samples = np.unique(labels)
    partitions = partition(samples)
    #tries all possibile partitions of sample list and compute the score
    #in other words, samples are combined in subsets and the most convenient
    #scheme is found
    scores = pd.DataFrame(columns=['samples_partition', 'score'])
    if verbose:
        print("Computing spatial segretation score")
    for p in partitions:
        if len(p) == 1:
            continue
        if len(p) < len(samples):
            #ex. p = [[1], [2], [3, 4]]
            # subsets = [1], [2], [3, 4]
            # cells from sample 3 and 4 will be assigned to the same sample
            new_labels = np.array(labels)
            #subset_dict = dict(zip(range(len(p)), p))   
            for idx, subset in enumerate(p):
                for l in subset:
                    #for each sample in this subset, aggregate
                    #find cells from sample l
                    # get positions in labels which value is l
                    new_labels[labels == l] = idx
            s = het_score_(data.X, new_labels, n_jobs)
            scores = scores.append(pd.DataFrame({'samples_partition':[p], 'score':[s]}), ignore_index=True)
        else:
            #compute overall score
            s = het_score_(data.X, labels, n_jobs)
            scores = scores.append(pd.DataFrame({'samples_partition':[p], 'score':[s]}), ignore_index=True)
        if verbose:
            print("Samples' partition: " + " vs ".join(str(pp) for pp in p) + " - score: " + str(s))
    return scores