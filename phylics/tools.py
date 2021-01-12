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
from .clustering.clustering_stats import silhouette_, davies_bouldin_, calinski_harabasz_
import numpy as np
import pandas as pd
from typing import Union, List, Optional, Sequence, Dict
from  sklearn.decomposition import PCA
from .utils import AnyRandom, _InitPos
from . import logging as logg

K_METHODS = ["kmeans", "agglomerative", "birch"]
LINKAGE = ["average", "complete", "single", "ward"]
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
def cluster(data: Union[CnvData, np.ndarray],method:Optional[str]="hdbscan", metric:Optional[str]="l1", 
            embeddings:Optional[Union[str, None]]=None, **kwargs):
            return NotImplemented

def clust_stats(data: Union[CnvData, np.ndarray], method:str, metric:Optional[Union[str, None]]=None, 
            linkage:Optional[Union[str, None]]=None, embeddings:Optional[Union[str, None]]=None, 
            min_k:Optional[int]=2, max_k:Optional[int]=15, index:Optional[Union[str, Sequence[str]]]="all"):
        """
        Calculate clustering statistics to propose to the user the best k for kmeans, agglomerative clustering and birch,
        according to three internal validation indices (silhoette, calinski-harabasz, davies-bouldin).

        It allows also to select the distance metric and the linkage method, for agglomerative clustering.
        """
        if method not in K_METHODS:
                raise ValueError("method must be one of [" + ', '.join(K_METHODS) +"]")
        if linkage is not None and not in LINKAGE:
                raise ValueError("linkage must be one of [" + ', '.join(LINKAGE) +"]")
        if metric 

        if method == "agglomerative" and linkage != "l2":
                logg.error("agglomerative clustering with ward linkage can accept only l2 distance")
                metric = "l2"
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
        elif isinstance(index, Sequence[str]):
                if "silhouette" in index:
                        silhouette = True 
                if "db" in index:
                        db = True 
                if "ch" in index:
                        ch = True
        indices = {"silhouette":{}, "db":{}, "ch":{}}
        if silhouette == True:
                indices["silhouette"] = silhouette_(data, method, metric, linkage, min_k, max_k)
        if db == True:
                indices["db"] = davies_bouldin_(data, method, metric, linkage, min_k, max_k)
        if ch == True:
                indices["ch"] = calinski_harabasz_(data, method, metric, linkage, min_k, max_k)

        return indices
        
    
