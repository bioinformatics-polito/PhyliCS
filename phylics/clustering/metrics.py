#!/usr/bin/env python
import os
import sys
import numpy as np
import pandas as pd 
from typing import Union, Optional 
from ..types import CnvData
from joblib import Parallel, delayed
from ..utils import clustering_func
from ..constants import K_METHODS
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score, adjusted_mutual_info_score, adjusted_rand_score, v_measure_score, fowlkes_mallows_score
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


def __index__(k, index, data, method, metric, linkage):
    if method == K_METHODS.AGGLOMERATIVE:
        clustering = clustering_func[method](n_clusters=k, affinity=metric, linkage=linkage).fit(data)
    else:
        clustering = clustering_func[method](n_clusters=k).fit(data)
    labels_ = clustering.labels_
    if index == silhouette_score:
        scores = index(data, labels_, metric=metric)
    else: 
        scores = index(data, labels_)
    return scores

def silhouette_(data: np.ndarray, method:str, metric:Optional[Union[str, None]]="euclidean", 
            linkage:Optional[Union[str, None]]=None, 
            min_k:Optional[int]=2, max_k:Optional[int]=15, n_jobs:Optional[int]=1):
    scores_dict = {}
    scores = Parallel(n_jobs=n_jobs)(delayed(__index__)(k, silhouette_score, data, method, metric, linkage) 
                                                                for k in np.arange(min_k, max_k+1))
    for k, score in zip (np.arange(min_k, max_k+1), scores):
        scores_dict[k] = score
    return scores_dict

def davies_bouldin_(data: np.ndarray, method:str, metric:Optional[Union[str, None]]="euclidean", 
            linkage:Optional[Union[str, None]]=None, 
            min_k:Optional[int]=2, max_k:Optional[int]=15, n_jobs:Optional[int]=1):
    scores_dict = {}
    scores = Parallel(n_jobs=n_jobs)(delayed(__index__)(k, davies_bouldin_score, data, method, metric, linkage) 
                                                                for k in np.arange(min_k, max_k+1))
    for k, score in zip (np.arange(min_k, max_k+1), scores):
        scores_dict[k] = score
    return scores_dict

def calinski_harabasz_(data: np.ndarray, method:str, metric:Optional[Union[str, None]]="euclidean", 
            linkage:Optional[Union[str, None]]=None, 
            min_k:Optional[int]=2, max_k:Optional[int]=15, n_jobs:Optional[int]=1):
    scores_dict = {}
    scores = Parallel(n_jobs=n_jobs)(delayed(__index__)(k, calinski_harabasz_score, data, method, metric, linkage) 
                                                                for k in np.arange(min_k, max_k+1))
    for k, score in zip (np.arange(min_k, max_k+1), scores):
        scores_dict[k] = score
    return scores_dict

def cluster_accuracy_(clusterer_:object, data:np.ndarray, true_labels:Union[list, np.array]):
    pred_labels = clusterer_.fit_predict(data)
    ari = adjusted_rand_score(true_labels, pred_labels)
    ami = adjusted_mutual_info_score(true_labels, pred_labels)
    vm = v_measure_score(true_labels, pred_labels)
    fm = fowlkes_mallows_score(true_labels, pred_labels)
    return pd.DataFrame({'method':[clusterer_.method], 'ari':[ari], 'ami':[ami], 'vm':[vm], 'fm':[fm]})