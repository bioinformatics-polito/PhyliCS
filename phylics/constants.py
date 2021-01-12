#!/usr/bin/env python
from .utils import dotdict

K_METHODS = dotdict({
    'KMEANS' : "kmeans", 
    'AGGLOMERATIVE' :"agglomerative", 
    'BIRCH' : "birch"})

LINKAGE = dotdict({
        'AVERAGE' : "average", 
        'COMPLETE' : "complete", 
        'SINGLE' : "single", 
        'WARD' : "ward"})

EMBEDDINGS = dotdict({
    'PCA' : "pca",
    'UMAP' : "umap"

})