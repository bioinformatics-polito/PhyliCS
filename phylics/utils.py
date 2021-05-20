import pandas as pd
from numpy import random
from ._compat import Literal
from typing import Union, Mapping, Sequence
import collections.abc as cabc 
from sklearn.cluster import KMeans, AgglomerativeClustering, Birch, AffinityPropagation, DBSCAN, OPTICS, SpectralClustering

__all__ = ["from_ginkgo_to_phylics", 
            "load_annotation_",
            "sanitize_annotation"]

AnyRandom = Union[None, int, random.Generator] 

_InitPos = Literal['spectral', 'random']

clustering_func = {
    "kmeans" : KMeans,
    "agglomerative": AgglomerativeClustering,
    "birch": Birch,
    "affinity": AffinityPropagation, 
    "dbscan": DBSCAN,
    #"hdbscan": HDBSCAN,
    "optics": OPTICS,
    "spectral":SpectralClustering
}

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def from_ginkgo_to_phylics(filepath:str): 
    df = pd.read_csv(filepath, sep="\t")
    X = df.drop(['CHR', 'START', 'END'], axis=1).transpose()
    boundaries = df[['CHR', 'START', 'END']]
    return X, boundaries

def load_annotation_(filepath:str):
    return pd.read_csv(filepath, header=None, squeeze=True, index_col=0, sep="\t")

def sanitize_annotation(ann:Union[pd.Series, Mapping[str, Union[float, int]]]):
    if isinstance(ann, pd.Series):
        return ann
    elif isinstance(ann, cabc.Mapping):
        return pd.Series(ann)
    else:
        raise ValueError("The annotation must be of type pd.Series or Mapping: got {})".format(type(ann)))

def partition(set_):
    if len(set_) == 0:
        yield []
        return
    for i in range(int(2**len(set_)/2)):
        parts = [[], []]
        for item in set_:
            parts[i&1].append(item)
            i >>= 1
        for b in partition(parts[1]):
            yield [parts[0]]+b






