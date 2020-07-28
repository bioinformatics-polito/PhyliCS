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
from .data_types import CnvData, VariableFeatures
from .preprocessing._highly_variant_features import _highly_variable_features
import numpy as np
import pandas as pd
from typing import Union, List, Optional, Sequence
from  sklearn.decomposition import PCA
import umap


def variable_features(X:CnvData, min_disp: Optional[float] = None, max_disp: Optional[float] = None,
    min_mean: Optional[float] = None, max_mean: Optional[float] = None, n_top_features: Optional[int] = None,
    n_bins: int = 20 ):
    df = _highly_variable_features(X = X, min_disp=min_disp, max_disp=max_disp, min_mean=min_mean,
            max_mean = max_mean, n_top_features = n_top_features, n_bins=n_bins)

    return df['highly_variable'], df['means'], df['dispersions'], df['dispersions_norm']
       
def informative_pcs(X):
    return NotImplemented

class Reducer:
    @staticmethod
    def umap_(X, features:Union[Sequence, pd.Series, np.array] = None, **kwargs):
        if features != None:
            X = X[[features]]
        reducer = umap.UMAP(**kwargs)
        return pd.DataFrame(reducer.fit_transform(X), index=X.index)
    
    @staticmethod
    def pca_(X, features:Union[Sequence, pd.Series, np.array] = None, **kwargs):
        if features != None:
            X = X[[features]]
        pca = PCA(**kwargs)
        return pd.DataFrame(pca.fit_transform(X), index=X.index)
    

    


#class Clusterer: