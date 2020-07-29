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
import numpy as np
import pandas as pd
from typing import Union, List, Optional, Sequence
from  sklearn.decomposition import PCA
import umap
from .utils import AnyRandom


def variable_features(X:CnvData, min_disp: Optional[float] = None, max_disp: Optional[float] = None,
    min_mean: Optional[float] = None, max_mean: Optional[float] = None, n_top_features: Optional[int] = None,
    n_bins: int = 20 ):
    df = _highly_variable_features(X = X, min_disp=min_disp, max_disp=max_disp, min_mean=min_mean,
            max_mean = max_mean, n_top_features = n_top_features, n_bins=n_bins)

    return df['highly_variable'], df['means'], df['dispersions'], df['dispersions_norm']
       

def pca(data: Union[CnvData, np.ndarray], n_comps: Optional[int] = None, svd_solver: str = 'arpack', random_state: AnyRandom = 0, 
            use_highly_variable: Optional[bool] = False):
    return _pca(data, n_comps, svd_solver, random_state, use_highly_variable)


    
