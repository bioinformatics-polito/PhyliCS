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
import pandas as pd 
import numpy as np
from sklearn.metrics import silhouette_score
from sklearn.metrics import pairwise_distances

def het_score_(data:pd.DataFrame, labels:np.array, n_jobs:int=1):
    X = pairwise_distances(data, metric='l1', n_jobs=n_jobs)
    return silhouette_score(X, labels, metric='precomputed')