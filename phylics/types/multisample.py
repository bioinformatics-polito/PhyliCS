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
# sample.py: this module implements all methods to manage samples
# ==========================================================================

__all__ = ['MultiSample']


from .sample import Sample
from .cnvdata import CnvData 
from ..tools import segregation_score
import pandas as pd
from ..drawer import Drawer
import numpy as np
from typing import Union, Tuple, List, Optional, Mapping, Any, Iterable, Sequence, Dict
from .. import logging as logg
from ..tools import umap

class MultiSample(Sample):
    def __init__(self, 
                cnv_data:CnvData, 
                sample_names:str="MultiSample"):      
        super().__init__(cnv_data, sample_names)

    @classmethod
    def from_list(cls, samples:list):
        sample_names = ""
        boundaries = None
        X = None
        sample_df = pd.DataFrame(columns=["sample"])
        for s in samples:
            sample_names += s.name
            if boundaries is None:
                boundaries = s.get_boundaries()
            if X is None:
                X = s.get_cnv_dataframe()
            else:
                X = X.append(s.get_cnv_dataframe(), ignore_index=True)
            sample_list =  [s.name] * s.shape()[0]
            sample_df = sample_df.append(pd.DataFrame({"sample":sample_list}), ignore_index=True)
        return cls(cnv_data=CnvData(X, feat=boundaries, obs=sample_df, obs_names='row_names', feat_names='col_names'), sample_names=sample_names)
    
    def get_sample_labels(self):
        return self.get_annotation("sample")

    def segregation_score(self, n_jobs:int=1, verbose:bool=False):
        return segregation_score(self.cnv_data, n_jobs, verbose)
    
    def plot_samples(self, outpath:str=None, 
                figsize:Tuple[int, int]=None, **kwargs):
        if 'umap' in self.cnv_data.uns:
            projection = self.get_umap('X')
        else:
            projection = umap(self.cnv_data)
        Drawer.draw('scatter', data=projection, title = 'Multi sample projection', x_label='X', y_label='Y', outpath=outpath,  
                    figsize = figsize, labels=self.cnv_data.obs['sample'].values, legend=True, **kwargs)


    def plot_dendrogram(self, outpath:str=None, clusters:bool=False, figsize:Tuple[int, int]=None, **kwargs):
        #TODO manage multi-layer labels
        if clusters == True: 
            if 'label' in self.cnv_data.obs.columns:
                model =  self.get_clusterer()
                if hasattr(model, 'children_'):
                    #agglomerative clustering has been computed
                    counts = np.zeros(model.children_.shape[0])
                    n_samples = len(model.labels_)
                    for i, merge in enumerate(model.children_):
                        current_count = 0
                        for child_idx in merge:
                            if child_idx < n_samples:
                                current_count += 1  # leaf node
                            else:
                                current_count += counts[child_idx - n_samples]
                        counts[i] = current_count

                    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)
                    Drawer.draw('heatmap', data=self.get_cnv_dataframe(), boundaries=self.get_boundaries(), linkage=linkage_matrix,  title = 'Multi dataset cluster dendrogram & heatmap', outpath=outpath,  
                        labels=self.cnv_data.obs[['sample', 'label']], legend=True, **kwargs)
                else:
                    logg.error("Cluster model has no attribute 'children'. Consider executing agglomerative clustering.")
            else:
                logg.error("{} object has no column 'label'".format(self.cnv_data.feat))
        else:
            Drawer.draw('heatmap', data=self.get_cnv_dataframe(), boundaries=self.get_boundaries(), row_cluster=True,  title = 'Multi dataset dendrogram & heatmap',  outpath=outpath,  
                        labels=self.cnv_data.obs['sample'].values, legend=True, **kwargs)
