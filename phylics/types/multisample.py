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
        super().__init__(cnv_data, "MultiSample")
        self.sample_names = sample_names

    @classmethod
    #TODO fix arguments: undefined list of arguments
    def from_list(cls, *samples:tuple):
        sample_names = ""
        boundaries = None
        X = None
        sample_df = pd.DataFrame(columns=["sample"])
        for s in samples:
            if sample_names == "":
                sample_names += s.name
            else:
                sample_names += "_" + s.name
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

    def SHscores(self, n_jobs:int=1, verbose:bool=False):
        scores = segregation_score(self.cnv_data, n_jobs, verbose)
        self.cnv_data.uns["scores"] = scores
        return scores

    def get_SHscores(self):
        return self.cnv_data.uns["scores"]

    def plot_SHscores(self, outpath:str=None, xticks_orientation:int=0, figsize:tuple=(18, 7)):
        if 'scores' not in self.cnv_data.uns.keys():
            raise ValueError(
                        'Did not find cnv.uns[\'scores\']. '
                        'Consider running `MultiSampleSample.SHscore()` first.'
                    )
        scores = self.cnv_data.uns["scores"]
        partitions = scores["samples_partition"].values
        data =  scores["score"].values
        yticks = []
        for pp in partitions:
            str = ''
            for p in pp:
                if str == "":
                    str = '+'.join(p)
                else:
                    str += " vs " + "+".join(p)
            yticks.append(str)
        Drawer.draw('dots', data=data, yticks=yticks, title="SHscore", 
                    x_label="sample partition", y_label="score", figsize=figsize, outpath=outpath)
    
    def plot_samples(self, outpath:str=None, 
                figsize:Tuple[int, int]=None, **kwargs):
        if 'umap' in self.cnv_data.uns.keys():
            projection = self.get_umap('X')
        else:
            projection = umap(self.cnv_data)
        Drawer.draw('scatter', data=projection, title = 'Multi sample projection', x_label='X', y_label='Y', outpath=outpath,  
                    figsize = figsize, labels=self.cnv_data.obs['sample'].values, legend=True, **kwargs)


    def plot_dendrogram(self, outpath:str=None, clusters:bool=False, figsize:Tuple[int, int]=None, **kwargs):
        if clusters == True: 
            if 'cluster' in self.cnv_data.obs.columns:
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
                    Drawer.draw('heatmap', data=self.get_cnv_dataframe(), row_cluster=True, boundaries=self.get_boundaries(), linkage=linkage_matrix,  title = 'Cluster dendrogram & heatmap', outpath=outpath,  
                        labels=self.cnv_data.obs[['sample','cluster']], legend=False, **kwargs)
                else:
                    logg.error("Cluster model has no attribute 'children'. Consider executing agglomerative clustering.")
            else:
                logg.error("{} object has no column 'cluster'".format(self.cnv_data.feat))
        else:
            Drawer.draw('heatmap', data=self.get_cnv_dataframe(), row_cluster=True, boundaries=self.get_boundaries(),  title = 'Dataset dendrogram & heatmap', outpath=outpath, labels=self.cnv_data.obs['sample'], 
                         legend=True, **kwargs)
    


