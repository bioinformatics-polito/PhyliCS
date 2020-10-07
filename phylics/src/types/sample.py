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

__all__ = ['Sample']

from .cnvdata import CnvData 
from ..utils import *
from ..tools import variable_features, pca, umap
from ..drawer import Drawer
import random
import pandas as pd
import numpy as np
from typing import Union, Tuple, List, Optional, Mapping, Any, Iterable, Sequence, Dict
from ..utils import AnyRandom, _InitPos
from ..plotting._utils import _Dimensions
from .. import logging as logg


_FILTER_SINGLE_METHODS_ = {  
                                'EQ' :'eq', 
                                'LT' : 'lt',
                                'GT' : 'gt', 
                                'LT_EQ' : 'lt_eq', 
                                'GT_EQ' : 'gt_eq',
                                }
_FILTER_INTERVAL_METHODS_ = {
                                'IN' : 'in', 
                                'OUT' : 'out', 
                                'IN_EQ' : 'in_eq', 
                                'OUT_EQ' : 'out_eq'}

_FILTER_PERCENTILE_METHODS_ = {
                                'LT' : 'lt',
                                'GT' : 'gt', 
                                'LT_EQ' : 'lt_eq', 
                                'GT_EQ' : 'gt_eq',

                                }
_FILTER_LIST_METHODS_ = {
                                'IN' : 'in', 
                                'OUT' : 'out', }


class Sample:

    def __init__(self, 
                cnv_data:CnvData, 
                sample_name:str="sample"):      
        self.cnv_data = cnv_data
        self.name = sample_name

    @classmethod
    def from_file(cls, filepath:str, flavour:str="ginkgo", sample_name:str="sample"):
        if flavour == "ginkgo":
            X, boundaries = from_ginkgo_to_phylics(filepath)
            return cls(cnv_data=CnvData(X, feat=boundaries, obs_names='row_names', feat_names='col_names'), sample_name=sample_name)
        else:
            raise ValueError("Implemented only for ginkgo-like files")
        
    def shape(self):
        return self.cnv_data.shape
        
    def add_annotation(self, ann:Union[pd.Series, Mapping[str, Union[float, int]]], key:str, axis:int=0):
        ann = sanitize_annotation(ann)
        if axis == 0:
            if (ann.index.equals(self.cnv_data.obs.index) == False):
                raise ValueError("indices mismatch")
            self.cnv_data.obs[key] = ann
        elif axis == 1:
            if (ann.index.equals(self.cnv_data.feat.index) == False):
                raise ValueError("indices mismatch")
            self.cnv_data.feat[key] = ann
        else:
            raise ValueError("IllegalArgument: axis = {}. Accepted values are 0 and 1")
    
    def load_annotation(self, filepath:str, key:str, axis:int=0):
        ann = load_annotation_(filepath)
        self.add_annotation(ann, key, axis)
    
    def get_cnv_dataframe(self):
        return self.cnv_data.to_df()
    
    def get_cnv_matrix(self):
        return self.cnv_data.X
    
    def get_cell(self):
        return self.cnv_data.obs_names

    def get_boundaries(self):
        return self.cnv_data.feat[['CHR', 'START', 'END']]

    def get_annotations(self, axis:int=0):
        if axis == 0:
            return self.cnv_data.obs
        elif axis == 1:
            return self.cnv_data.feat
        else:
            raise ValueError("IllegalArgument: axis = {}. Accepted values are 0 and 1")
    

    def get_annotation(self, key:str, axis:int=0):
        if axis == 0:
            return self.cnv_data.obs[key]
        elif axis == 1:
            return self.cnv_data.feat[key]
        else:
            raise ValueError("IllegalArgument: axis = {}. Accepted values are 0 and 1")

    def count(self):
        return self.cnv_data.n_obs
    
    def drop_cells(self, cells:list, inplace:bool=False):
        cnv_data = self.cnv_data.drop_obs(select=cells, inplace=inplace)
        return cnv_data

    def mad(self):
        mad = self.cnv_data.to_df().mad(axis=1)
        self.add_annotation(mad, 'mad')

    def variable_features(self, min_disp: Optional[float] = None, max_disp: Optional[float] = None, min_mean: Optional[float] = None, 
        max_mean: Optional[float] = None, n_top_features: Optional[int] = None, n_bins: int = 20 ):
    
        highly_variable, means, dispersions, dispersions_norm = variable_features(X = self.cnv_data.to_df(), min_disp=min_disp, max_disp=max_disp, min_mean=min_mean,
            max_mean = max_mean, n_top_features = n_top_features, n_bins=n_bins)
        self.cnv_data.feat['highly_variable'] = highly_variable
        self.cnv_data.feat['means'] = means
        self.cnv_data.feat['dispersions'] = dispersions
        self.cnv_data.feat['dispersions_norm'] = dispersions_norm

    def get_variable_features(self):
        return self.cnv_data.feat[['highly_variable', 'means', 'dispersions', 'dispersions_norm']]
    
    def pca(self, n_comps: Optional[int] = None, jackstraw_perms: Optional[int] = None, 
        svd_solver: str = 'arpack', random_state: AnyRandom = 0, 
        use_highly_variable: Optional[bool] = False):
            pca_result = pca(self.cnv_data, n_comps,  jackstraw_perms, svd_solver, random_state, use_highly_variable)
            self.cnv_data.uns['pca'] = {}
            self.cnv_data.uns['pca']['X'] = pca_result[0]
            self.cnv_data.uns['pca']['variance'] = pca_result[1]
            self.cnv_data.uns['pca']['variance_ratio'] = pca_result[2]
            self.cnv_data.uns['pca']['components'] = pca_result[3]
            self.cnv_data.uns['pca']['perm_variance_ratio'] = pca_result[4]

    def get_pca(self, sub_field: Union[str, None]=None):
        if sub_field is None:
            return self.cnv_data.uns['pca']
        else:
            if sub_field == 'X':
                return self.cnv_data.uns['pca']['X']
            elif sub_field == 'variance':
                return self.cnv_data.uns['pca']['variance']
            elif sub_field == 'variance_ratio':
                return self.cnv_data.uns['pca']['variance_ratio']
            elif sub_field == 'components':
                return self.cnv_data.uns['pca']['components']
            elif sub_field == 'perm_variance_ratio':
                return self.cnv_data.uns['pca']['perm_variance_ratio']
            else:
                logg.warning(
                    f"InvalidArgument: sub_field={sub_field} is not valid. "
                )
                return self.cnv_data.uns['pca']
    

    def jack_straw(self):
        return NotImplemented

    def informative_pcs(self, method:str='jackstraw'):
       
        return NotImplemented 

    def umap(self, n_neighbors: int = 15, n_components: int = 2, metric: str = 'euclidean', metric_kwds: Dict = None,
                min_dist: float = 0.5, spread: float = 1.0, maxiter: Optional[int] = None, alpha: float = 1.0, gamma: float = 1.0,
                fast: Union[bool, None] = False, negative_sample_rate: int = 5, local_connectivity: Union[int, None] = 1,
                init_pos: Union[_InitPos, np.ndarray, None] = 'spectral', random_state: AnyRandom = 0, a: Optional[float] = None,
                b: Optional[float] = None):
                umap_result = umap(self.cnv_data, n_neighbors, n_components, metric, metric_kwds, min_dist, spread, maxiter, alpha, gamma,
                                fast, negative_sample_rate, local_connectivity, init_pos, random_state, a, b)
                self.cnv_data.uns['umap'] = {}
                self.cnv_data.uns['umap']['X'] = umap_result
    def get_umap(self, sub_field: Union[str, None]=None):
        if sub_field is None:
            return self.cnv_data.uns['umap']
        else:
            if sub_field == 'X':
                return self.cnv_data.uns['umap']['X']
            else:
                logg.warning(
                    f"InvalidArgument: sub_field={sub_field} is not valid. "
                )
                return self.cnv_data.uns['umap']


    def clusters(self):
        return NotImplemented 


    #Filter functions
    def filter(self, key:str, method:str,
                    value:Union[int, float, Sequence[Union[float, int]], Tuple[Union[int, float], Union[int, float]]], 
                    percentile:bool=False, inplace:bool=False):
        if key not in self.cnv_data.obs.columns:
            raise ValueError('Did not find cnv_data.obs[\'{}\']. '.format(key))
        if isinstance(value, tuple):
            return self._filter_interval(key, method, value, inplace=inplace)
        elif isinstance(value, Union[int, float]):
            return self._filter_value(key, method, value, percentile=percentile, inplace=inplace)
        elif isinstance(value, cabc.Sequence):
            return self._filter_list(key, method, value, inplace=inplace)
        else:
            raise TypeError("value: expected int, string, sequence or tuple but got {}".format(type(value)))

    def _filter_interval(self, key:str, method:str, value: Tuple[Union[int, float], Union[int, float]], inplace:bool=False):
        if method not in list(_FILTER_INTERVAL_METHODS_.values()):
            raise ValueError("'method' must be one of %r"% list(_FILTER_INTERVAL_METHODS_.values()))
        if value[0] >= value[1]:
            raise ValueError("Wrong interval of values: value_0 must be smaller than value_1")
        annotation = self.get_annotation(key)
        if method == _FILTER_INTERVAL_METHODS_['IN']:
            subset = self.cnv_data[annotation > value[0] and annotation < value[1], :]
        elif method ==  _FILTER_INTERVAL_METHODS_['IN_EQ']:
            subset = self.cnv_data[annotation >= value[0] and annotation <= value[1], :]
        elif method ==  _FILTER_INTERVAL_METHODS_['OUT']:
            subset = self.cnv_data[annotation < value[0] or annotation > value[1], :]
        elif method ==  _FILTER_INTERVAL_METHODS_['OUT_EQ']:
            subset = self.cnv_data[annotation < value[0] or annotation >= value[1], :]
        
        if inplace == True:
            self.cnv_data = subset
        else:
            return Sample(subset, self.name)

    def _filter_value(self, key:str, method:str, value:Union[int, float], percentile:bool=False, inplace:bool=False):
        annotation = self.get_annotation(key)
        if percentile == False: 
            if method not in list(_FILTER_SINGLE_METHODS_.values()):
                raise ValueError("'method' must be one of %r"% list(_FILTER_SINGLE_METHODS_.values()))
            if method == _FILTER_SINGLE_METHODS_["EQ"]:
                subset = self.cnv_data[annotation == value, :]
            elif method == _FILTER_SINGLE_METHODS_["LT"]:
                subset = self.cnv_data[annotation < value, :]
            elif method == _FILTER_SINGLE_METHODS_["LT_EQ"]:
                subset = self.cnv_data[annotation <= value, :]
            elif method == _FILTER_SINGLE_METHODS_["GT"]:
                subset = self.cnv_data[annotation > value, :]
            elif method == _FILTER_SINGLE_METHODS_["GT_EQ"]:
                subset = self.cnv_data[annotation >= value, :]
        else:
            if method not in list(_FILTER_PERCENTILE_METHODS_.values()):
                raise ValueError("'method' must be one of %r"% list(_FILTER_PERCENTILE_METHODS_.values()))
            threshold = annotation.quantile(value)
            if method == _FILTER_SINGLE_METHODS_["LT"]:
                subset = self.cnv_data[annotation < threshold, :]
            elif method == _FILTER_SINGLE_METHODS_["LT_EQ"]:
                subset = self.cnv_data[annotation <= threshold, :]
            elif method == _FILTER_SINGLE_METHODS_["GT"]:
                subset = self.cnv_data[annotation > threshold, :]
            elif method == _FILTER_SINGLE_METHODS_["GT_EQ"]:
                subset = self.cnv_data[annotation >= threshold, :]
        if inplace == True:
            self.cnv_data = subset
        else:
            return Sample(subset, self.name)

    def _filter_list(self, key:str, method:str, value: Sequence[Union[int, float]], inplace:bool=False):
        if method not in list(_FILTER_LIST_METHODS_.values()):
            raise ValueError("'method' must be one of %r"% list(_FILTER_LIST_METHODS_.values()))
        annotation = self.get_annotation(key)
        if method == _FILTER_LIST_METHODS_['IN']:
            subset = self.cnv_data[annotation.isin(value), :]
        elif method ==  _FILTER_LIST_METHODS_['OUT']:
            subset = self.cnv_data[~annotation.isin(value), :]

        if inplace == True:
            self.cnv_data = subset
        else:
            return Sample(subset, self.name)    

    def plot_annotation_dist(self, ann:str, axis:int=0, grid:bool=False, quantiles:List[float]=None, figsize:Tuple[int, int]=None, 
                                outpath:str=None, **kwargs):
        Drawer.draw('dist', self.get_annotation(key=ann, axis=axis), grid, quantiles, figsize, outpath, **kwargs)

    def heatmap(self, method:str ='ward', metric:str ='euclidean', outpath:str=None,
                    vmin:int = 0, vmax:int = 12, vcenter:int=2, figsize:Tuple[int, int]=(37, 21), fontsize:int=16):
        Drawer.draw('heatmap', self.get_cnv_matrix(), self.get_boundaries(), method, metric,  outpath=outpath, sample=self.name, 
            vmin = vmin, vmax=vmax, vcenter=vcenter, figsize=figsize, fontsize=fontsize)

    """
    def scatter_plot(self, outpath:str=None, figsize:Tuple[int, int]=None, **kwargs):
        embeddings = Reducer.umap_(self.cnv_data.get_cnvs())
        Drawer.draw('scatter', X=embeddings, figsize=figsize, outpath=outpath, **kwargs)
    """
    
    def plot_variable_features(self, log: bool = False, outpath: str = None):
        if hasattr(self.cnv_data.feat, 'highly_variable'):
            Drawer.draw('variable_features', X=self.get_variable_features(), log=log, outpath=outpath)
        else:
            raise AttributeError("{} object has no attribute 'highly_variable'".format(self.cnv_data.feat))
    
    def plot_pca(self, projection: _Dimensions = "2d", outpath:str=None, 
                figsize:Tuple[int, int]=None, **kwargs):
        if 'pca' in self.cnv_data.uns:
            if projection == '3d':
                z_label = 'PC3'
            else:
                z_label = None
            Drawer.draw('scatter', data=self.get_pca('X'), title = 'PCA projection', x_label='PC1', y_label='PC2', z_label=z_label,
                projection = projection, outpath=outpath,  
                figsize = figsize, **kwargs)
        else:
            logg.error("{} object has no attribute 'pca'".format(self.cnv_data.uns))

    
    def plot_explained_variance(self, n_pcs:int = 50, figsize:Tuple[int, int]=None,
                show: bool=False, outpath: str = None):
        if 'pca' in self.cnv_data.uns:
            Drawer.draw('jackstraw', explained_variance_ratio= self.get_pca('variance_ratio'), 
                        mean_expl_var_ratio_perm=self.get_pca('perm_variance_ratio'), n_pcs = n_pcs, figsize=figsize,
                        show = show, outpath = outpath)
        else:
            logg.error("{} object has no attribute 'pca'".format(self.cnv_data.uns))
        
    def plot_umap(self, projection: _Dimensions = "2d", outpath:str=None, 
                figsize:Tuple[int, int]=None, **kwargs):
        if 'umap' in self.cnv_data.uns:
            if projection == '3d':
                z_label = 'Z'
            else:
                z_label = None
            Drawer.draw('scatter', data=self.get_umap('X'), title = 'UMAP embeddings', x_label='X', y_label='Y', z_label=z_label,
                projection = projection, outpath=outpath,  
                figsize = figsize, **kwargs)
        else:
            logg.error("{} object has no attribute 'umap'".format(self.cnv_data.uns))

    
    



    
    
    


