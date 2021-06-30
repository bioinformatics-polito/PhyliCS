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
#from .multisample import MultiSample
from .. import utils
from ..tools import variable_features, pca, umap, nk_clust, clust_accuracy_ranking, cluster, zscore
from ..drawer import Drawer
import random
import pandas as pd
import numpy as np
from typing import Union, Tuple, List, Optional, Mapping, Any, Iterable, Sequence, Dict
from ..utils import AnyRandom, _InitPos
from ..plotting._utils import _Dimensions
from .. import logging as logg
import collections.abc as cabc
from ..clustering.utils import ClusterConfig
from ..constants import EMBEDDINGS


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
            X, boundaries = utils.from_ginkgo_to_phylics(filepath)
            return cls(cnv_data=CnvData(X, feat=boundaries, obs_names='row_names', feat_names='col_names'), sample_name=sample_name)
        else:
            raise ValueError("Implemented only for ginkgo-like files")
        
    def shape(self):
        return self.cnv_data.shape

    def copy(self) -> "Sample":
            return copy.deepcopy(self)
        
    def add_annotation(self, ann:Union[pd.Series, Mapping[str, Union[float, int]]], key:str, axis:str="obs"):
        ann = utils.sanitize_annotation(ann)
        if axis == "obs":
            if (ann.index.equals(self.cnv_data.obs.index) == False):
                raise ValueError("indices mismatch")
            self.cnv_data.obs[key] = ann
        elif axis == "feat":
            if (ann.index.equals(self.cnv_data.feat.index) == False):
                raise ValueError("indices mismatch")
            self.cnv_data.feat[key] = ann
        else:
            raise ValueError("IllegalArgument: axis = {}. Accepted values are obs and feat")
    
    def load_annotation(self, filepath:str, key:str, axis:str="obs"):
        ann = utils.load_annotation_(filepath)
        self.add_annotation(ann, key, axis)
    
    def get_cnv_dataframe(self):
        return self.cnv_data.to_df()
    
    def get_cnv_matrix(self):
        return self.cnv_data.X
    
    def get_cell(self):
        return self.cnv_data.obs_names

    def get_boundaries(self):
        return self.cnv_data.feat[['CHR', 'START', 'END']]

    def get_annotations(self, axis:str="obs"):
        if axis == "obs":
            return self.cnv_data.obs
        elif axis == "feat":
            return self.cnv_data.feat
        else:
            raise ValueError("IllegalArgument: axis = {}. Accepted values are obs and feat")
    

    def get_annotation(self, key:str, axis:str="obs"):
        if axis == "obs":
            return self.cnv_data.obs[key]
        elif axis == "feat":
            return self.cnv_data.feat[key]
        else:
            raise ValueError("IllegalArgument: axis = {}. Accepted values are obs and feat")

    def count(self):
        return self.cnv_data.n_obs

    #da sistemare: basterebbe usare il meccanismo delle views
    def sort_rows(self, by:Union[list, np.array]):
        return Sample(self.cnv_data.sort_rows(by), self.name)
    
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
                b: Optional[float] = None, use_highly_variable: Optional[bool] = False):
                umap_result = umap(self.cnv_data, n_neighbors, n_components, metric, metric_kwds, min_dist, spread, maxiter, alpha, gamma,
                                fast, negative_sample_rate, local_connectivity, init_pos, random_state, a, b, use_highly_variable)
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
    # normalization
    def zscore(self, zero_center: bool = True, max_value: Optional[float] = None):
        X_, mean, var = zscore(self.cnv_data.X, zero_center, max_value)
        self.cnv_data.uns['zscore'] = {}
        self.cnv_data.uns['zscore']['X'] = X_
        self.cnv_data.feat['mean'] = mean
        self.cnv_data.feat['var'] = var

    def get_zscore(self):
        return self.cnv_data.uns['zscore']['X']

    # clustering
    def nk_clust(self, method:str, metric:Optional[Union[str, None]]="euclidean", 
            linkage:Optional[Union[str, None]]=None, embeddings:Optional[Union[str, None]]=None, n_comps: Optional[int] = None,
            min_k:Optional[int]=2, max_k:Optional[int]=15, index:Optional[Union[str, List[str]]]="all",
            n_jobs:Optional[int]=1):
        """
        if embeddings is not None:
            if embeddings == EMBEDDINGS.PCA:
                if "pca" not in self.cnv_data.uns.keys():
                    raise ValueError(
                        'Did not find cnv.uns[\'pca\']. '
                        'Consider running `Sample.pca()` first.'
                    )
                else:
                    data = self.get_pca("X")
            elif embeddings == EMBEDDINGS.UMAP:
                if "umap" not in self.cnv_data.uns.keys():
                    raise ValueError(
                        'Did not find cnv.uns[\'umap\']. '
                        'Consider running `Sample.umap()` first.'
                    )
                else:
                    data = self.get_umap("X")
            else:
                raise ValueError("Accepted values for embeggings are 'umap', 'pca'")
        else: 
            data = (self.cnv_data.X)
        """
        results = nk_clust(self.cnv_data , method, metric, linkage, embeddings, n_comps, min_k, max_k, index, n_jobs) 
        scores_df = pd.DataFrame(results)
        scores_df.index.name = "k"
        self.cnv_data.uns['nk_clust'] = scores_df
        return scores_df

    def cluster_benchmark(self, methods:Sequence[object], cluster_configurations:ClusterConfig, labels:Union[list, np.array], 
                        embeddings:Optional[Union[str, None]]=None, n_comps: Optional[int] = None, n_jobs:Optional[int]=1):
        ranking = clust_accuracy_ranking(self.cnv_data, cluster_configurations, labels, embeddings, n_comps, n_jobs)
        self.cnv_data.uns['cluster_benchmark'] = ranking
        return ranking

    def cluster(self, method:str, embeddings:Optional[Union[str, None]]=None, n_comps: Optional[int] = None, **kwargs):
        res = cluster(self.cnv_data, method, embeddings, n_comps, **kwargs)
        self.cnv_data.uns['cluster_model'] = res
        self.cnv_data.obs['cluster'] = res.labels_
        return res

    def get_clusterer(self):
        return self.cnv_data.uns['cluster_model']
    
    def get_clusters(self):
        return self.cnv_data.obs['cluster']

    def get_cluster_labels(self):
        return self.cnv_data.obs['cluster'].values

    #Multi-sample analysis
    #def co_cluster(self, samples:list, method:str, embeddings:Optional[Union[str, None]]=None, **kwargs):
    #    multisample = MultiSample.from_list(samples)
    #    multisample.cluster(method, embeddings, **kwargs)
    #    return multisample
    
    #def multi_sample_score(self, samples:list, n_jobs:int=1, verbose:bool=False):
    #    multisample = MultiSample.from_list(samples)
    #    return multisample.het_score(n_jobs, verbose)

    #Filter functions
    def filter(self, key:str, method:str,
                    value:Union[int, float, Sequence[Union[float, int]], Tuple[Union[int, float], Union[int, float]]], 
                    percentile:bool=False, inplace:bool=False):
        if key not in self.cnv_data.obs.columns:
            raise ValueError('Did not find cnv_data.obs[\'{}\']. '.format(key))
        if isinstance(value, tuple):
            return self._filter_interval(key, method, value, inplace=inplace)
        elif isinstance(value, int) or isinstance(value, float):
            return self._filter_value(key, method, value, percentile=percentile, inplace=inplace)
        elif isinstance(value, list) or isinstance(value, np.array):
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
            subset = self.cnv_data[(annotation > value[0]) & (annotation < value[1]), :]
        elif method ==  _FILTER_INTERVAL_METHODS_['IN_EQ']:
            subset = self.cnv_data[(annotation >= value[0]) & (annotation <= value[1]), :]
        elif method ==  _FILTER_INTERVAL_METHODS_['OUT']:
            subset = self.cnv_data[(annotation < value[0]) | (annotation > value[1]), :]
        elif method ==  _FILTER_INTERVAL_METHODS_['OUT_EQ']:
            subset = self.cnv_data[(annotation <= value[0]) | (annotation >= value[1]), :]
        #index = subset.obs_names
        #annotations = self.get_annotations()
        #new_annotations = annotations[annotations.index.isin(index)]
        #s = Sample(subset, self.name)
        #for c in new_annotations.columns:
        #    s.add_annotation(new_annotations[c], c)
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

    def plot_annotation_dist(self, ann:str, axis:str="obs", grid:bool=False, quantiles:List[float]=None, figsize:Tuple[int, int]=None, 
                                outpath:str=None, **kwargs):
        Drawer.draw('dist', self.get_annotation(key=ann, axis=axis), grid, quantiles, figsize, outpath, **kwargs)

    """
    def heatmap(self, method:str ='ward', metric:str ='euclidean', outpath:str=None,
                    vmin:int = 0, vmax:int = 12, vcenter:int=2, figsize:Tuple[int, int]=(37, 21), fontsize:int=16):
        Drawer.draw('heatmap', self.get_cnv_matrix(), self.get_boundaries(),  outpath=outpath, sample=self.name, 
            vmin = vmin, vmax=vmax, vcenter=vcenter, figsize=figsize, fontsize=fontsize)

    
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
        
    def plot_umap(self, outpath:str=None, 
                figsize:Tuple[int, int]=None, **kwargs):
        if 'umap' in self.cnv_data.uns:
            Drawer.draw('scatter', data=self.get_umap('X'), title = 'UMAP embeddings', x_label='X', y_label='Y', 
                outpath=outpath, figsize = figsize, **kwargs)
        else:
            logg.error("{} object has no attribute 'umap'".format(self.cnv_data.uns))

    def plot_zscore(self, outpath:str=None, figsize:Tuple[int, int]=None, **kwargs):
        if 'zscore' in self.cnv_data.uns:
            Drawer.draw('heatmap', data=self.get_zscore(), boundaries=self.get_boundaries(),  title = 'Z-Score', outpath=outpath, row_cluster=True, 
                     vcenter=0, figsize=figsize, **kwargs)
        else:
            logg.error("{} object has no attribute 'zscore'".format(self.cnv_data.uns))

    def plot_clusters(self, plot:str="scatter", outpath:str=None, figsize:Tuple[int, int]=None, **kwargs):
        if 'cluster' in self.cnv_data.obs.columns:  
            if plot == "scatter":
                if 'umap' in self.cnv_data.uns:
                    projection = self.get_umap('X')
                else:
                    projection = umap(self.cnv_data)
                Drawer.draw('scatter', data=projection, title = 'Clusters', x_label='X', y_label='Y', outpath=outpath,  
                    figsize = figsize, labels=self.get_clusters(), legend=True, **kwargs)
            elif plot == "heatmap":
                s = self.sort_rows(by="cluster")
                Drawer.draw('heatmap', data=s.get_cnv_dataframe(), boundaries=s.get_boundaries(),  title = 'Cluster heatmap', outpath=outpath,  
                    labels=s.cnv_data.obs['cluster'], figsize=figsize, legend=True, **kwargs)
        else:
            logg.error("{} object has no column 'cluster'".format(self.cnv_data.feat))

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
                        labels=self.cnv_data.obs['cluster'], legend=True, **kwargs)
                else:
                    logg.error("Cluster model has no attribute 'children'. Consider executing agglomerative clustering.")
            else:
                logg.error("{} object has no column 'cluster'".format(self.cnv_data.feat))
        else:
            Drawer.draw('heatmap', data=self.get_cnv_dataframe(), row_cluster=True, boundaries=self.get_boundaries(),  title = 'Dataset dendrogram & heatmap', outpath=outpath,  
                         **kwargs)
    



    
    
    


