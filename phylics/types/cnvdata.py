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
# cnvdata.py: this module implements phylics data type
# ==========================================================================

from ._utils import _normalize_indices, Index,  Index1D
import collections.abc as cabc
import numpy as np
import pandas as pd
from copy import copy
from typing import Union, List, Tuple, Optional, Mapping, Iterable, Any, Sequence

from .views import (
    ArrayView,
    DictView,
    DataFrameView,
    as_view,
    _resolve_idxs,
)

from anndata.compat import (
    _slice_uns_sparse_matrices,)

__all__ = ["CnvData"]

class CnvData:
    def __init__(self, X:Union[np.ndarray, pd.DataFrame, "CnvData"], 
            obs_names:Optional[Union[Sequence[Union['str','int']], str]] = 'auto', 
            feat_names:Optional[Union[Sequence[Union['str','int']], str]] = 'auto', 
            obs:Optional[Union[pd.DataFrame, Mapping[str, Iterable[Any]]]]=None, 
            feat:Optional[Union[pd.DataFrame, Mapping[str, Iterable[Any]]]]=None, 
            uns:Optional[Mapping[str, Any]]=None, raw:Optional[Mapping[str, Any]]=None, 
            asview: bool = False,
            oidx: Index1D = None,
            fidx: Index1D = None
            ):
        if asview == True:
            if not isinstance(X, CnvData):
                raise ValueError("`X` has to be an CnvData object.")
            self._init_as_view(X, oidx, fidx)
        else:
            self._is_view = False
            if isinstance(obs_names, str):
                if obs_names == 'auto':
                    self.obs_names = np.arange(0, X.shape[0])
                elif obs_names == 'row_names':
                    if isinstance(X, pd.DataFrame):
                        self.obs_names = list(X.index.values)
                    else:
                        raise ValueError("obs_names = 'raw_names' requires X of type pd.DataFrame")
                elif obs_names == 'obs_names':
                    if obs != None:
                        if isinstance(obs, pd.DataFrame):
                            self.obs_names = list(obs.index.values)
                        elif isinstance(obs, cabc.Mapping):
                            self.obs_names = list(obs.keys())
                        else:
                            raise ValueError("obs must be one of type pd.DataFrame or Mappin")
                    else:
                        raise ValueError("obs_names = 'obs_names' requires obs of type pd.DataFrame or Mapping") 
                else:
                    raise ValueError("Illegal argument for obs_names ({})".format(obs_names))               
            elif isinstance(obs_names, cabc.Sequence):
                self.obs_names = obs_names
            else:
                raise ValueError("obs_names must be of type Sequence or string")
            
            if isinstance(feat_names, str):
                if feat_names == 'auto':
                    self.feat_names = np.arange(0, X.shape[1])
                elif feat_names == 'col_names':
                    if isinstance(X, pd.DataFrame):
                        self.feat_names = list(X.columns.values)
                    else:
                        raise ValueError("feat_names = 'col_names' requires X of type pd.DataFrame")
                elif feat_names == 'feat_names':
                    if feat != None:
                        if isinstance(feat, pd.DataFrame):
                            self.feat_names = list(obs.index.values)
                        elif isinstance(feat, cabc.Mapping):
                            self.feat_names = list(feat.keys())
                        else:
                            raise ValueError("feat must be of type pd.DataFrame or Mapping")
                    else:
                        raise ValueError("feat_names = 'obs_names' requires obs of type pd.DataFrame or Mapping") 
                else:
                    raise ValueError("Illegal argument for feat_names ({})".format(feat_names))               
            elif isinstance(feat_names, cabc.Sequence):
                self.feat_names = feat_names
            else:
                raise ValueError("feat_names must be of type Sequence or string")
            

            if len(self.obs_names) != X.shape[0] or len(self.feat_names) != X.shape[1]:
                raise ValueError("X: expected shape ({} {}) but got ({})".format(len(self.obs_names), len(self.feat_names), X.shape))

            if isinstance(X, pd.DataFrame):
                X_arr = X.to_numpy()
                self.X = X_arr
            elif isinstance(X, np.ndarray):
                self.X = X 
            else:
                raise ValueError("X must be one of (np.ndarray, pd.DataFrame)")

            if obs is None:
                self.obs = pd.DataFrame(index=self.obs_names)
            elif isinstance(obs, pd.DataFrame):
                self.obs = obs
            elif isinstance(obs, cabc.Mapping):
                self.obs = pd.DataFrame.from_dict(obs, orient='index')
            else:
                raise ValueError("obs must be one of (pd.DataFrame, Mapping[str, Iterable[Any]])")
        
            if feat is None:
                self.feat = pd.DataFrame(index=self.feat_names)
            elif isinstance(feat, pd.DataFrame):
                self.feat = feat
            elif isinstance(feat, cabc.Mapping):
                self.feat = pd.DataFrame.from_dict(feat, orient='index')
            else:
                raise ValueError("feat must be one of (pd.DataFrame, Mapping[str, Iterable[Any]])")

            if len(self.obs) != X.shape[0] or len(self.feat) != X.shape[1]:
                raise ValueError("X: expected shape ({} {}) but got ({})".format(len(self.obs), len(self.feat), X.shape))

            self.n_obs = len(self.obs_names)
            self.n_feat = len(self.feat_names)
            if uns == None:
                self.uns = {}
            else:
                self.uns = uns 
            if raw == None:
                self.raw = {}
                self.raw['X'] = self.X
                self.raw['obs'] = self.obs
                self.raw['feat'] = self.feat
            else:
                self.raw = raw
            self.shape = (self.n_obs, self.n_feat)
            

    def _init_as_view(self, data_ref:"CnvData" , oidx: Index, fidx: Index):
        self._is_view = True
        if isinstance(oidx, (int, np.integer)):
            oidx = slice(oidx, oidx + 1, 1)
        if isinstance(fidx, (int, np.integer)):
            fidx = slice(fidx, fidx + 1, 1)
        if data_ref._is_view:
            prev_oidx, prev_fidx = data_ref._oidx, data_ref._fidx
            data_ref = data_ref._adata_ref
            oidx, fidx = _resolve_idxs((prev_oidx, prev_fidx), (oidx, fidx), data_ref)
        self.oidx = oidx
        self.idx = fidx
    
        obs_sub = data_ref.obs.iloc[oidx]
        feat_sub = data_ref.feat.iloc[fidx]
        self.X = data_ref.to_df().iloc[oidx, fidx].to_numpy()
        uns_new = _slice_uns_sparse_matrices(
            copy(data_ref.uns), self.oidx, data_ref.n_obs
        )
        # fix categories
        self._remove_unused_categories(data_ref.obs, obs_sub, uns_new)
        self._remove_unused_categories(data_ref.feat, feat_sub, uns_new)
        # set attributes
        self.obs = DataFrameView(obs_sub, view_args=(self, "obs"))
        self.obs_names = list(self.obs.index.values)
        self.feat = DataFrameView(feat_sub, view_args=(self, "feat"))
        self.feat_names = list(self.feat.index.values)
        self.uns = DictView(uns_new, view_args=(self, "uns"))
        self.n_obs = len(self.obs)
        self.n_feat = len(self.feat)

        if data_ref.raw is not None:
            self.raw = data_ref.raw
        else:
            self.raw = None
        self.shape = (self.n_obs, self.n_feat)
    

    def __repr__(self) -> str:
        descr = f"CnvData object with n_obs × n_vars = {self.n_obs} × {self.n_feat}"
        for attr in ["obs", "feat", "uns","raw"]:
            x = getattr(self, attr)
            keys = []
            if isinstance(x, pd.DataFrame):
                keys = x.to_dict().keys()
            else:
                keys = x.keys()
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(list(keys))[1:-1]}"
        return descr

    def __eq__(self, other):
        """Equality testing"""
        raise NotImplementedError(
            "Equality comparisons are not supported for CnvData objects, "
            "instead compare the desired attributes."
        )
    
    def __getitem__(self, index: Index) -> "CnVData":
        """Returns a sliced view of the object."""
        oidx, fidx = self._normalize_indices(index)
        return CnvData(self, oidx=oidx, fidx=fidx, asview=True)
    
    def _normalize_indices(self, index: Optional[Index]) -> Tuple[slice, slice]:
        return _normalize_indices(index, self.obs.index, self.feat.index)
    
    def _remove_unused_categories(self, df_full, df_sub, uns):
        from pandas.api.types import is_categorical

        for k in df_full:
            if not is_categorical(df_full[k]):
                continue
            all_categories = df_full[k].cat.categories
            df_sub[k].cat.remove_unused_categories(inplace=True)
            # also correct the colors...
            color_key = f"{k}_colors"
            if color_key not in uns:
                continue
            color_vec = uns[color_key]
            if np.array(color_vec).ndim == 0:
                # Make 0D arrays into 1D ones
                uns[color_key] = np.array(color_vec)[(None,)]
            elif len(color_vec) != len(all_categories):
                # Reset colors
                del uns[color_key]
            else:
                idx = np.where(np.in1d(all_categories, df_sub[k].cat.categories))[0]
                uns[color_key] = np.array(color_vec)[(idx,)]
    
    def drop_obs(self, select:Union[Sequence[Union[int, str, np.integer]], int, str], inplace:bool=False):
        X_df = pd.DataFrame(self.X, index=self.obs_names, columns=self.feat_names)
        X_df = X_df.drop(select, axis=0)
        obs_filtered = self.obs.drop(select, axis=0)
        if inplace == True:
            self.X = X_df.to_numpy()
            self.obs = obs_filtered
            self.obs_names = self.obs.index.values
            self.n_obs = len(self.obs_names)
            return self
        else:
            return CnvData(X=X_df, obs_names = list(obs_filtered.index.values), feat_names= 'col_names', obs=obs_filtered, feat=self.feat, uns = self.uns)

    def drop_feat(self, select:Union[Sequence[Union[int, str, np.integer]], int, str], inplace:bool=False):
        X_df = pd.DataFrame(self.X, index=self.obs_names, columns=self.feat_names)
        X_df = X_df.drop(select, axis=1)
        feat_filtered = self.feat.drop(select)
        if inplace == True:
            self.X = X_df.to_numpy()
            self.feat = feat_filtered
            self.feat_names = self.feat.index.values
            self.n_feat = len(self.feat_names)
            return self
        else:
            return CnvData(X=X_df, feat_names= list(feat_filtered.index.values), obs_names='row_names', obs=self.obs, feat=feat_filtered, uns = self.uns)

    def to_df(self):
        return pd.DataFrame(self.X, index=self.obs_names, columns=self.feat_names)
    
    def select_obs(self, select:Sequence[Union[int, str, np.integer]], inplace:bool=False):
        X_df = pd.DataFrame(self.X, index=self.obs_names, columns=self.feat_names)
        X_df = X_df.loc[select]
        obs_filtered = self.obs.loc[select]
        shape = X_df.shape   
        if inplace == True:
            self.X = X_df.to_numpy()
            self.obs = obs_filtered
            self.obs_names = self.obs.index.values
            self.n_obs = len(self.obs_names)
            self.shape = shape
            return self
        else:
            return CnvData(X=X_df, obs_names = list(obs_filtered.index.values), feat_names= 'col_names', obs=obs_filtered, feat=self.feat, uns = self.uns)
        
    def copy(self) -> "CnvData":
            return CnvData(
                X=self.X.copy(),
                obs=self.obs.copy(),
                feat=self.feat.copy(),
                obs_names=self.obs_names,
                feat_names=self.feat_names,
                uns=self.uns.copy() if self.uns is not None else None,
                raw=self.raw.copy() if self.raw is not None else None,
            )
