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

import collections.abc as cabc
import numpy as np
import pandas as pd
from typing import Union, List, Tuple, Optional, Mapping, Iterable, Any, Sequence

__all__ = ["CnvData"]

class CnvData:
    def __init__(self, X:Union[np.ndarray, pd.DataFrame], 
            obs_names:Optional[Union[Sequence[Union['str','int']], str]] = 'auto', 
            feat_names:Optional[Union[Sequence[Union['str','int']], str]] = 'auto', 
            obs:Optional[Union[pd.DataFrame, Mapping[str, Iterable[Any]]]]=None, 
            feat:Optional[Union[pd.DataFrame, Mapping[str, Iterable[Any]]]]=None, 
            uns:Optional[Mapping[str, Any]]=None, raw:Optional[Mapping[str, Any]]=None, 
            shape:Tuple[int, int]=None):
        
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

        if shape == None:
            shape = X.shape
        self.shape = shape
    
    def drop_obs(self, select:Union[Sequence[Union[int, str, np.integer]], int, str], inplace:bool=False):
        X_df = pd.DataFrame(self.X, index=self.obs_names, columns=self.feat_names)
        X_df = X_df.drop(select, axis=0)
        obs_filtered = self.obs.drop(select, axis=0)
        shape = X_df.shape   
        if inplace == True:
            self.X = X_df.to_numpy()
            self.obs = obs_filtered
            self.obs_names = self.obs.index.values
            self.n_obs = len(self.obs_names)
            self.shape = shape
            return self
        else:
            return CnvData(X=X_df, obs_names = list(obs_filtered.index.values), feat_names= 'col_names', obs=obs_filtered, feat=self.feat, uns = self.uns, shape=shape)

    def drop_feat(self, select:Union[Sequence[Union[int, str, np.integer]], int, str], inplace:bool=False):
        X_df = pd.DataFrame(self.X, index=self.obs_names, columns=self.feat_names)
        X_df = X_df.drop(select, axis=1)
        feat_filtered = self.feat.drop(select)
        shape = X_df.shape   
        if inplace == True:
            self.X = X_df.to_numpy()
            self.feat = feat_filtered
            self.feat_names = self.feat.index.values
            self.n_feat = len(self.feat_names)
            self.shape = shape
            return self
        else:
            return CnvData(X=X_df, feat_names= list(feat_filtered.index.values), obs_names='row_names', obs=self.obs, feat=feat_filtered, uns = self.uns, shape=shape)

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
            return CnvData(X=X_df, obs_names = list(obs_filtered.index.values), feat_names= 'col_names', obs=obs_filtered, feat=self.feat, uns = self.uns, shape=shape)
        


   

