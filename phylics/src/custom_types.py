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
# custom_types.py: this module implements phylics custom types
# ==========================================================================
__all__ = ['CNVS']

import numpy as np
import pandas as pd
from typing import Union, List

import umap

class CNVS:
    def __init__(self, cnvs_df: pd.DataFrame):
        self.cnvs_df = cnvs_df
        self.cnvs = self.cnvs_df.drop(['CHR', 'START', 'END'], axis=1).transpose()
        self.boundaries = self.cnvs_df[['CHR', 'START', 'END']].copy()
        self.cells = self.cnvs.index.values

    
    def __repr__(self):
        return repr(self.cnvs_df)

    def __str__(self):
        return repr(self.cnvs_df)
    
    
    @classmethod
    def from_file(cls, cnvs_df:str):
        return cls(pd.read_csv(cnvs_df, sep="\t"))

    def get_dataframe(self):
        return self.cnvs_df
    def head(self, n:int = 5):
        print(self.cnvs_df.head(n))
    def get_cnvs(self):
        return self.cnvs
    def get_boundaries(self):
        return self.boundaries
    def get_cell_names(self):
        return self.cells
    def get_mads(self):
        return self.cnvs.mad(axis=1).to_dict()
    def get_cell_cnvs(self, cellid:str):
        return self.cnvs_df[cellid].values
    def transpose_cnvs(self):
        return self.cnvs.transpose()   
    
    def get_cells(self, cells:List[str]):
        cnvs_df = self.cnvs_df[np.append(['CHR', 'START', 'END'], cells)]
        cnvs = CNVS(cnvs_df)
        return cnvs
    
    def drop_cells(self, cells:List[str]):
        cnvs_df = self.cnvs_df.drop(cells, axis=1)
        cnvs = CNVS(cnvs_df)
        return cnvs
    
class LookUpTable:
    def __init__(self, mapping:dict):
        self.lut = mapping
    
    def __repr__(self):
        return repr(pd.Series(self.lut).to_frame().T)

    def __str__(self):
        return repr(pd.Series(self.lut).to_frame().T)

    @classmethod
    def from_lists(cls, ordered_ids:List, ordered_values:List):
        lut = dict(zip(ordered_ids, ordered_values))
        return cls(lut)
    
    def get(self):
        return self.lut

    def items(self):
        return self.lut.items()
    
    def keys(self):
        return list(self.lut.keys())

    def values(self):
        return list(self.lut.values())
    
    def get_by_key(self, key:str):
        return self.lut[key]
    
    def get_by_value(self, value:str):
        values = []
        for v in self.lut.values():
            if v == value:
                values.append(v)
        return values

    def get_by_keys(self, keys:list):
        return {cell : v for cell, v in self.lut.items() if cell in keys}

    def get_by_values(self, values:list):
        return {cell : v for cell, v in self.lut.items() if v in values}
    
    def drop_by_key(self, key:Union[str, int]):
        return {cell : v for cell, v in self.lut.items() if cell != key}
    
    def drop_by_value(self, value:Union[str, int, float]):
        return {cell : v for cell, v in self.lut.items() if v != value}

    def drop_by_keys(self, keys:list):
        return {cell : v for cell, v in self.lut.items() if cell not in keys}

    def drop_by_values(self, values:list):
        return {cell : v for cell, v in self.lut.items() if v not in values}

class Reducer:
    @staticmethod
    def umap_(X, **kwargs):
        reducer = umap.UMAP(**kwargs)
        return reducer.fit_transform(X)
    
    @staticmethod
    def pca_(X, **kwargs):
        return X
        