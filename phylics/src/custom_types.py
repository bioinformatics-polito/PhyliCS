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


import pandas as pd
from typing import Union, List

class CNVS:
    def __init__(self, cnvs_df: pd.DataFrame):
        self.cnvs_df = cnvs_df
        self.cnvs = self.cnvs_df.drop(['CHR', 'START', 'END'], axis=1).transpose()
        self.boundaries = self.cnvs_df[['CHR', 'START', 'END']].copy()
        self.cells = self.cnvs.index.values
    
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
    def get_cells(self):
        return self.cells
    def get_mads(self):
        return self.cnvs.mad(axis=1).to_dict()

    def transpose_cnvs(self):
        return self.cnvs.transpose()   
    
    def drop_cells(self, cells:List[str]):
        cnvs_df = self.cnvs_df.drop(cells, axis=1)
        cnvs = CNVS(cnvs_df)
        return cnvs
    
