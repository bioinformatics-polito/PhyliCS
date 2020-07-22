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

from .custom_types import CNVS
from .drawer import Drawer
import random
import pandas as pd
import numpy as np
from typing import Union, Tuple, List

_FILTER_BY_SINGLE_PLOIDY_METHODS_ = {  
                                'EQ' :'eq', 
                                'LT' : 'lt',
                                'GT' : 'gt', 
                                'LT_EQ' : 'lt_eq', 
                                'GT_EQ' : 'gt_eq'}
_FILTER_BY_RANGE_PLOIDY_METHODS_ = {
                                'IN' : 'in', 
                                'OUT' : 'out', 
                                'IN_EQ' : 'in_eq', 
                                'OUT_EQ' : 'out_eq'}
_FILTER_METHODS_ = {
    'PERC':'percentile',
    'VALUE':'value'}



class Sample:

    def __init__(self, cnvs_dataframe:CNVS, cell_mad:Union[dict, str]=None, cell_ploidy:Union[dict, str]=None, 
                        cell_coverage:Union[dict, str]=None, sample_name:str="sample"):      
        self.df = cnvs_dataframe
        self.cnvs = self.df.get_cnvs()
        self.boundaries = self.df.get_boundaries()
        self.cells = self.df.get_cells()
        self.name = sample_name
        
        if cell_ploidy != None:
            if isinstance(cell_ploidy, dict):
                self.cell_ploidy = cell_ploidy
            elif isinstance(cell_ploidy, str):
                self.cell_ploidy = pd.read_csv(cell_ploidy, header=None, squeeze=True, index_col=0, sep="\t").to_dict()
        else:
            self.cell_ploidy = None
        
        if cell_coverage != None:
            if isinstance(cell_coverage, dict):
                self.cell_coverage = cell_coverage
            elif isinstance(cell_ploidy, str):
                self.cell_coverage = pd.read_csv(cell_coverage, header=None, squeeze=True, index_col=0, sep="\t").to_dict()
        else:
            self.cell_coverage = None
        
        if cell_mad == None:
            self.cell_mad = self.df.get_mads()
        elif isinstance(cell_mad, dict):
            self.cell_mad = cell_mad
        elif isinstance(cell_mad, str):
            self.cell_mad = pd.read_csv(cell_mad, header=None, squeeze=True, index_col=0, sep="\t").to_dict()

    @classmethod
    def from_file(cls, cnvs_dataframe:str, cell_mad:Union[dict, str]=None, cell_ploidy:Union[dict, str]=None, 
                        cell_coverage:Union[dict, str]=None, sample_name:str="sample"):
        return cls(CNVS.from_file(cnvs_dataframe), cell_mad, cell_ploidy, cell_coverage, sample_name)
        

    def set_cell_mad(self, cell_mad:Union[dict, str]):
        if isinstance(cell_mad, dict):
            self.cell_mad = cell_mad
        elif isinstance(cell_mad, str):
            self.cell_mad = pd.read_csv(cell_mad, header=0, squeeze=True, index_col=0, sep="\t").to_dict()

    def set_cell_ploidy(self, cell_ploidy:Union[dict, str]):
        if isinstance(cell_ploidy, dict):
            self.cell_ploidy = cell_ploidy
        elif isinstance(cell_ploidy, str):
            self.cell_ploidy = pd.read_csv(cell_ploidy, header=0, squeeze=True, index_col=0, sep="\t").to_dict()

    def set_cell_coverage(self, cell_coverage:Union[dict, str]):
        if isinstance(cell_coverage, dict):
            self.cell_coverage = cell_coverage
        elif isinstance(cell_coverage, str):
            self.cell_coverage = pd.read_csv(cell_coverage, header=0, squeeze=True, index_col=0, sep="\t").to_dict()        
    
    def set_name(self, sample_name:str):
        self.name = sample_name
    
    def get_dataframe(self):
        return self.df 
    
    def get_cnvs(self):
        return self.cnvs
    
    def get_cells(self):
        return self.cells

    def get_boundaries(self):
        return self.boundaries

    def get_name(self):
        return self.name

    def get_cell_ploidy(self):
        return self.cell_ploidy
    
    def get_cell_coverage(self):
        return self.cell_coverage
    
    def get_cell_mad(self):
        return self.cell_mad

    def count(self):
        return len(self.cells)

    def plot_ploidy_dist(self, kde:bool=True, rug:bool=False, vertical:bool=False, grid:bool=False,
                quantiles:List[float]=None, axlabel:str=None, label:str=None, figsize:Tuple[int, int]=None, outpath:str=None):
        Drawer.draw('dist', list(self.cell_ploidy.values()), kde, rug, vertical, grid, quantiles, axlabel, label, figsize, outpath)
    
    def plot_mad_dist(self, kde:bool=True, rug:bool=False, vertical:bool=False, grid:bool=False,
                quantiles:List[float]=None, axlabel:str=None, label:str=None, figsize:Tuple[int, int]=None, outpath:str=None):
        Drawer.draw('dist', list(self.cell_mad.values()), kde, rug, vertical, grid, quantiles, axlabel, label, figsize, outpath)

    def plot_cell_coverage_dist(self, kde:bool=True, rug:bool=False, vertical:bool=False, grid:bool=False,
                quantiles:List[float]=None, axlabel:str=None, label:str=None, figsize:Tuple[int, int]=None, outpath:str=None):
        Drawer.draw('dist', list(self.cell_coverage.values()), kde, rug, vertical, grid, quantiles, axlabel, label, figsize, outpath)

    def heatmap(self, method:str ='ward', metric:str ='euclidean', outpath:str=None,
                    vmin:int = 0, vmax:int = 12, vcenter:int=2, figsize:Tuple[int, int]=(37, 21), fontsize:int=16):
        Drawer.draw('heatmap', self.cnvs, self.boundaries, method, metric,  outpath=outpath, sample=self.name, 
            vmin = vmin, vmax=vmax, vcenter=vcenter, figsize=figsize, fontsize=fontsize)
    
    #Filter functions

    def filter_by_ploidy(self, method:str='eq', ploidy:Union[float, int, Tuple]=2.0):
        """
        Filter out group of cells based on their ploidy 
        -----------------------------------------------------
        method: the method according which the cells are filtered.
                    When a single value is specified, acceptable values are:
                        - 'eq', 'lt', 'gt', 'lt_eq', 'gt_eq' , to filter out cells with ploidy equal, 
                        less than, greater than, less or equal than, greater or equal then the provided values
                    When a ploidy interval is provided, acceptable values are:
                        - 'in', 'out', 'in_eq', 'out_eq', to filter out cells with ploidy whithin or out of
                        the interval (extremes excluded and included, respectively).
        ploidy: single value or interval of cell_ploidy
        """
        if (method not in list(_FILTER_BY_SINGLE_PLOIDY_METHODS_.values())) and (method not in list(_FILTER_BY_RANGE_PLOIDY_METHODS_.values())):
            raise ValueError("filter_by_ploidy: method must be one of %r or %r."
            % list(_FILTER_BY_SINGLE_PLOIDY_METHODS_.values()), list(_FILTER_BY_RANGE_PLOIDY_METHODS_.values()))
        cells = []
        if isinstance(ploidy, float) or isinstance(ploidy, int):
            if method == _FILTER_BY_SINGLE_PLOIDY_METHODS_['EQ']:
                for cell, p in self.cell_ploidy.items():
                    if p == ploidy:
                        cells.append(cell)
            elif method == _FILTER_BY_SINGLE_PLOIDY_METHODS_['LT']:
                for cell, p in self.cell_ploidy.items():
                    if p < ploidy:
                        cells.append(cell)
            elif method == _FILTER_BY_SINGLE_PLOIDY_METHODS_['LT_EQ']:
                for cell, p in self.cell_ploidy.items():
                    if p <= ploidy:
                        cells.append(cell)
            elif method == _FILTER_BY_SINGLE_PLOIDY_METHODS_['GT']:
                for cell, p in self.cell_ploidy.items():
                    if p > ploidy:
                        cells.append(cell)
            elif method == _FILTER_BY_SINGLE_PLOIDY_METHODS_['GT_EQ']:
                for cell, p in self.cell_ploidy.items():
                    if p >= ploidy:
                        cells.append(cell)
            else:
                raise ValueError("filter_by_ploidy: with a single value, method must be one of  %r." % list(_FILTER_BY_SINGLE_PLOIDY_METHODS_.values()))
        elif isinstance(ploidy, Tuple):
            if ploidy[0] >= ploidy[1]:
                raise ValueError("filter_by_ploidy: the first element of the ploidy tuple must be less than the second element.")
            if method == _FILTER_BY_RANGE_PLOIDY_METHODS_['IN']:
                for cell, p in self.cell_ploidy.items():
                    if p > ploidy[0] and p < ploidy[1]:
                        cells.append(cell)
            elif method ==  _FILTER_BY_RANGE_PLOIDY_METHODS_['IN_EQ']:
                for cell, p in self.cell_ploidy.items():
                    if p >= ploidy[0] and p <= ploidy[1]:
                        cells.append(cell)
            elif method ==  _FILTER_BY_RANGE_PLOIDY_METHODS_['OUT']:
                for cell, p in self.cell_ploidy.items():
                    if p < ploidy[0] or p > ploidy[1]:
                        cells.append(cell)
            elif method ==  _FILTER_BY_RANGE_PLOIDY_METHODS_['OUT_EQ']:
                for cell, p in self.cell_ploidy.items():
                    if p <= ploidy[0] or p >= ploidy[1]:
                        cells.append(cell)
            else:
                raise ValueError("filter_by_ploidy: with a range of values, method must be one of  %r." % list(_FILTER_BY_RANGE_PLOIDY_METHODS_.values()))
        else:
            raise TypeError("filter_by_ploidy: the ploidy must be of type int or float or tuple (e.g. (1.7, 2.0))")
        cnvs_df = self.df.drop_cells(cells)
        filtered_cell_mad = { cell: self.cell_mad[cell] for cell in cells }
        if self.cell_ploidy != None:
            filtered_cell_ploidy = { cell: self.cell_ploidy[cell] for cell in cells }
        else:
            filtered_cell_ploidy = None
        if self.cell_coverage != None:
            filtered_cell_coverage = { cell: self.cell_coverage[cell] for cell in cells } 
        else:
            filtered_cell_coverage = None
        sample = Sample(cnvs_df, filtered_cell_mad, filtered_cell_ploidy, filtered_cell_coverage, self.name)
            
        return sample
    
    def filter_by_mad(self, method:str='percentile', threshold:float='0.9'):
        """
        Filter out cells with MAD above a given threshold
        -----------------------------------------------------
        method: the method according which the cells are filtered. 
                        'percentile' and 'value' are accepted values.
        threshold: the MAD value above which cells are filtered out. 
                    This values has two interpretations: if method is 'percentile', 
                    means that cells with MAD up to that percentile are kept; if 
                    method is 'value', then this parameter is interpreted as the exact
                    MAD value above which cells must be removed.
        """
        if method not in list(_FILTER_METHODS_.values()):
            raise ValueError("filter_by_mad: method must be one of %r." % list(_FILTER_METHODS_.values()))
        if method == _FILTER_METHODS_['PERC']:
            perc = np.quantile(list(self.cell_mad.values()), float(threshold))
            filtered_cell_mad = { cell: mad for cell, mad in self.cell_mad.items() if mad <= perc}
        elif method == _FILTER_METHODS_['VALUE']:
            filtered_cell_mad = { cell: mad for cell, mad in self.cell_mad.items() if mad <= threshold}
        cells = filtered_cell_mad.keys()
        if self.cell_ploidy != None:
            filtered_cell_ploidy = { cell: self.cell_ploidy[cell] for cell in cells }
        else:
            filtered_cell_ploidy = None
        if self.cell_coverage != None:
            filtered_cell_coverage = { cell: self.cell_coverage[cell] for cell in cells }
        else:
            filtered_cell_coverage = None
        #TODO
        # fix = cells must be kept and not removed
        cnvs_df = self.df.drop_cells(cells)
        sample = Sample(cnvs_df, filtered_cell_mad, filtered_cell_ploidy, filtered_cell_coverage, self.name)
        
        return sample

    def filter_by_coverage(self, method:str='percentile', threshold:float='0.9'):
        """
        Filter out cells with a read count below a given threshold
        -----------------------------------------------------
        method: the method according which the cells are filtered. 
                        'percentile' and 'value' are accepted values.
        threshold: the read count below which cells are filtered out. 
                    This values has two interpretations: if method is 'percentile', 
                    means that cells with read count up to that percentile are filtered out; if 
                    method is 'value', then this parameter is interpreted as the exact
                    read count below which cells must be removed.
        """
        if method not in list(_FILTER_METHODS_.values()):
            raise ValueError("filter_by_coverage: method must be one of %r." % list(_FILTER_METHODS_.values()))
        if method == _FILTER_METHODS_['PERC']:
            perc = np.quantile(list(self.cell_coverage.values()), float(threshold))
            filtered_cell_coverage = { cell: mad for cell, mad in self.cell_coverage.items() if mad <= perc}
        elif method == _FILTER_METHODS_['VALUE']:
            filtered_cell_coverage = { cell: mad for cell, mad in self.cell_coverage.items() if mad <= threshold}
        cells = filtered_cell_coverage.keys()
        filtered_cell_mad = { cell: self.cell_mad[cell] for cell in cells }
        if self.cell_ploidy != None:
            filtered_cell_ploidy = { cell: self.cell_ploidy[cell] for cell in cells }        
        else:
            filtered_cell_ploidy = None        
        cnvs_df = self.df.drop_cells(cells)
        sample = Sample(cnvs_df, filtered_cell_mad, filtered_cell_ploidy, filtered_cell_coverage, self.name)
        return sample

    
    
    


