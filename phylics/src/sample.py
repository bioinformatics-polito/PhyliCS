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


from .custom_types import CNVS
from .drawer import Drawer
import random
import pandas as pd
import numpy as np
from typing import Union, Tuple

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

    def __init__(self, cnv_matrix:Union[CNVS, str], mads:Union[dict, str]=None, ploidies:Union[dict, str]=None, 
                        coverages:Union[dict, str]=None, sample_name:str="sample"):
        if isinstance(cnv_matrix, CNVS):
            self.df = cnv_matrix
        elif isinstance(cnv_matrix, str):
            self.df = CNVS(cnv_matrix)
        
        self.cnvs = self.df.get_cnvs()
        self.boundaries = self.df.get_boundaries()
        self.cells = self.df.get_cells()
        self.name = sample_name

        if ploidies != None:
            if isinstance(ploidies, dict):
                self.ploidies = ploidies
            elif isinstance(ploidies, str):
                self.ploidies = pd.read_csv(ploidies, header=0, squeeze=True, index_col=0, sep="\t").to_dict()
        else:
            self.ploidies = None
        
        if coverages != None:
            if isinstance(coverages, dict):
                self.coverages = coverages
            elif isinstance(ploidies, str):
                self.coverages = pd.read_csv(coverages, header=0, squeeze=True, index_col=0, sep="\t").to_dict()
        else:
            self.coverages = None
        
        if mads == None:
            self.mads = self.df.get_mads()
        elif isinstance(mads, dict):
            self.mads = mads
        elif isinstance(mads, str):
            self.mads = pd.read_csv(mads, header=0, squeeze=True, index_col=0, sep="\t").to_dict()
    
    def set_mads(self, mads:Union[dict, str]):
        if isinstance(mads, dict):
            self.mads = mads
        elif isinstance(mads, str):
            self.mads = pd.read_csv(mads, header=0, squeeze=True, index_col=0, sep="\t").to_dict()

    def set_ploidies(self, ploidies:Union[dict, str]):
        if isinstance(ploidies, dict):
            self.ploidies = ploidies
        elif isinstance(ploidies, str):
            self.ploidies = pd.read_csv(ploidies, header=0, squeeze=True, index_col=0, sep="\t").to_dict()

    def set_coverages(self, coverages:Union[dict, str]):
        if isinstance(coverages, dict):
            self.coverages = coverages
        elif isinstance(coverages, str):
            self.coverages = pd.read_csv(coverages, header=0, squeeze=True, index_col=0, sep="\t").to_dict()        
    
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

    def get_ploidies(self):
        return self.ploidies
    
    def get_coverages(self):
        return self.coverages
    
    def get_mads(self):
        return self.mads

    def count(self):
        return len(self.cells)

    def plot_ploidy_dist(self, kde:bool=True, rug:bool=False, vertical:bool=False, 
                axlabel:str=None, label:str=None, figsize:Tuple[int, int]=None, outpath:str=None):
        Drawer.draw('dist', list(self.ploidies.values()), kde, rug, vertical, axlabel, label, figsize, outpath)
    
    def plot_mad_dist(self, kde:bool=True, rug:bool=False, vertical:bool=False, 
                axlabel:str=None, label:str=None, figsize:Tuple[int, int]=None, outpath:str=None):
        Drawer.draw('dist', list(self.mads.values()), kde, rug, vertical, axlabel, label, figsize, outpath)

    def plot_cell_coverage_dist(self, kde:bool=True, rug:bool=False, vertical:bool=False, 
                axlabel:str=None, label:str=None, figsize:Tuple[int, int]=None, outpath:str=None):
        Drawer.draw('dist', list(self.coverages.values()), kde, rug, vertical, axlabel, label, figsize, outpath)

    def heatmap(self, method:str ='ward', metric:str ='euclidean', outpath:str=None,
                    vmin:int = 0, vmax:int = 12, vcenter:int=2, figsize:Tuple[int, int]=(37, 21), fontsize:int=16):
        Drawer.draw('heatmap', self.cnvs, self.boundaries, method, metric,  outpath=outpath, sample=self.name, 
            vmin = vmin, vmax=vmax, vcenter=vcenter, figsize=figsize, fontsize=fontsize)
    
    #Filter functions

    def filter_by_ploidy(self, method:str='eq', ploidy:Union[float, int, Tuple[float, float]]=2.0):
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
        ploidy: single value or interval of ploidies
        """
        if (method not in list(_FILTER_BY_SINGLE_PLOIDY_METHODS_.values())) and (method not in list(_FILTER_BY_RANGE_PLOIDY_METHODS_.values())):
            raise ValueError("filter_by_ploidy: method must be one of %r or %r."
            % list(_FILTER_BY_SINGLE_PLOIDY_METHODS_.values()), list(_FILTER_BY_RANGE_PLOIDY_METHODS_.values()))
        cells = []
        if isinstance(ploidy, float) or isinstance(ploidy, int):
            if method == _FILTER_BY_SINGLE_PLOIDY_METHODS_['EQ']:
                for cell, p in self.ploidies.items():
                    if p == ploidy:
                        cells.append(cell)
            elif method == _FILTER_BY_SINGLE_PLOIDY_METHODS_['LT']:
                for cell, p in self.ploidies.items():
                    if p < ploidy:
                        cells.append(cell)
            elif method == _FILTER_BY_SINGLE_PLOIDY_METHODS_['LT_EQ']:
                for cell, p in self.ploidies.items():
                    if p <= ploidy:
                        cells.append(cell)
            elif method == _FILTER_BY_SINGLE_PLOIDY_METHODS_['GT']:
                for cell, p in self.ploidies.items():
                    if p > ploidy:
                        cells.append(cell)
            elif method == _FILTER_BY_SINGLE_PLOIDY_METHODS_['GT_EQ']:
                for cell, p in self.ploidies.items():
                    if p >= ploidy:
                        cells.append(cell)
            else:
                raise ValueError("filter_by_ploidy: with a single value, method must be one of  %r." % list(_FILTER_BY_SINGLE_PLOIDY_METHODS_.values()))
        elif isinstance(ploidy, Tuple):
            if ploidy[0] >= ploidy[1]:
                raise ValueError("filter_by_ploidy: the first element of the ploidy tuple must be less than the second element.")
            if method == _FILTER_BY_RANGE_PLOIDY_METHODS_['IN']:
                for cell, p in self.ploidies.items():
                    if p > ploidy[0] and p < ploidy[1]:
                        cells.append(cell)
            elif method ==  _FILTER_BY_RANGE_PLOIDY_METHODS_['IN_EQ']:
                for cell, p in self.ploidies.items():
                    if p >= ploidy[0] and p <= ploidy[1]:
                        cells.append(cell)
            elif method ==  _FILTER_BY_RANGE_PLOIDY_METHODS_['OUT']:
                for cell, p in self.ploidies.items():
                    if p < ploidy[0] or p > ploidy[1]:
                        cells.append(cell)
            elif method ==  _FILTER_BY_RANGE_PLOIDY_METHODS_['OUT_EQ']:
                for cell, p in self.ploidies.items():
                    if p <= ploidy[0] or p >= ploidy[1]:
                        cells.append(cell)
            else:
                raise ValueError("filter_by_ploidy: with a range of values, method must be one of  %r." % list(_FILTER_BY_RANGE_PLOIDY_METHODS_.values()))
    
        cnvs_df = self.df.drop_cells(cells)
        filtered_mads = { cell: self.mads[cell] for cell in cells }
        if self.ploidies != None:
            filtered_ploidies = { cell: self.ploidies[cell] for cell in cells }
        else:
            filtered_ploidies = None
        if self.coverages != None:
            filtered_coverages = { cell: self.coverages[cell] for cell in cells } 
        else:
            filtered_coverages = None
        sample = Sample(cnvs_df, filtered_mads, filtered_ploidies, filtered_coverages, self.name)
            
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
            mads = pd.Series(self.mads)
            perc = mads.percentile(threshold)
            filtered_mads = mads[mads <= perc]
        elif method == _FILTER_METHODS_['VALUE']:
            mads = pd.Series(self.mads)
            filtered_mads = mads[mads <= threshold]
        cells = filtered_mads.index.values
        if self.ploidies != None:
            filtered_ploidies = { cell: self.ploidies[cell] for cell in cells }
        else:
            filtered_ploidies = None
        if self.coverages != None:
            filtered_coverages = { cell: self.coverages[cell] for cell in cells }
        else:
            filtered_coverages = None
        cnvs_df = self.df.drop_cells(cells)
        sample = Sample(cnvs_df, filtered_mads.to_dict(), filtered_ploidies, filtered_coverages, self.name)
        
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
            coverages = pd.Series(self.coverages)
            perc = coverages.percentile(threshold)
            filtered_coverages = coverages[coverages <= perc]
        elif method == _FILTER_METHODS_['VALUE']:
            coverages = pd.Series(self.coverages)
            filtered_coverages = coverages[coverages <= threshold]
        cells = filtered_coverages.index.values
        filtered_mads = { cell: self.mads[cell] for cell in cells }
        if self.ploidies != None:
            filtered_ploidies = { cell: self.ploidies[cell] for cell in cells }        
        else:
            filtered_ploidies = None        
        cnvs_df = self.df.drop_cells(cells)
        sample = Sample(cnvs_df, filtered_mads, filtered_ploidies, filtered_coverages.to_dict, self.name)
        return sample

    
    
    


