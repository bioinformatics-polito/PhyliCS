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
        
        if coverages != None:
            if isinstance(coverages, dict):
                self.coverages = coverages
            elif isinstance(ploidies, str):
                self.coverages = pd.read_csv(coverages, header=0, squeeze=True, index_col=0, sep="\t").to_dict()
        
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

    def filter_by_ploidy(self, ploidy:Union[float, Tuple[float, float]]):
        #filter out cells with/within the given ploidy/ploidies
    
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
    def filter_by_coverage(self, method:str='percentile', threshold:float='0.9')
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

    def highly_variant_features(self):
        return x
    
    def explained_variance_by_pcs(self):
        return x
    


