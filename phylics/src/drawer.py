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
# drawer.py: this module implements all methods to draw figures
# ==========================================================================

__all__ = ['Drawer']


import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib 
import matplotlib.cm as cm
from matplotlib import rcParams
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.transforms as transforms
from typing import Tuple, Union, List, Iterable, Collection, Sequence, Optional
from .data_types import CnvData, VariableFeatures
from ._compat import Literal
from . import logging as logg

def heatmap(cnvs, boundaries, method='ward', metric='euclidean', outpath=None, verbose=False, sample=None,
                vmin:int = 0, vmax:int = 12, vcenter:int=2, figsize:Tuple[int, int]=(37, 21), fontsize:int=16): 
    divnorm = colors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    
    if method == 'ward' and metric != 'euclidean':
        raise ValueError("metric must be 'euclidean' if method is 'ward'")


    chr_limits = []
    for k,v in boundaries.groupby('CHR', sort=False)['END'].max().items():
        chr_limits.append(boundaries.index[(boundaries['CHR'] == k) & (boundaries['END'] == v)].tolist()[0])
    chr_boundaries = np.append(0, chr_limits)
    chr_list = boundaries['CHR'].unique().tolist()
    chrN_list = []

    for x in chr_list:
        x = x[3:] #remove 'chr' for readability
        chrN_list.append(x)

    #compute the position where chromosome labels will be placed on the plots
    start = 0
    pos_list = []
    for end in chr_limits:
        pos_list.append((start+end)/2)
        start = end+1

    
    cbar_kws={"ticks":np.arange(vmin,vmax+1,1)}
    h = sns.clustermap(cnvs, method=method, metric=metric, col_cluster=False, yticklabels = False,  
                cmap='RdBu_r', vmin=vmin, vmax=vmax,norm=divnorm, figsize=figsize, cbar_kws=cbar_kws)    
    h.cax.set_position([0.05, .2, .03, .45])

    ax = h.ax_heatmap
    #place vertical lines to identify chromosomes
    for pos in chr_limits:
        ax.axvline(x=pos, color='black')

    #place chromosome ticks at the right position
    ax.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter((chrN_list)))
    ax.tick_params(axis='x', rotation=0, labelsize=fontsize)
    
    ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_boundaries))
    ax.tick_params(axis='x', length=20, which='minor')

    ax.set_xlabel("Chromosomes", fontweight='bold', fontsize=fontsize)
    ax.set_ylabel("cells", fontsize=16, fontweight='bold')

    #plt.gcf().set_size_inches(37, 21)
    plt.gcf().suptitle("CNV heatmap - Sample = " + sample, fontsize=fontsize+2, fontweight='bold')

    if outpath != None:
        plt.savefig(outpath)
    else:
        plt.show()
    plt.clf()
    return h

def dist(a:Union[list, np.array, pd.Series], grid:bool=False, quantiles:List[float]=None, figsize:Tuple[int, int]=None, outpath:str=None, **kwargs):

    ax = sns.distplot(a, **kwargs)
    
    if quantiles != None:
        for q in quantiles:
            pos = np.quantile(a, q)
            plt.axvline(pos, color='red', linestyle='--', label=q)    
        plt.legend()
    if grid == True:
        plt.grid()
    if figsize != None:
        plt.gcf().set_size_inches(figsize)
    plt.subplots_adjust(wspace=0.3)
    if outpath != None:
        plt.savefig(outpath)
    else:
        plt.show()
    plt.clf()
    return ax

def scatter(x:np.array=None, y:np.array=None, 
                     X:np.ndarray=None, outpath:str=None, figsize:Tuple[int, int]=None, **kwargs):
    if X.any() != None:
        ax = sns.scatterplot(X[:,0], X[:,1], **kwargs)
    else: 
        ax = sns.scatterplot(x, y, **kwargs)
    if figsize != None:
        plt.gcf().set_size_inches(figsize)
    plt.subplots_adjust(wspace=0.3)
    if outpath != None:
        plt.savefig(outpath)
    else:
        plt.show()
    return ax

def highly_variable_features(X: VariableFeatures, log: bool = False, outpath: str = None):
    means = X.means
    dispersions = X.dispersions
    dispersions_norm = X.dispersions_norm
    feat_subset = X.highly_variable
    size = rcParams['figure.figsize']
    plt.figure(figsize=(2*size[0], size[1]))
    plt.subplots_adjust(wspace=0.3)
    for idx, d in enumerate([dispersions_norm, dispersions]):
        plt.subplot(1, 2, idx + 1)
        for label, color, mask in zip(['highly variable features', 'other features'],
                                      ['red', 'grey'],
                                      [feat_subset, ~feat_subset]):
            if False: means_, disps_ = np.log10(means[mask]), np.log10(d[mask])
            else: means_, disps_ = means[mask], d[mask]
            plt.scatter(means_, disps_, label=label, c=color, s=1)
        if log:  # there's a bug in autoscale
            plt.xscale('log')
            plt.yscale('log')
            min_dispersion = np.min(dispersions)
            y_min = 0.95*min_dispersion if min_dispersion > 0 else 1e-1
            plt.xlim(0.95*np.min(means), 1.05*np.max(means))
            plt.ylim(y_min, 1.05*np.max(dispersions))
        if idx == 0: plt.legend()
        plt.xlabel(('$log_{10}$ ' if False else '') + 'mean copy-number of bin')
        plt.ylabel(('$log_{10}$ ' if False else '') + 'dispersion of copy-numbers'
                  + (' (normalized)' if idx == 0 else ' (not normalized)'))
    if outpath != None:
        plt.savefig(outpath)
    else: 
        plt.show()
    ax = plt.gca()
    plt.clf()
    return ax




_DRAWING_FUNCTIONS_ = {
    'heatmap' : heatmap,
    'dist': dist,
    'scatter':scatter,
    'variable_features':highly_variable_features
}


class Drawer:
    @staticmethod
    def draw(plot_name:str, *args, **kwargs):
        func = _DRAWING_FUNCTIONS_[plot_name]
        return func(*args, **kwargs)