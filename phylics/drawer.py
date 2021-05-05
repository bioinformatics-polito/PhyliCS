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
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
from copy import copy
import matplotlib 
import matplotlib.cm as cm
from matplotlib.cm import get_cmap
from matplotlib import rcParams
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.transforms as transforms
from matplotlib.patches import Patch
from typing import Tuple, Union, List, Iterable, Collection, Sequence, Optional
from .types import CnvData
from ._compat import Literal
from . import logging as logg
from .plotting._utils import make_projection_available, ColorLike, _Dimensions

#matplotlib.rcParams.update({'font.size': 24})

def dots_plot(data:Union[np.array, list, pd.Series], yticks:Union[np.array, list, pd.Series]=None,  
               title:str=None, x_label:str="X", y_label:str="Y", figsize:tuple=(18, 7), outpath:str=None):
    fig2, ax1 = plt.subplots()
    fig2.set_size_inches(figsize)
            
    #yticks = yticks if yticks is not None else np.range(len(data))

    matplotlib.rcParams.update({'font.size': 16})
    ax1.scatter(yticks, data)

    for i, d in enumerate(data):
        txt =  "{:.4f}".format(d)
        ax1.annotate(txt, xy= (i+0.003, d+0.003))
   
    ax1.grid(True)
    plt.xticks(rotation=30)
    plt.subplots_adjust(hspace=0, bottom=0.3)
    ax1.set_xlabel(x_label, fontsize=16)
    ax1.set_ylabel(y_label, fontsize=16)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    if title is not None:
        fig2.suptitle(title, fontweight='bold')        
   
    if outpath != None:
        fig2.savefig(outpath)
    else:
        fig2.show()
    plt.close('all')

    return ax1


def clustermap2(data, boundaries, labels:Union[pd.Series, pd.DataFrame]=None,  row_cluster:bool=False, vmin:int = 0, vmax:int = 12, vcenter:int=2, 
            linkage:Union[np.ndarray, None]=None, outpath=None, legend: bool=False,  title: str=None, 
                 figsize:Tuple[int, int]=(37, 21)):
    warnings.filterwarnings("ignore")
    sns.set(font_scale=3.5)
    
    if figsize == None:
        figsize = (37, 21)
    divnorm = colors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    cbar_kws={"ticks":np.arange(vmin,vmax+1,1)}

    cluster_colors = None
    cluster_lut = {}
    if labels is not None:
        if isinstance(labels, pd.Series):
            color_palette = sns.color_palette("hls", len(np.unique(labels)))
            cluster_lut = dict(zip(np.unique(labels), color_palette))
            cluster_colors=pd.Series(labels, index=data.index, name=labels.name).map(cluster_lut)
        elif isinstance(labels, pd.DataFrame):
            cluster_colors = pd.DataFrame(columns=labels.columns)
            for c in labels.columns:
                color_palette = sns.color_palette("hls", len(np.unique(labels[c])))
                cluster_lut = dict(zip(np.unique(labels[c]), color_palette))
                cluster_colors[c] = pd.Series(labels[c], index=data.index, name=c).map(cluster_lut)
                
    h = sns.clustermap(data, row_cluster=row_cluster, col_cluster=False, row_linkage=linkage, yticklabels = False,  
                row_colors=cluster_colors, cmap='RdBu_r', vmin=vmin, vmax=vmax,norm=divnorm, cbar_kws=cbar_kws, 
                figsize=figsize
                )   
    
    if legend == True: 
        # Draw the legend bar for the classes                 
        # label in np.unique(labels):
        #    h.ax_row_dendrogram.bar(0, 0, color=cluster_lut[label],
        #                label=label, linewidth=0)
        #h.ax_row_dendrogram.legend(loc='best', bbox_to_anchor=(0.05, .2, .03, .45), ncol=len(np.unique(labels)), 
                #title_fontsize='x-large', bbox_transform=plt.gcf().transFigure)
        #plt.figlegend(title=labels.name, loc='center', bbox_to_anchor=(0.5,0.9))      
        handles = [Patch(facecolor=cluster_lut[l]) for l in cluster_lut]
        plt.legend(handles, cluster_lut, 
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')

        #h.ax_heatmap.legend(handles, cluster_lut, title=labels.name, bbox_to_anchor=(0.5, 0.9), ncol=len(np.unique(labels)), bbox_transform=plt.gcf().transFigure, loc='upper right')
    

    h.fig.tight_layout()
    h.fig.subplots_adjust(right=0.9, left=0.0)
    h.ax_cbar.set_position((0.95, .2, .03, .4))

    ax = h.ax_heatmap

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

    #place vertical lines to identify chromosomes
    for pos in chr_limits:
        ax.axvline(x=pos, color='black')

    #place chromosome ticks at the right position
    ax.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter((chrN_list)))
    ax.tick_params(axis='x',  width=2, rotation=90)
    
    ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_boundaries))
    ax.tick_params(axis='x', length=20, which='minor')

    ax.set_xlabel("Chromosomes",)
    ax.set_ylabel("Cells")

    #hm = h.ax_heatmap.get_position()
    #plt.setp(h.ax_heatmap.xaxis.get_majorticklabels(), fontsize=6)
    #h.ax_heatmap.set_position([hm.x0, hm.y0, hm.width*0.25, hm.height])
    #row = h.ax_row_dendrogram.get_position()
    #h.ax_row_dendrogram.set_position([row.x0, row.y0, row.width*0.25, row.height*0.5])


    #plt.gcf().set_size_inches(37, 21)
    if title is not None:
        plt.gcf().suptitle(title, fontsize=40)

    if outpath != None:
        plt.savefig(outpath)
    else:
        plt.show()
    plt.close('all')
    sns.set(font_scale=1)
    return h

def clustermap(data, boundaries, labels:Union[np.array, None]=None, outpath=None, title: str=None, legend: bool=False, linkage:Union[np.ndarray, None]=None,
                 vmin:int = 0, vmax:int = 12, vcenter:int=2, figsize:Tuple[int, int]=(37, 21), fontsize:int=16, row_cluster:bool=False): 
    warnings.filterwarnings("ignore")
    divnorm = colors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    cbar_kws={"ticks":np.arange(vmin,vmax+1,1)}
    if labels is not None:
        if linkage is None:
            #sort rows by cluster label
            data['label'] = labels
            data = data.sort_values(by='label')
            labels = data["label"].values
            data = data.drop(['label'], axis=1)
        
        # Cluster colors
        color_palette = sns.color_palette("hls", len(np.unique(labels)))
        cluster_lut = dict(zip(np.unique(labels), color_palette))
        cluster_colors=pd.Series(labels, index=data.index).map(cluster_lut)
        if linkage is None:
            if row_cluster == False:
                h = sns.clustermap(data, row_cluster=False, col_cluster=False, yticklabels = False,  row_colors=cluster_colors,
                    cmap='RdBu_r', vmin=vmin, vmax=vmax,norm=divnorm, figsize=figsize, cbar_kws=cbar_kws)   
            else:
                 h = sns.clustermap(data, method="ward", col_cluster=False, yticklabels = False,  row_colors=cluster_colors,
                    cmap='RdBu_r', vmin=vmin, vmax=vmax,norm=divnorm, figsize=figsize, cbar_kws=cbar_kws)
        else:
            h = sns.clustermap(data, row_linkage=linkage, col_cluster=False, yticklabels = False,  row_colors=cluster_colors,
                    cmap='RdBu_r', vmin=vmin, vmax=vmax, norm=divnorm, figsize=figsize, cbar_kws=cbar_kws)

    else:
        if row_cluster == False:
            h = sns.clustermap(data, method="ward", row_cluster=False, col_cluster=False, yticklabels = False,  
                    cmap='RdBu_r', vmin=vmin, vmax=vmax,norm=divnorm, figsize=figsize, cbar_kws=cbar_kws)
        else:
            h = sns.clustermap(data, method="ward", col_cluster=False, yticklabels = False,  
                    cmap='RdBu_r', vmin=vmin, vmax=vmax,norm=divnorm, figsize=figsize, cbar_kws=cbar_kws)
    
    if legend == True: 
        # Draw the legend bar for the classes                 
        for label in np.unique(labels):
            h.ax_row_dendrogram.bar(0, 0, color=cluster_lut[label],
                        label=label, linewidth=0)
        h.ax_row_dendrogram.legend(loc='center', bbox_to_anchor=(0.5,0.9), ncol=len(np.unique(labels)), fontsize=14, 
                title_fontsize='x-large', bbox_transform=plt.gcf().transFigure)
        
    h.cax.set_position([0.05, .2, .03, .45])
    ax = h.ax_heatmap

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
    if title is not None:
        plt.gcf().suptitle(title)

    if outpath != None:
        plt.savefig(outpath)
    else:
        plt.show()
    plt.close('all')
    return h

def dist(a:Union[list, np.array, pd.Series], grid:bool=False, quantiles:List[float]=None, figsize:Tuple[int, int]=None, outpath:str=None, **kwargs):

    ax = sns.displot(a, **kwargs)
    
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
    plt.close('all')
    return ax

_ScatterBasis = Literal['pca', 'umap', 'X']


def scatter(data:np.ndarray, projection: _Dimensions = "2d", outpath:str=None, title: str=None, 
                x_label: str=None, y_label: str=None, z_label: str=None,
                alpha: float = 1.0, bw: bool = False, cmap: str = "Paired", 
                na_color: ColorLike = "lightgray", labels: Union[np.array, pd.Series, None] = None, legend: bool=False,
                wspace : Union[float, None]= 0.3, figsize:Tuple[int, int]=None, show: bool=False, **kwargs):
    
    warnings.filterwarnings('error')
    n_dimensions = data.shape[1]
    if n_dimensions < 2:
        logg.error(
                        f"TypeError: Input data format is not valid. At least a 2D array is required. "
                    )
        ax =  None

    else:
        if projection == "2d":
            X = data[:, [0,1]] 
        elif projection == "3d":
            if n_dimensions < 3:
                logg.warning(
                        f"TypeError: Input data has only 2 projection. "
                    )
                X = data[:, [0,1]]
            else:
                X = data[:, [0,1,2]]
        else:
            logg.warning(
                        f"InvalidArgument: projection={projection} is not valid. "
                    )
            projection = "2d"
            X = data[:, [0,1]]    
 
        fig = plt.figure()
        if projection == "3d":
            ax = fig.add_subplot(111, projection=projection)
        else:
            ax = fig.add_subplot(111)

        if labels is not None:
            if isinstance(labels, pd.Series):
                labels = labels.values
            cmap = copy(get_cmap(cmap, len(np.unique(labels))))
            cmap.set_bad(na_color)
            na_color = colors.to_hex(na_color, keep_alpha=True)
            print(na_color)
            clusters_dict = {k: v for v, k in enumerate(np.unique(labels))}
            cluster_colors = [cmap.colors[clusters_dict[x]] if x >= 0
                      else na_color
                      for x in labels] 
            colors_legend = [cmap.colors[clusters_dict[x]] if x >= 0
                      else na_color
                      for x in np.unique(labels)]
            kwargs["c"] = cluster_colors
            if labels is not None and legend==True:
                #ax.legend()
                markers = [plt.Line2D([0,0],[0,0],color=color, marker='o', linestyle='') for color in colors_legend ]
                fig.legend(markers, clusters_dict.keys(), numpoints=1)
        else:
            cmap = copy(get_cmap(cmap))
            cmap.set_bad(na_color)
            na_color = colors.to_hex(na_color, keep_alpha=True)
            kwargs["cmap"] = cmap
    
        
        if title is None:
            ax.set_title('')
        else:
            ax.set_title(title, fontsize=18)

        if projection == "2d":
            cax = ax.scatter(X[:,0], 
                    X[:,1], 
                    alpha=alpha,
                    marker=".",
                    **kwargs
                    )
        else:
            cax = ax.scatter(X[:,0], 
                    X[:,1], 
                    X[:,2],
                    alpha=alpha,
                    marker=".",
                    **kwargs
                    )
        if x_label is not None:
            ax.set_xlabel(x_label, fontsize=16)
        if y_label is not None:
            ax.set_ylabel(y_label, fontsize=16)
        if projection == '3d' and z_label is not None:
            ax.set_zlabel(z_label, labelpad=-7)
        ax.autoscale_view()

        if figsize is None:
            figsize = rcParams['figure.figsize']
        if wspace is None:
            #  try to set a wspace that is not too large or too small given the
            #  current figure size
            wspace = 0.75 / rcParams['figure.figsize'][0] + 0.0
    
        fig.set_size_inches(figsize)
        fig.subplots_adjust(wspace=wspace)
        if outpath != None:
            fig.savefig(outpath)
        else:
            plt.show()
        plt.close('all')
    
    return ax 

def highly_variable_features(X: CnvData, log: bool = False, outpath: str = None):
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
    plt.close('all')
    return ax

def jackstraw(explained_variance_ratio: np.array, mean_expl_var_ratio_perm: np.array, n_pcs:int = 50, figsize:Tuple[int, int]=None,
                show: bool=False, outpath: str = None):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(explained_variance_ratio, color='green', label='Explained by PCs')
    ax.plot(mean_expl_var_ratio_perm, color='red', label='Explained by chance')

    ax.set_xlabel('number of components')
    ax.set_ylabel('cumulative explained variance')

    plt.xticks(np.arange(0, n_pcs, 1))

    ax.legend()
    ax.grid(True)

    fig.suptitle("Explained variance ratio")

    if figsize is None:
        figsize = rcParams['figure.figsize']
    fig.set_size_inches(figsize)
    wspace = 0.75 / rcParams['figure.figsize'][0] + 0.0
    fig.subplots_adjust(wspace=wspace)
    if outpath != None:
        fig.savefig(outpath)
    else:
        if show == True:
            fig.show()
    plt.close('all')

    return ax


_DRAWING_FUNCTIONS_ = {
    'heatmap' : clustermap2,
    'dist': dist,
    'scatter':scatter,
    'variable_features':highly_variable_features,
    'jackstraw': jackstraw,
    'dots': dots_plot
}


class Drawer:
    @staticmethod
    def draw(plot_name:str, *args, **kwargs):
        func = _DRAWING_FUNCTIONS_[plot_name]
        return func(*args, **kwargs)