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
# single_sample_post_analysis.py: Single-sample analysis module
# ==========================================================================


from .funcs import *
from .check_funcs import *

#from dicttoxml import dicttoxml
#from xml.dom.minidom import parseString
#import json
import random
import pandas as pd
import numpy as np

class SingleSampleAnalysis:

    def __init__(self, sample: str = 'sample',  method: str ='ward', 
            metric:str='euclidean', verbose:bool =False, seed:int=None):
        self.sample = sample
        self.method = method
        self.metric = metric
        self.verbose = verbose
        self.seed = seed
    
        if self.seed != None:
            random.seed(int(seed))

    def heatmap_(self, cnvsf: str, results: str, outpath:str=None):
        df = pd.read_csv(cnvsf, sep="\t", usecols = lambda column : column not in ['Unnamed: 103', 'Unnamed: 113', 'Unnamed: 32'])

        cnvs = df.drop(['CHR', 'START', 'END'], axis=1).transpose()
        boundaries = df[['CHR', 'START', 'END']].copy()
        heatmap(cnvs, boundaries, self.method, self.metric, False, outpath=outpath, sample=self.sample)
    
    def clusters_(self):
        sys.exit()

    def het_score(self):
        sys.exit
        