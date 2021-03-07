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
# custom_types.py: this module implements phylics worker types
# ==========================================================================
from ..utils import clustering_func
__all__ = ["ClusterConfig"]

class ClusterConfig:
    """
        Utility class which holds the clustering schema to be used with a given algorithm.        
    """
    def __init__(self, method, **kwargs):
        self.method_ = method
        for attr in kwargs.keys():
            self.__dict__[attr] = kwargs[attr]
            self.clusterer_ = clustering_func[method](**kwargs)
    def get_method(self):
        return self.method_
    def get_attr_dict(self):
        return self.__dict__
    def get_clusterer(self):
        return self.clusterer_
    