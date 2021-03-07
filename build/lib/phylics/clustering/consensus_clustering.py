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

import numpy as np
from itertools import combinations
import bisect

#TODO
# - fit() has to compute the consensus matrix for all configurations and return the best one

class ConsensusCluster:
    """
      Implementation of Consensus clustering, following the paper
      https://link.springer.com/content/pdf/10.1023%2FA%3A1023949509487.pdf
      Args:
        * cluster -> clustering class
        * NOTE: the class is to be instantiated with parameter `n_clusters`,
          and possess a `fit_predict` method, which is invoked on data.
        * H -> number of resamplings for each cluster config
        * resample_proportion -> percentage to sample
        * configurations -> list of instances of the class ClusterConfig to be tested
        * M -> consensus matrices for each ClusterConfig (shape =(n_configs,data.shape[0],data.shape[0]))
                (NOTE: every consensus matrix is retained, like specified in the paper)
        * A -> area under CDF for each ClusterConfig 
                (see paper: section 3.3.1. Consensus distribution.)
        * delta -> changes in areas under CDF
                (see paper: section 3.3.1. Consensus distribution.)
        * self.best -> best clustering configuration
      """

    def __init__(self, cluster, configurations, H, resample_proportion=0.5):
        assert 0 <= resample_proportion <= 1, "proportion has to be between 0 and 1"
        self.cluster_ = cluster
        self.resample_proportion_ = resample_proportion
        self.configurations = configurations
        self.H_ = H
        self.M = None
        self.A = None
        self.delta = None
        self.best = None
        

    def _internal_resample(self, data, proportion):
        """
        Args:
          * data -> (examples,attributes) format
          * proportion -> percentage to sample
        """
        resampled_indices = np.random.choice(
            range(data.shape[0]), size=int(data.shape[0]*proportion), replace=False)
        return resampled_indices, data[resampled_indices, :]

    def fit(self, data, verbose=False):
        """
        Fits a consensus matrix for each number of clusters
        Args:
          * data -> (examples,attributes) format
          * verbose -> should print or not
        """
        # pylint: disable=unsubscriptable-object
        M = np.zeros((len(self.configurations), data.shape[0], data.shape[0]))
        Is = np.zeros((data.shape[0],)*2)
        conf_lut = {k: v for v, k in enumerate(self.configurations)}
        for i_, conf in enumerate(self.configurations):  # for each number of clusters
            if verbose:
                print("Iteration = %d" % (i_))
            for h in range(self.H_):  # resample H times
                if verbose:
                    print("\tAt resampling h = %d" % (h))
                resampled_indices, resample_data = self._internal_resample(
                    data, self.resample_proportion_)
                Mh = self.cluster_(**conf.get_attr_dict()).fit_predict(resample_data)
                # find indexes of elements from same clusters with bisection
                # on sorted array => this is more efficient than brute force search
                id_clusts = np.argsort(Mh)
                sorted_ = Mh[id_clusts]
                for i in np.unique(sorted_):  # for each cluster
                    ia = bisect.bisect_left(sorted_, i)
                    ib = bisect.bisect_right(sorted_, i)
                    is_ = id_clusts[ia:ib]
                    ids_ = np.array(list(combinations(is_, 2))).T
                    # sometimes only one element is in a cluster (no combinations)
                    if ids_.size != 0:     
                        M[i_, ids_[0], ids_[1]] += 1 
                # increment counts
                ids_2 = np.array(list(combinations(resampled_indices, 2))).T
                Is[ids_2[0], ids_2[1]] += 1
            M[i_] /= Is+1e-8  # consensus matrix
            # Mk[i_] is upper triangular (with zeros on diagonal), we now make it symmetric
            M[i_] += M[i_].T
            M[i_, range(data.shape[0]), range(
                data.shape[0])] = 1  # always with self
            Is.fill(0)  # reset counter
        self.M = M
        # fits areas under the CDFs
        self.A = np.zeros(len(self.configurations))
        for i, m in enumerate(M):
            hist, bins = np.histogram(m.ravel(), density=True)
            self.A[i] = np.sum(h*(b-a)
                             for b, a, h in zip(bins[1:], bins[:-1], np.cumsum(hist)))
        # fits differences between areas under CDFs
        self.delta = np.array([(Ab-Aa)/Aa if i > 2 else Aa
                                for Ab, Aa, i in zip(self.A[1:], self.A[:-1], np.arange(len(self.configurations)))])
        self.best = np.argmax(self.delta) 

    