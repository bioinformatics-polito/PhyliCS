from typing import Optional, Union, Dict

import numpy as np
from anndata import AnnData
from sklearn.utils import check_random_state, check_array

from .. import logging as logg
from .._settings import settings
from .._compat import Literal
from ..utils import AnyRandom, _InitPos
from ..types import CnvData
from umap import UMAP


def _umap(
    data: Union[CnvData, np.ndarray],
    n_neighbors: int = 15, 
    n_components: int = 2, 
    metric: str = 'euclidean', 
    metric_kwds: Dict = None,
    min_dist: float = 0.5,
    spread: float = 1.0,
    maxiter: Optional[int] = None,
    alpha: float = 1.0,
    gamma: float = 1.0,
    fast: Union[bool, None] = False,
    negative_sample_rate: int = 5,
    local_connectivity: Union[int, None] = 1,
    init_pos: Union[_InitPos, np.ndarray, None] = 'spectral',
    random_state: AnyRandom = 0,
    a: Optional[float] = None,
    b: Optional[float] = None,
) -> np.ndarray:
    """\
    Embed the neighborhood graph using UMAP [McInnes18]_.
    UMAP (Uniform Manifold Approximation and Projection) is a manifold learning
    technique suitable for visualizing high-dimensional data. Besides tending to
    be faster than tSNE, it optimizes the embedding such that it best reflects
    the topology of the data, which we represent throughout Scanpy using a
    neighborhood graph.  We use the
    implementation of `umap-learn <https://github.com/lmcinnes/umap>`__
    [McInnes18]_.
    Parameters
    ----------
    data
        cnv data matrix.
    n_neighbors: float (optional, default 15)
 |      The size of local neighborhood (in terms of number of neighboring
 |      sample points) used for manifold approximation. Larger values
 |      result in more global views of the manifold, while smaller
 |      values result in more local data being preserved. In general
 |      values should be in the range 2 to 100.
 |  
 |  n_components: int (optional, default 2)
 |      The dimension of the space to embed into. This defaults to 2 to
 |      provide easy visualization, but can reasonably be set to any
 |      integer value in the range 2 to 100.
 |  
 |  metric: string or function (optional, default 'euclidean')
 |      The metric to use to compute distances in high dimensional space.
 |      If a string is passed it must match a valid predefined metric. If
 |      a general metric is required a function that takes two 1d arrays and
 |      returns a float can be provided. For performance purposes it is
 |      required that this be a numba jit'd function. Valid string metrics
 |      include:
 |          * euclidean
 |          * manhattan
 |          * chebyshev
 |          * minkowski
 |          * canberra
 |          * braycurtis
 |          * mahalanobis
 |          * wminkowski
 |          * seuclidean
 |          * cosine
 |          * correlation
 |          * haversine
 |          * hamming
 |          * jaccard
 |          * dice
 |          * russelrao
 |          * kulsinski
 |          * ll_dirichlet
 |          * hellinger
 |          * rogerstanimoto
 |          * sokalmichener
 |          * sokalsneath
 |          * yule
 |      Metrics that take arguments (such as minkowski, mahalanobis etc.)
 |      can have arguments passed via the metric_kwds dictionary. At this
 |      time care must be taken and dictionary elements must be ordered
 |      appropriately; this will hopefully be fixed in the future.
    min_dist
        The effective minimum distance between embedded points. Smaller values
        will result in a more clustered/clumped embedding where nearby points on
        the manifold are drawn closer together, while larger values will result
        on a more even dispersal of points. The value should be set relative to
        the ``spread`` value, which determines the scale at which embedded
        points will be spread out. The default of in the `umap-learn` package is
        0.1.
    spread
        The effective scale of embedded points. In combination with `min_dist`
        this determines how clustered/clumped the embedded points are.
    maxiter
        The number of iterations (epochs) of the optimization. Called `n_epochs`
        in the original UMAP.
    alpha
        The initial learning rate for the embedding optimization.
    gamma
        Weighting applied to negative samples in low dimensional embedding
        optimization. Values higher than one will result in greater weight
        being given to negative samples.
    fast: bool (optional, default False)
 |      For some datasets the nearest neighbor computation can consume a lot of
 |      memory. If you find that UMAP is failing due to memory constraints
 |      consider setting this option to True. This approach is more
 |      computationally expensive, but avoids excessive memory use.
    negative_sample_rate
        The number of negative edge/1-simplex samples to use per positive
        edge/1-simplex sample in optimizing the low dimensional embedding.
    local_connectivity: int (optional, default 1)
 |      The local connectivity required -- i.e. the number of nearest
 |      neighbors that should be assumed to be connected at a local level.
 |      The higher this value the more connected the manifold becomes
 |      locally. In practice this should be not more than the local intrinsic
 |      dimension of the manifold.
    init_pos
        How to initialize the low dimensional embedding. Called `init` in the
        original UMAP. Options are:
        * 'spectral': use a spectral embedding of the graph.
        * 'random': assign initial embedding positions at random.
        * A numpy array of initial embedding positions.
    random_state
        If `int`, `random_state` is the seed used by the random number generator;
        If `RandomState` or `Generator`, `random_state` is the random number generator;
        If `None`, the random number generator is the `RandomState` instance used
        by `np.random`.
    a
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    b
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    copy
        Return a copy instead of writing to adata.
    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.
    **X_umap** : `adata.obsm` field
        UMAP coordinates of data.
    """
    data = data if isinstance(data, CnvData) else CnvData(data)
    
    start = logg.info('computing UMAP') 
    from umap.umap_ import find_ab_params, simplicial_set_embedding
    if a is None or b is None:
        a, b = find_ab_params(spread, min_dist)
    else:
        a = a
        b = b
    random_state = check_random_state(random_state)
    
    # the data matrix X is really only used for determining the number of connected components
    # for the init condition in the UMAP embedding
    n_epochs = 0 if maxiter is None else maxiter
    umap_ = UMAP(
            n_neighbors = n_neighbors, 
            n_components = n_components, 
            metric=metric, 
            metric_kwds=metric_kwds, 
            n_epochs=n_epochs, 
            learning_rate=alpha, 
            init = init_pos, 
            min_dist = min_dist, 
            spread = spread, 
            low_memory = fast, 
            local_connectivity=local_connectivity, 
            repulsion_strength=gamma, 
            negative_sample_rate=negative_sample_rate, 
            a=a, 
            b=b, 
            random_state=random_state,  
            verbose=settings.verbosity > 3,
    )
    X_umap = umap_.fit_transform(data.X)
    logg.info(
        '    finished',
        time=start,
        deep=(
            'added\n'
            "    'X_umap', UMAP coordinates"
        ),
    )
    return X_umap