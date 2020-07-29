from ..types import CnvData
import numpy as np
import pandas as pd
from typing import Union, List, Optional, Sequence
from  sklearn.decomposition import PCA
from sklearn.utils import check_array, check_random_state
from .. import logging as logg
from .._settings import settings
from ..utils import AnyRandom

def _pca(
    data: Union[CnvData, np.ndarray],
    n_comps: Optional[int] = None,
    svd_solver: str = 'arpack',
    random_state: AnyRandom = 0,
    use_highly_variable: Optional[bool] = None
) -> np.ndarray:
    if svd_solver in {'auto', 'randomized'}:
        logg.info(
            'Note that scikit-learn\'s randomized PCA might not be exactly '
            'reproducible across different computational platforms. For exact '
            'reproducibility, choose `svd_solver=\'arpack\'.`'
        )
    
    cnvdata = data if isinstance(data, CnvData) else CnvData(data)

    if use_highly_variable is True and 'highly_variable' not in cnvdata.feat.columns:
        raise ValueError(
            'Did not find cnv.feat[\'highly_variable\']. '
            'Either your data already only consists of highly-variable features '
            'or consider running `Sample.variable_features()` first.'
        )
    if use_highly_variable:
        logg.info('    on highly variable features')
    data_comp = (
        cnvdata[:, cnvdata.feat['highly_variable']] if use_highly_variable else cnvdata
    )

    if n_comps is None:
        min_dim = min(data_comp.n_feats, data_comp.n_obs)
        if settings.N_PCS >= min_dim:
            n_comps = min_dim - 1
        else:
            n_comps = settings.N_PCS

    logg.info(f'    with n_comps={n_comps}')

    random_state = check_random_state(random_state)

    X = data_comp.X
    pca_ = PCA(n_components=n_comps, svd_solver=svd_solver, random_state=random_state)
    X_pca = pca_.fit_transform(X)
    return (X_pca, pca_.components_, pca_.explained_variance_ratio_, pca_.explained_variance_)
