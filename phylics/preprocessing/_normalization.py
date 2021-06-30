import warnings
from typing import Union, Optional
import numpy as np
import pandas as pd
from ..types import CnvData
from .. import logging as logg
from .._settings import settings
from ._utils import _get_mean_var


def scale(
    data: Union[np.ndarray, pd.DataFrame],
    zero_center: bool = True,
    max_value: Optional[float] = None
):
    """\
    Scale data to unit variance and zero mean (z-score)
    .. note::
        Variables (bins) that do not display any variation (are constant across
        all observations) are retained and (for zero_center==True) set to 0
        during this operation.
    Parameters
    ----------
    X
        The data matrix of shape `n_cells` × `n_bins`.
        Rows correspond to cells and columns to bins.
    zero_center
        If `False`, omit zero-centering variables
    max_value
        Clip (truncate) to this value after scaling. If `None`, do not clip.
    copy
        Whether this function should be performed inplace. If an AnnData object
        is passed, this also determines if a copy is returned.

    Returns
    -------
    Depending on `copy` returns or updates `cnvdata` with a scaled `cnvdata.X`,
    annotated with `'mean'` and `'std'` in `cnvdata.feat`.
    """
    X = data if isinstance(data, np.ndarray) else data.to_numpy() 
    X = X.copy()
    X = X.astype(float)

    mean, var = _get_mean_var(X)
    std = np.sqrt(var)
    std[std == 0] = 1
    if zero_center:
        X -= mean
    X /= std

    # do the clipping
    if max_value is not None:
        X[X > max_value] = max_value

    return X, mean, var
