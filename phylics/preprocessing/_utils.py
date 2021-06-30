import numpy as np
from pandas import DataFrame
from scipy import sparse
import numba
from typing import Union

def _get_mean_var(X: Union[DataFrame, np.ndarray], *, axis=0):
    if isinstance(X, DataFrame):
        mean, var = _mean_variance_axis(X, axis=axis)
    else:
        mean = np.mean(X, axis=axis, dtype=np.float64)
        mean_sq = np.multiply(X, X).mean(axis=axis, dtype=np.float64)
        var = mean_sq - mean ** 2
    # enforce R convention (unbiased estimator) for variance
    var *= X.shape[axis] / (X.shape[axis] - 1)
    return mean, var

def _mean_variance_axis(mtx: DataFrame, axis: int):
    assert axis in (0, 1)
    if isinstance(mtx, DataFrame) == False:
        raise ValueError("This function only works on pandas DataFrames")
    var = mtx.var(axis = axis)
    mean = mtx.mean(axis = axis)
    return mean, var

