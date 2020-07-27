import warnings
from typing import Optional

import numpy as np
import pandas as pd
from ..data_types import CnvData

from ._utils import _get_mean_var
from ._distributed import materialize_as_ndarray


def _highly_variable_features(
    cnvs: CnvData,
    min_disp: Optional[float] = None,
    max_disp: Optional[float] = None,
    min_mean: Optional[float] = None,
    max_mean: Optional[float] = None,
    n_top_features: Optional[int] = None,
    n_bins: int = 20,
) -> pd.DataFrame:
    """\
    Returns
    -------
    A DataFrame that contains the columns
    `highly_variable`, `means`, `dispersions`, and `dispersions_norm`.
    """

    if min_disp is None: min_disp = 0.5
    if min_mean is None: min_mean = 0.0125
    if max_mean is None: max_mean = 3
    if max_disp is None: max_disp = np.inf

    X = cnvs.get_cnvs()
    '''
    if flavor == 'seurat':
        if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
            X *= np.log(adata.uns['log1p']['base'])
        X = np.expm1(X)
    '''
    mean, var = materialize_as_ndarray(_get_mean_var(X))
    # now actually compute the dispersion
    mean[mean == 0] = 1e-12  # set entries equal to zero to small value
    dispersion = var / mean
    """
    if flavor == 'seurat':  # logarithmized mean as in Seurat
        dispersion[dispersion == 0] = np.nan
        dispersion = np.log(dispersion)
        mean = np.log1p(mean)
    """
    # all of the following quantities are "per-gene" here
    df = pd.DataFrame()
    df['means'] = mean
    df['dispersions'] = dispersion
    """
    if flavor == 'seurat':
        df['mean_bin'] = pd.cut(df['means'], bins=n_bins)
        disp_grouped = df.groupby('mean_bin')['dispersions']
        disp_mean_bin = disp_grouped.mean()
        disp_std_bin = disp_grouped.std(ddof=1)
        # retrieve those genes that have nan std, these are the ones where
        # only a single gene fell in the bin and implicitly set them to have
        # a normalized disperion of 1
        one_gene_per_bin = disp_std_bin.isnull()
        gen_indices = np.where(one_gene_per_bin[df['mean_bin'].values])[0].tolist()
        if len(gen_indices) > 0:
            logg.debug(
                f'Gene indices {gen_indices} fell into a single bin: their '
                'normalized dispersion was set to 1.\n    '
                'Decreasing `n_bins` will likely avoid this effect.'
            )
        # Circumvent pandas 0.23 bug. Both sides of the assignment have dtype==float32,
        # but there’s still a dtype error without “.value”.
        disp_std_bin[one_gene_per_bin.values] = disp_mean_bin[one_gene_per_bin.values].values
        disp_mean_bin[one_gene_per_bin.values] = 0
        # actually do the normalization
        df['dispersions_norm'] = (
            (
                df['dispersions'].values  # use values here as index differs
                - disp_mean_bin[df['mean_bin'].values].values
            ) / disp_std_bin[df['mean_bin'].values].values
        )
    
    elif flavor == 'cell_ranger':
    """
    from statsmodels import robust
    df['mean_bin'] = pd.cut(df['means'], np.r_[
        -np.inf,
        np.percentile(df['means'], np.arange(10, 105, 5)),
        np.inf
    ])
    disp_grouped = df.groupby('mean_bin')['dispersions']
    disp_median_bin = disp_grouped.median()
    # the next line raises the warning: "Mean of empty slice"
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        disp_mad_bin = disp_grouped.apply(robust.mad)
        df['dispersions_norm'] = (df['dispersions'].values
        - disp_median_bin[df['mean_bin'].values].values
        ) / disp_mad_bin[df['mean_bin'].values].values
    
    dispersion_norm = df['dispersions_norm'].values.astype('float32')
    if n_top_features is not None:
        dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
        dispersion_norm[::-1].sort()  # interestingly, np.argpartition is slightly slower
        if n_top_features > len(cnvs.features):
            #logg.info(f'`n_top_genes` > `adata.n_var`, returning all genes.')
            n_top_features =len(cnvs.features)
        disp_cut_off = dispersion_norm[n_top_features-1]
        feat_subset = np.nan_to_num(df['dispersions_norm'].values) >= disp_cut_off
        """
        logg.debug(
            f'the {n_top_genes} top genes correspond to a '
            f'normalized dispersion cutoff of {disp_cut_off}'
        )
        """
    else:
        dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
        feat_subset = np.logical_and.reduce((
            mean > min_mean, mean < max_mean,
            dispersion_norm > min_disp,
            dispersion_norm < max_disp,
        ))

    df['highly_variable'] = feat_subset
    return df

