import warnings
import collections.abc as cabc
from abc import ABC
from functools import lru_cache
from typing import Union, List, Sequence, Tuple, Collection, Optional
from .._compat import Literal

import numpy as np
from matplotlib import pyplot as pl
from matplotlib import rcParams, ticker, gridspec, axes
from matplotlib.axes import Axes
from matplotlib.colors import is_color_like
from matplotlib.figure import SubplotParams as sppars, Figure
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from cycler import Cycler, cycler

ColorLike = Union[str, Tuple[float, ...]]
_IGraphLayout = Literal['fa', 'fr', 'rt', 'rt_circular', 'drl', 'eq_tree', ...]
_FontWeight = Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
_FontSize = Literal[
    'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'
]

_Dimensions = Literal["2d", "3d"]


def make_projection_available(projection):
    avail_projections = {'2d', '3d'}
    if projection not in avail_projections:
        raise ValueError(f'choose projection from {avail_projections}')
    if projection == '2d':
        return

    from io import BytesIO
    from matplotlib import __version__ as mpl_version
    from mpl_toolkits.mplot3d import Axes3D

    fig = Figure()
    ax = Axes3D(fig)

    circles = PatchCollection([Circle((5, 1)), Circle((2, 2))])
    ax.add_collection3d(circles, zs=[1, 2])

    buf = BytesIO()
    try:
        fig.savefig(buf)
    except ValueError as e:
        if not 'operands could not be broadcast together' in str(e):
            raise e
        raise ValueError(
            'There is a known error with matplotlib 3d plotting, '
            f'and your version ({mpl_version}) seems to be affected. '
            'Please install matplotlib==3.0.2 or wait for '
            'https://github.com/matplotlib/matplotlib/issues/14298'
        )

