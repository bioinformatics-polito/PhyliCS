import pandas as pd
from numpy import random
from typing import Union, Mapping
import collections.abc as cabc 

__all__ = ["from_ginkgo_to_phylics", 
            "load_annotation_",
            "sanitize_annotation"]

AnyRandom = Union[None, int, random.Generator] 

def from_ginkgo_to_phylics(filepath:str): 
    df = pd.read_csv(filepath, sep="\t")
    X = df.drop(['CHR', 'START', 'END'], axis=1).transpose()
    boundaries = df[['CHR', 'START', 'END']]
    return X, boundaries

def load_annotation_(filepath:str):
    return pd.read_csv(filepath, header=None, squeeze=True, index_col=0, sep="\t")

def sanitize_annotation(ann:Union[pd.Series, Mapping[str, Union[float, int]]]):
    if isinstance(ann, pd.Series):
        return ann
    elif isinstance(ann, cabc.Mapping):
        return pd.Series(ann)
    else:
        raise ValueError("The annotation must be of type pd.Series or Mapping: got {})".format(type(ann)))



