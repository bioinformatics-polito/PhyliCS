#|/usr/bin/envs python

import argparse
import pandas as pd
import numpy as np
import umap
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description="umap")


parser.add_argument("input", metavar='SegCopy', action='store',
        help='cnvs file.', type=str)
parser.add_argument("outdir", metavar='outdir', action='store',
        help='Output directory path.', type=str)

args=parser.parse_args()

input_file = args.input
outdir = args.outdir

df = pd.read_csv(input_file, sep="\t", usecols = lambda column : column not in ['CHR', 'START', 'END']).transpose()
                        
standard_embedding = pd.DataFrame(data=umap.UMAP(random_state=42).fit_transform(df.values), index=df.index)

standard_embedding.to_csv(outdir+"/umap_reduction.csv", sep="\t")

