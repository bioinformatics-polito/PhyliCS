#!/usr/bin/env python

# ==========================================================================
#                                  PhyliCS
# ==========================================================================
# This file is part of PhyliCS.
#
# TOOL is Free Software: you can redistribute it and/or modify it
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
# valid_cells.py: Cell filtering module
# ==========================================================================


import os
import argparse
import numpy as np
import pandas as pd
from check_funcs import *
from funcs import print_msg


parser = argparse.ArgumentParser(description="Ploidy interval filtering.")

parser.add_argument("sample", metavar='sample_name', action='store',
        help='Sample name',
        type=str)


parser.add_argument("results", metavar='results.txt', action='store',
                help='Path to stats file',
                        type=str)

parser.add_argument("segcopy", metavar='SegCopy', action='store',
                help='Path to cnvs file',
                        type=str)


parser.add_argument("outliers", metavar='interval/value', action=check_valid_outliers(),
        help='Interval of values (p1-p2) or single values (p) to be filtered out. At least one interval/value must be specified.\nIntervals must be compliant to the following format: \n\tp1-p2.\nValues must be compliant to the following format:\n\tp.',
        nargs="+", type=str)

parser.add_argument("outdir", metavar='out_dir', action='store',
                help='Path to the output directory', type=str)
parser.add_argument("--verbose", action='store_true', help='Verbose execution.')


args=parser.parse_args()

results = args.results
segcopy = args.segcopy
sample = args.sample
outliers = args.outliers
outdir = args.outdir
verbose = False

if args.verbose:
   verbose = True

df = pd.read_csv(results, sep="\t").dropna()

tot_cells = len(df)

for outlier in outliers:
    if is_interval(outlier):
        ploidy = outlier.split('-')
        low = float(ploidy[0])
        high = float(ploidy[1])
        df = df.loc[(df["Copy_Number"] < low) | (df["Copy_Number"] > high)]
    elif is_number(outlier):
            p = float(outlier)
            df = df.loc[(df["Copy_Number"] != p)]


df.to_csv(os.path.join(outdir, "results.txt"), sep='\t', index=False)

valid_cells = df.Sample.tolist()

remaining_cells = len(valid_cells)
removed_cells = tot_cells - remaining_cells


header_cols = ['CHR', 'START', 'END']
valid_cols = np.append(valid_cells, header_cols)

df = pd.read_csv(segcopy, sep="\t", usecols = valid_cols)
df.to_csv(os.path.join(outdir,'SegCopy'), sep='\t', index=False)

print_msg("Initial cells: {}".format(tot_cells), 0, verbose)
print_msg("Filtered out cells: {}".format(removed_cells), 0, verbose)
print_msg("Remaining cells: {}".format(remaining_cells), 0, verbose)






