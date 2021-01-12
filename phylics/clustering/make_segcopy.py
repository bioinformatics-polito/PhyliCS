import os
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "From intersectBed to ginkgo beds")

parser.add_argument("beds_dir", metavar="input/dir/path", action="store", help="Input directory path", type=str)
parser.add_argument("out_dir", metavar="output/dir/path", action="store", help="Output directory path", type=str)

args = parser.parse_args()

beds = args.beds_dir
outdir = args.out_dir
segcopy = pd.DataFrame(columns=["CHR", "START", "END"])
for filename in os.listdir(beds):
    if filename.endswith(".bed"):
         # print(os.path.join(directory, filename))
        bed = os.path.join(beds, filename)
        cellid = os.path.splitext(filename)[0]
        
        df = pd.read_csv(bed, sep="\t", header=None)

        if len(segcopy) == 0: #first cell
            segcopy["CHR"] = df[0]
            segcopy["START"] = df[1]
            segcopy["END"] = df[2]
        
        segcopy[cellid] = df[3]

segcopy.to_csv(os.path.join(outdir, "SegCopy"), sep="\t", index=False)
