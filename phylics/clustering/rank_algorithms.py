#!/usr/bin/envs python

import sys, os
import argparse
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns

#TODO
#retrieve rankings from clusterings, pca, umap and compute an average over the ranks

def main():
    parser = argparse.ArgumentParser(description="Over all clustering ranking.")

    parser.add_argument("-i", "--input_folder", required="True", help="Input folder path")
    parser.add_argument("-o", "--output_folder", required="True", help="output folder path")

    args = parser.parse_args()

    
    exp = ["Exp100", "Exp200", "Exp400"]
    preproc = ["clusterings", "clusterings_pca", "clusterings_umap"]

    for p in preproc:
        merged_ranking = pd.DataFrame()
        for e in exp:
            df = pd.read_csv(os.path.join(args.input_folder, e, p, "clusterings_ranking.csv"), sep="\t", index_col=0)
            merged_ranking[e] = df['rank'].astype(float)
        print(p)
        print(merged_ranking.mean(axis=1).sort_values(ascending=False))
        print("\n")
    
                
            
if __name__ == "__main__":
    sys.exit(main())