#!/usr/bin/envs python
import os, sys
import argparse
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import minmax_scale



def main():
    parser = argparse.ArgumentParser(description="Over all clustering ranking.")

    parser.add_argument("-i", "--input_folder", required="True", help="Input folder path")
    parser.add_argument("-o", "--output_folder", required="True", help="output folder path")

    args = parser.parse_args()

    
    exp = ["Exp100", "Exp200", "Exp400"]
    preproc = ["clusterings", "clusterings_pca", "clusterings_umap"]

    avg = pd.DataFrame()
    for p in preproc:
        avg_score = pd.DataFrame()
        for e in exp:
            df = pd.read_csv(os.path.join(args.input_folder, e, p, "clusterings_ranking.csv"), sep="\t", index_col=0)
            df = df.drop('rank', axis=1)
            df['ari'] = minmax_scale(df['ari'], feature_range=(0, 1))
            avg_score[e] = df.mean(axis=1)
        avg[p] = avg_score.mean(axis=1)
    avg['rank'] = avg[preproc].rank().mean(axis=1)
    avg = avg.sort_values(by='rank', ascending=False)
    avg.to_csv(os.path.join(args.output_folder, "cross_merge.csv"), sep="\t")
if __name__ == "__main__":
    sys.exit(main())