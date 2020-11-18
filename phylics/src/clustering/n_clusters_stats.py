#!/usr/bin/envs python
import os, sys
import argparse
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    parser = argparse.ArgumentParser(description="Compute stats of clustering metrics.")

    parser.add_argument("-f", "--file", help="Input file.", required=True)
    parser.add_argument("-o", "--out_path", help="Output path.", required=True)

    args =  parser.parse_args()

    df = pd.read_csv(args.file, sep="\t", index_col=0)
    #Compute statistics (count, mean, std, min, max, 25th perc, 50th perc and 75th perc.) for all clusterings
    stats_df = df.describe()
    stats_df.to_csv(os.path.join(args.out_path, "comp_time_stats.csv"))

    #Plot mean scores + std deviations
    plt.suptitle("Mean computation time and standard deviation\n",fontsize=16)
    sns.barplot(data=df)
    plt.gcf().set_size_inches(25, 12)
    plt.xlabel('clusterings', fontsize=14)
    plt.ylabel('mean clusters number', fontsize=14)
    plt.savefig(os.path.join(args.out_path, "clusters_n_means.png"))
    plt.clf()


if __name__ == "__main__":
    sys.exit(main())
