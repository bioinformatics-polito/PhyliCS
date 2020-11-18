#!/usr/bin/envs python
import os, sys
import argparse
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    parser = argparse.ArgumentParser(description="Compute stats of clustering metrics.")

    parser.add_argument("-f", "--file", help="Input file.", required=True)
    parser.add_argument("-o", "--out_path", help="Output path.", required=True)

    args =  parser.parse_args()

    pd.set_option('use_inf_as_na', True)
    df = pd.read_csv(args.file, sep="\t", index_col=0)
    

    #df.replace([np.inf, -np.inf], np.nan)
    

    #Compute statistics (count, mean, std, min, max, 25th perc, 50th perc and 75th perc.) for all clusterings
    stats_df = df.describe()
    stats_df.to_csv(os.path.join(args.out_path, "dunn_stats.csv"))

    #Plot mean scores + std deviations
    plt.suptitle("Mean Dunn Index and standard deviation\n",fontsize=16)
    sns.barplot(data=df)
    plt.gcf().set_size_inches(25, 12)
    plt.xlabel('clusterings')
    plt.ylabel('mean dunn index')
    plt.savefig(os.path.join(args.out_path, "dunn_means.png"))
    plt.clf()
    
    #Kruscal-wallis test
    h, p = stats.kruskal(*df.transpose().values, nan_policy='omit')

    #violin plot + bar plot of missing values
    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(25, 18),  gridspec_kw={'height_ratios': [2, 1]})
    ax1.set_title("Dunn Index \nfor different clustering algorithms\nKruskal-Wallis H-test result: s={}, p={}".format(h, p))
    sns.violinplot(data=df, ax=ax1)
    ax1.set_xlabel('clusterings')
    ax2.set_ylabel('dunn index')
    #plt.gcf().set_size_inches(25, 12)
    ax2.set_title("Bad clustering occurrencies")
    max_y = df.isna().sum().max() + 1
    ax2.set_ylim(0, max_y)
    splot = sns.barplot(data=pd.DataFrame(df.isna().sum()).transpose(), ax=ax2)
    for p in splot.patches:
        splot.annotate(format(p.get_height(), '.2f'), (p.get_x() + p.get_width() / 2., p.get_height()), 
        ha = 'center', va = 'center', xytext = (0, 10), textcoords = 'offset points')
    ax2.set_xlabel('clusterings')
    ax2.set_ylabel('Occurrencies')
    fig.savefig(os.path.join(args.out_path, "dunn_violin.png"))
    fig.clf()

if __name__ == "__main__":
    sys.exit(main())
