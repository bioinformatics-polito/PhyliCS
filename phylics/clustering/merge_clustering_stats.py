#!/usr/bin/envs python
import os, sys, glob
import argparse
import pandas as pd 
from math import nan
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

def compute_stats_and_save(df, index_name, out_path):
    #Compute statistics (count, mean, std, min, max, 25th perc, 50th perc and 75th perc.) for all clusterings
    stats_df = df.describe()
    stats_df.to_csv(os.path.join(out_path, index_name + "_stats.csv"))

    #Plot mean scores + std deviations
    plt.suptitle("Mean " + index_name + " and standard deviation\n",fontsize=16)
    sns.barplot(data=df)
    plt.gcf().set_size_inches(25, 12)
    plt.xlabel('clusterings')
    plt.ylabel('mean ' + index_name)
    plt.savefig(os.path.join(out_path, index_name + "_means.png"))
    plt.clf()

    #Kruscal-wallis test
    h, p = stats.kruskal(*df.transpose().values, nan_policy='omit')

    with open(os.path.join(out_path, index_name + "_kruskal.csv"), 'w+') as f:
        f.write("kruskal_wallis_p\t" + str(p) + "\n")

    #violin plot 
    plt.suptitle(index_name + " distribution for different clustering algorithms\nKruskal-Wallis H-test result: s={}, p={}".format(h, p), fontsize=16)
    sns.violinplot(data=df)
    plt.xlabel('clusterings')
    plt.ylabel(index_name)
    plt.gcf().set_size_inches(25, 12)
    plt.savefig(os.path.join(out_path, index_name + "_violin.png"))
    plt.clf()

def compute_and_plot_heatmap(stats_dict, clusterings, out_path):
    stats_df = pd.DataFrame(index=stats_dict.keys(), columns = clusterings)
    for index in stats_dict.keys():
        means = stats_dict[index].mean()
        stats_df.loc[index] = means

    stats_df = stats_df.astype(float).round(3)
    stats_df = stats_df.transpose()

    # Compute a ranking  of clustering algorithms along each ext index
    # Then sort them, in descending order, according to the average rank
    # computed on all indices
    stats_df['rank'] = stats_df[['ari', 'ami', 'fmi', 'vm']].rank().mean(axis=1)
    stats_df = stats_df.sort_values(by='rank', ascending=False)
    stats_df.to_csv(os.path.join(out_path, "clusterings_ranking.csv"), sep='\t')
    stats_df = stats_df.drop('rank', axis=1)
    sns.heatmap(stats_df, cmap="coolwarm", vmin=0, vmax=1, center=0.5, fmt='.3g',
        linewidth=0,annot=True, square=True)
    #cg.ax_row_dendrogram.set_visible(False)
    #cg.ax_col_dendrogram.set_visible(False)
    plt.savefig(os.path.join(out_path,  "validation_indices_heatmap.png"),  bbox_inches = 'tight',
    pad_inches = 0)
    plt.clf()

def main():
    parser = argparse.ArgumentParser(description="Merge clustering statistics.")
    parser.add_argument("-d", "--data_folder", help="Data folder path.", required=True)
    parser.add_argument("-s","--simulations", help="Simulation number.", required=True, type=int)
    parser.add_argument("-c", "--clusterings", nargs="+", help="Clustering algorithm list.", required=True)
    parser.add_argument("-o", "--out_folder", help="Output folder.", required=True)

    args = parser.parse_args()

    data_folder = args.data_folder
    simulations = args.simulations
    clusterings = args.clusterings
    out_folder = args.out_folder


    samples = []
    for i in range(1, simulations + 1):
        samples.append("Sample" + str(i))
    ari = pd.DataFrame(columns=clusterings, index=samples)
    ami = pd.DataFrame(columns=clusterings, index=samples)
    fmi = pd.DataFrame(columns=clusterings, index=samples)
    vm = pd.DataFrame(columns=clusterings, index=samples)

    apn = pd.DataFrame(columns=clusterings, index=samples)
    time = pd.DataFrame(columns=clusterings, index=samples)

    for c in clusterings:
        ari_ = pd.Series(index=samples)
        ami_ = pd.Series(index=samples)
        fmi_ = pd.Series(index=samples)
        vm_ = pd.Series(index=samples)

        apn_ = pd.Series(index=samples)
        time_ = pd.Series(index=samples)

        prefix = os.path.join(data_folder, c)
        for s in samples:
            filepath = os.path.join(prefix, s + "_" + c + "_ext_validation.csv")
            stats_ext = pd.read_csv(filepath, index_col=0, sep="\t", header=None).transpose()
            stats_ext = stats_ext.loc[:,~stats_ext.columns.duplicated()]
            ari_.loc[s] = stats_ext["ARI"].values[0]
            ami_.loc[s] = stats_ext["AMI"].values[0]
            fmi_.loc[s] = stats_ext["FMI"].values[0]
            vm_.loc[s] = stats_ext["VM"].values[0]

            filepath = os.path.join(prefix, s + "_" + c + "_performance.csv")
            stats = pd.read_csv(filepath, index_col=0, sep="\t", header=None).transpose()
            stats = stats.loc[:,~stats_ext.columns.duplicated()]
            apn_.loc[s] = stats["APN"].values[0]
            time_.loc[s] = stats["computation_time"].values[0]
        ari[c] = ari_
        ami[c] = ami_
        fmi[c] = fmi_
        vm[c] = vm_

        apn[c] = apn_
        time[c] = time_

    ari.to_csv(os.path.join(out_folder, "ARI.csv"), sep="\t")
    ami.to_csv(os.path.join(out_folder, "AMI.csv"), sep="\t")
    fmi.to_csv(os.path.join(out_folder, "FMI.csv"), sep="\t")
    vm.to_csv(os.path.join(out_folder, "VM.csv"), sep="\t")

    apn.to_csv(os.path.join(out_folder, "apn_scores.csv"), sep="\t")
    time.to_csv(os.path.join(out_folder, "computation_times.csv"), sep="\t")

    stats_dict = {"ari":ari, "ami":ami, "fmi":fmi, "vm":vm}
    sns.set(font_scale=1.8)
    # Compute stats over all indices
    
    compute_stats_and_save(ari, "ARI", out_folder)
    compute_stats_and_save(ami, "AMI", out_folder)
    compute_stats_and_save(fmi, "FMI", out_folder)
    compute_stats_and_save(vm, "VM", out_folder)
    #compute_stats_and_save(silhouette_scores, "silhouette", out_folder)
    
    # Plot heatmap
    compute_and_plot_heatmap(stats_dict, clusterings, out_folder)
    
    
if __name__ == "__main__":
    sys.exit(main())