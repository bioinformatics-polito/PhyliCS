#!/usr/bin/env python

import glob
import sys
import os
import argparse
import pandas as pd 
from math import nan 
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy.stats import mannwhitneyu, shapiro, norm


#https://towardsdatascience.com/prepare-dinner-save-the-day-by-calculating-confidence-interval-of-non-parametric-statistical-29d031d079d0

def non_param_unpaired_CI(sample1, sample2, conf):
  n1 = len(sample1)  
  n2 = len(sample2)  
  alpha = 1-conf      
  N = norm.ppf(1 - alpha/2) 

  # The confidence interval for the difference between the two population
  # medians is derived through the n x m differences.
  diffs = sorted([i-j for i in sample1 for j in sample2])
  
  # the Kth smallest to the Kth largest of the n x m differences then determine 
  # the confidence interval, where K is:
  k = np.math.ceil(n1*n2/2 - (N * (n1*n2*(n1+n2+1)/12)**0.5))
  
  CI = (round(diffs[k-1],3), round(diffs[len(diffs)-k],3))
  return CI

def collect_het_scores(root, exps, samples):
    het_scores = pd.DataFrame(columns=["het", "hom"]) #het and hom are inverted with respect to the filesystem
    for exp in exps:
        prefix = os.path.join(root, exp, "data")
        for sample in samples:
            het_score = nan 
            hom_score = nan
            subsample_dir = os.path.join(prefix, sample, "subsamples")
            f =  glob.glob("{dir}/hom/*_postCNV/heterogeneity_scores.csv".format(dir=subsample_dir))
            if len(f) != 0:
                het_score = pd.read_csv(f[0], sep = "\t").iloc[0]["het_score"].astype(float)
            f = glob.glob("{dir}/het/*_postCNV/heterogeneity_scores.csv".format(dir=subsample_dir))
            if len(f) != 0:
                hom_score = pd.read_csv(f[0], sep = "\t").iloc[0]["het_score"].astype(float)
            het_scores = het_scores.append({"het":het_score, "hom":hom_score}, ignore_index=True)
    return het_scores 

def compute_stats_and_plot(het_scores, outpath):
    stat, p1 =  mannwhitneyu(het_scores['het'].dropna(), het_scores['hom'].dropna())
    _, p2 = shapiro( het_scores['hom'].dropna())
    CI = non_param_unpaired_CI(het_scores["hom"].dropna(), het_scores["het"].dropna(), 0.05)
    CI_hom = []
    CI_hom.append(het_scores["hom"].dropna().quantile(0.5))
    CI_hom.append(het_scores["hom"].dropna().quantile(0.95))
    print("Hom CI=({l},{r})".format(l=round(CI_hom[0],4), r=round(CI_hom[1],4)))
    CI_het = []
    CI_het.append(het_scores["het"].dropna().quantile(0.5)) 
    CI_het.append(het_scores["het"].dropna().quantile(0.95))
    print("Het CI=({l},{r})".format(l=round(CI_het[0], 4), r=round(CI_het[1], 4)))
    print("Shapiro test for 'hom' normality p: {p2}".format(p2=p2))
    print("CI={CI}".format(CI=CI))
    sns.violinplot(data=het_scores)
    textstr = '\n'.join((
        r'Mann-Whitney-U test pval=%.3E' % (p1, ),
        r'Shapiro test for hom normality pval=%.5f' % (p2, ),
        r'CI=%s' % (CI, )))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
    plt.savefig(os.path.join(outpath, "aggregate_het_scores_violin.png"))

def main():
    parser = argparse.ArgumentParser(description="Aggregates heterogeneity score results.")
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory")

    args = parser.parse_args()

    EXPS = ["Exp1", "Exp2", "Exp3"]

    SAMPLES=["Sample1", "Sample2", "Sample3", "Sample4", "Sample5",
        "Sample6", "Sample7", "Sample8", "Sample9", "Sample10",
        "Sample11", "Sample12", "Sample13", "Sample14", "Sample15",
        "Sample16", "Sample17", "Sample18", "Sample19", "Sample20",
        "Sample21", "Sample22", "Sample23", "Sample24", "Sample25",
        "Sample26", "Sample27", "Sample28", "Sample29", "Sample30",
        "Sample31", "Sample32", "Sample33", "Sample34", "Sample35",
        "Sample36", "Sample37", "Sample38", "Sample39", "Sample40",
        "Sample41", "Sample42", "Sample43", "Sample44", "Sample45",
        "Sample46", "Sample47", "Sample48", "Sample49", "Sample50",
        ] 
    het_scores = collect_het_scores(args.input_dir, EXPS, SAMPLES)
    compute_stats_and_plot(het_scores, args.output_dir)
if __name__ == "__main__":
    sys.exit(main())