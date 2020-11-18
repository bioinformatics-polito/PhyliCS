#!/usr/bin/envs python
from math import nan
import sys, os
import argparse
import pandas as pd
import itertools
import itertools

from sklearn import metrics

def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text  
def main():
    parser = argparse.ArgumentParser(description="External clustering validation")

    parser.add_argument("-t", "--true_labels", required=True, help="True labels file",type=str)
    parser.add_argument("-p", "--pred_labels", required=True, help="True labels file", type=str)
    parser.add_argument("-o", "--out_prefix", required=True, help="Output file prefix", type=str)

    args = parser.parse_args()

    labels_true = pd.read_csv(args.true_labels, sep="\t", index_col=0)
    labels_pred = pd.read_csv(args.pred_labels, sep="\t", index_col=0)

    names = []
    for n in labels_pred.index:
        names.append(int(remove_prefix(n, "cell")))
    
    labels_pred.index = names

    labels_true = labels_true[labels_true["ClusterNumber"] != -1]
    labels_pred = labels_pred[labels_pred["cluster"] != -1]

    labels = pd.DataFrame({'labels_true':labels_true["ClusterNumber"]}, index = labels_true.index)
    labels["labels_pred"] = labels_pred["cluster"]
    labels = labels.dropna() #to drop rows of outliers, and check precision only of assigned points
    
    # Precision score

    """
        The precision is intuitively the ability of the classifier not to label 
        as positive a sample that is negative.
        The best value is 1 and the worst value is 0.

        weighted': Calculate metrics for each label, and find their average weighted by 
        support (the number of true instances for each label). 

    """

    
    try:
        prec = metrics.precision_score(labels["labels_true"].values, labels["labels_pred"].values, average='weighted')
    except:
        prec = nan
   
    # Recall score
    """
        The recall is intuitively the ability of the classifier to find all the positive samples.

        The best value is 1 and the worst value is 0.
    """
    try:
        rec = metrics.recall_score(labels["labels_true"].values, labels["labels_pred"].values, average='weighted')
    except:
        rec = nan
   

    #Adjusted Rand index

    """
    Random (uniform) label assignments have a ARI score close to 0.0

    No assumption is made on the cluster structure.

    Bounded range [-1, 1]: negative values are bad (independent labelings), 
    similar clusterings have a positive ARI, 1.0 is the perfect match score.

    """
    try:
        ari = metrics.adjusted_rand_score(labels["labels_true"].values, labels["labels_pred"].values)
    except:
        ari = nan
    # Mutual Information based scores
    """
    Perfect labeling is scored 1.0
    Bad (e.g. independent labelings) have non-positive scores
    """

    try:
        ami = metrics.adjusted_mutual_info_score(labels["labels_true"].values, labels["labels_pred"].values) 
    except:
        ami=nan
    # Fowlkes-Mallows scores
    """
        The Fowlkes-Mallows score FMI is defined as the geometric 
        mean of the pairwise precision and recall.

        The Fowlkes-Mallows score FMI is defined as the geometric mean 
        of the pairwise precision and recall.
    """
    try:
        fmi = metrics.fowlkes_mallows_score(labels["labels_true"].values, labels["labels_pred"].values)
    except:
        fmi = nan
    #V-measure
    """
        Bounded scores: 0.0 is as bad as it can be, 1.0 is a perfect score
    """
    try:
        vm = metrics.v_measure_score(labels["labels_true"].values, labels["labels_pred"].values)
    except:
        vm = nan

    """
    print("Precision\t{}".format(prec))
    print("Recall\t{}".format(rec))
    print("ARI\t{}".format(ari))
    print("AMI\t{}".format(ami))
    print("FMI\t{}".format(fmi))
    print("VM\t{}".format(vm))
    """
    with open(args.out_prefix+'_ext_validation.csv', 'w+') as f:
        f.write("Precision\t{}\n".format(prec))
        f.write("Recall\t{}\n".format(rec))
        f.write("ARI\t{}\n".format(ari))
        f.write("AMI\t{}\n".format(ami))
        f.write("FMI\t{}\n".format(fmi))
        f.write("VM\t{}\n".format(vm))
    
    
if __name__ == "__main__":
    sys.exit(main())