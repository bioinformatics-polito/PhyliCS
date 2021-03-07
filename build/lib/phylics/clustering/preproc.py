#|/usr/bin/envs python

import argparse
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.feature_selection import VarianceThreshold

parser = argparse.ArgumentParser(description="pca")


parser.add_argument("input", metavar='SegCopy', action='store',
        help='cnvs file.', type=str)
parser.add_argument("outdir", metavar='outdir', action='store',
        help='Output directory path.', type=str)

args=parser.parse_args()

input_file = args.input
outdir = args.outdir

df = pd.read_csv(input_file, sep="\t", usecols = lambda column : column not in ['CHR', 'START', 'END']).transpose()

# Low variance feature filtering 
selector = VarianceThreshold() #filters out features with 0 variance from cell to cell (not informative for clustering)
df_filtered = pd.DataFrame(data=selector.fit_transform(df), index=df.index)


# PCA (keep only significant features)
#n_samples = len(df)
#n_features = len(df.columns)
PCs = 50
pca = PCA(n_components=PCs)
pca_result = pd.DataFrame(data=pca.fit_transform(df_filtered.values), index=df.index)

expl_var_ratio = pca.explained_variance_ratio_
expl_var_ratio_df = pd.DataFrame(data=pca.explained_variance_ratio_.reshape(-1, len(pca.explained_variance_ratio_)))

expl_var_ratio_perm = pd.DataFrame() #columns = components, rows = permutations
n_perm = 50
for k in range(n_perm):
    df_perm = df_filtered.apply(np.random.permutation, 0)
    pca_perm = pca.fit_transform(df_perm.values)
    expl_var_ratio_single_df = pd.DataFrame(data=pca.explained_variance_ratio_.reshape(-1, len(pca.explained_variance_ratio_)))
    expl_var_ratio_perm = expl_var_ratio_perm.append(expl_var_ratio_single_df)
mean_expl_var_ratio_perm = expl_var_ratio_perm.mean().values


plt.figure()
plt.plot(expl_var_ratio, color='green', label='Explained by PCs')
plt.plot(mean_expl_var_ratio_perm, color='red', label='Explained by chance')

plt.gca().set_xlabel('number of components',  fontsize=18)
plt.gca().set_ylabel('cumulative explained variance',  fontsize=18)

plt.xticks(np.arange(0, PCs, 1))

plt.gca().legend(fontsize=18)
plt.gca().grid(True)


plt.suptitle("Explained variance ratio",  fontsize=22)
plt.gca().tick_params(axis="x", labelsize=12)
plt.gca().tick_params(axis="y", labelsize=12)
plt.gcf().set_size_inches(15, 12)
plt.savefig(outdir+"/expl_var_ratio.png")
plt.gcf().clf()


# compute a pvalue  of how the observed variance is different from the permuted variance for each PC

for c in expl_var_ratio_df.columns:
    expl_var_ratio_perm.loc[expl_var_ratio_perm[c] > expl_var_ratio_df.loc[0, c], c] = 1
    expl_var_ratio_perm.loc[expl_var_ratio_perm[c] <= expl_var_ratio_df.loc[0, c], c] = 0

pvals = expl_var_ratio_perm.mean()

#optPCs = pvals[pvals.gt(0.05)].index[0]

# find the interection point between the two curves
pcs = np.argwhere(np.diff(np.sign(expl_var_ratio - mean_expl_var_ratio_perm))).flatten()

if len(pcs) == 1:
    optPCs = pcs[0]
else:
    optPCs = optPCs[2] #when the noisy curves fluctuates 

plt.figure()
plt.plot(pvals, color='red')

plt.axvline(x=optPCs, color='blue', linestyle='--', label="optPCs")

plt.gca().set_xlabel('principal components', fontsize=18)
plt.gca().set_ylabel('pvalue', fontsize=18)

plt.xticks(np.arange(0, PCs, 1))
plt.gca().grid(True)
plt.legend(fontsize=18)
plt.suptitle("Significance of principal components", fontsize=22)
plt.gca().tick_params(axis="x", labelsize=12)
plt.gca().tick_params(axis="y", labelsize=12)
plt.gcf().set_size_inches(30, 12)
plt.savefig(outdir+"/pcs_pvals.png")
plt.gcf().clf()

#cumulative variance
plt.figure()

plt.plot(np.cumsum(expl_var_ratio))
plt.gcf().set_size_inches(30, 12)
plt.xticks(np.arange(0, PCs, 1))
plt.gca().grid(True)
plt.gca().tick_params(axis="x", labelsize=12)
plt.gca().tick_params(axis="y", labelsize=12)
plt.xlabel('Number of Components', fontsize=18)
plt.ylabel('Variance (%)', fontsize=18) #for each component
plt.title('Cumulative explained Variance ratio', fontsize=22)
plt.savefig(outdir+"/cum_expl_var_ratio.png")
plt.gcf().clf()

pca_result[pca_result.columns[0:optPCs]].to_csv(outdir+"/principal_components.csv", sep="\t")

