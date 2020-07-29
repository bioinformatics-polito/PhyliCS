import argparse
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="pca")


parser.add_argument("input", metavar='SegCopy', action='store',
        help='cnvs file.', type=str)
parser.add_argument("outdir", metavar='outdir', action='store',
        help='Output directory path.', type=str)

args=parser.parse_args()

input_file = args.input
outdir = args.outdir

df = pd.read_csv(input_file, sep="\t", usecols = lambda column : column not in ['CHR', 'START', 'END']).transpose()

PCs = 50
pca = PCA(n_components=PCs)
pca_result = pca.fit_transform(df.values)

expl_var_ratio = pca.explained_variance_ratio_

expl_var_ratio_perm = pd.DataFrame(columns=df.columns)

n_perm = 50
for k in range(n_perm):
    df_perm = df.apply(np.random.permutation, 0)
    pca_perm = pca.fit_transform(df_perm.values)
    expl_var_ratio_df = pd.DataFrame(data=pca.explained_variance_ratio_.reshape(-1, len(pca.explained_variance_ratio_)))
    expl_var_ratio_perm = expl_var_ratio_perm.append(expl_var_ratio_df)

mean_expl_var_ratio_perm = expl_var_ratio_df.mean().values


plt.figure()
plt.plot(expl_var_ratio, color='green', label='Explained by PCs')
plt.plot(mean_expl_var_ratio_perm, color='red', label='Explained by chance')

plt.gca().set_xlabel('number of components')
plt.gca().set_ylabel('cumulative explained variance')

plt.xticks(np.arange(0, PCs, 1))

plt.gca().legend()
plt.gca().grid(True)


plt.suptitle("Explained variance ratio")

plt.gcf().set_size_inches(21,21)
plt.savefig(outdir+"/expl_var_ratio.png")
plt.gcf().clf()

#cumulative variance
plt.figure()

plt.plot(np.cumsum(expl_var_ratio))
plt.gcf().set_size_inches(21,21)
plt.xticks(np.arange(0, PCs, 1))
plt.gca().grid(True)
plt.xlabel('Number of Components')
plt.ylabel('Variance (%)') #for each component
plt.title('Cumulative explained Variance rario')
plt.savefig(outdir+"/cum_expl_var_ratio2.png")
plt.gcf().clf()


