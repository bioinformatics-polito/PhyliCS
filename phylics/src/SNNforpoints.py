from numpy import *
from pylab import *
from Scientific import *
import math
import re
import string
import time
import  msvcrt
from math import sqrt
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics.pairwise import euclidean_distances

parser = argparse.ArgumentParser(description="Graph-based clustering")


parser.add_argument("input", metavar='SegCopy', action='store',
        help='cnvs file.', type=str)
parser.add_argument("outdir", metavar='outdir', action='store',
        help='Output directory path.', type=str)

parser.add_argument("n_pcs", metavar='n_pcs', action='store',
        help='Number of principal components.', type=int)



args=parser.parse_args()

input_file = args.input
outdir = args.outdir
npcs = args.n_pcs

df = pd.read_csv(input_file, sep="\t")

cnvs = df.drop(['CHR', 'START', 'END'], axis=1).transpose()
boundaries = df[['CHR', 'START', 'END']].copy()
chr_limits = boundaries.index[boundaries['END'].isin(boundaries.groupby('CHR', sort=False)['END'].max().values)].tolist()
chr_boundaries = np.append(0, chr_limits)
chr_list = boundaries['CHR'].unique().tolist()
chrN_list = []

for x in chr_list:
    x = x[3:] #remove 'chr' for readability
    chrN_list.append(x)

#compute the position where chromosome labels will be placed on the plots
start = 0
pos_list = []
for end in chr_limits:
    pos_list.append((start+end)/2)
    start = end+1

#preliminary dimensionality reduction
pca = PCA(n_components=npcs)
pca_results = pca.fit_transform(cnvs)

#tsne
N=len(cnvs)
optPerp=round(sqrt(N))
tsne = TSNE(n_components=2, perplexity=optPerp, n_iter=10000)
tsne_results = tsne.fit_transform(cnvs.values)

plt.scatter(*tsne_results.T, s=90, linewidth=0, c=cluster_colors, alpha=0.5)
plt.gcf().suptitle('dbscan clusters', fontsize=16, fontweight='bold' )
plt.gcf().set_size_inches(12,12)
plt.savefig(outdir+'/tsne_projection.png')
plt.clf()
         
"""
knn=zeros((n,n))
for i in range(n):
    for j in range(i,n):
        knn[i,j]=knn[j,i] = sqrt(dot(inputdata[i], inputdata[i]) - 2 * dot(inputdata[i], inputdata[j]) + dot(inputdata[j], inputdata[j]))
        #knn[i,j]=knn[j,i]=(inputdata[i,0]-inputdata[j,0])**2+(inputdata[i,1]-inputdata[j,1])**2
    #knn[i,i]=1000
"""
knn = euclidean_distances(cnvs.values, cnvs.values)

k=3

for i in range(n):
    kvalue=[knn[i,0]]
    krank=[0]
    for j in range(1,n):
        z=0
        while z<len(kvalue):
            if kvalue[z]>=knn[i,j]:
                kvalue[z:z]=[knn[i,j]]
                krank[z:z]=[j]
                break
            z=z+1
        if z==len(kvalue):
            kvalue[z:z]=[knn[i,j]]
            krank[z:z]=[j]
    #print krank
    for j in range(k):
        knn[i,krank[j]]=-1
    for j in range(n):
        if knn[i,j]>=0:
            knn[i,j]=0
        else:
            knn[i,j]=1
    #print "krank:",krank
    #print "kvalue:",kvalue
    kvalue=[]
    krank=[]

#print knn

snn=zeros((n,n))
for i in range(n-1):
    for j in range(i+1,n):
        num=0
        for z in range(n):
            if knn[i,z]==1 and knn[j,z]==1:
                num=num+1
        snn[i,j]=snn[j,i]=num

print('snn:\n',snn)

totalnum=zeros(n)
for i in range(n):
    for j in range(n):
        totalnum[i]+=snn[i,j]
        
print('totalnum:',totalnum)

noise=0
repre=3
strolink=2
clusternum=0
cluster=zeros(n)
merge=[]
for i in range(n):
    if totalnum[i]<=noise:
        for j in range(n):
            if snn[i,j]>0:
                snn[i,j]=0
                snn[j,i]=0
for i in range(n):
    if totalnum[i]>=repre and cluster[i]==0:
        print i
        clusternum+=1
        cluster[i]=clusternum
        for j in range(n):
            if snn[i,j]>=strolink:
                if cluster[j]<>0:
                    merge.append(cluster[j])
                    merge.append(clusternum)
                else:
                    cluster[j]=clusternum

print("cluster:",cluster)
print("merge:",merge)


while len(merge)>0:
    i=merge[0]
    j=merge[1]
    if i>j:
        i,j=j,i
    for z in range(n):
        if cluster[z]==j:
            cluster[z]=i
    merge.remove(i)
    merge.remove(j)
    for each in range(len(merge)):
        if merge[each]==j:
            merge[each]=i
    #print merge
print(cluster)
print(merge)
"""
figure(2)
for i in range(n):
    if cluster[i]==1:
        plot(inputdata[i,0],inputdata[i,1],'ob')
    else:
        plot(inputdata[i,0],inputdata[i,1],'og')
xlim(0,12)
ylim(0,10)        

show()
"""        
        
        
        
        
        
        



