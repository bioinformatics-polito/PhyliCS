import phylics
import pandas as pd
breast_a = phylics.Sample.from_file("/scratch/trcanmed/phylics/dataset/10x_breastA/SegCopy", sample_name="S_A")
df = pd.read_csv("/scratch/trcanmed/phylics/dataset/10x_breastA/results.txt", sep="\t", index_col=0).dropna()
ploidy = df['Copy_Number']
breast_a.add_annotation(ploidy, "mean_cn")
breast_a.mad()
breast_a.variable_features()
breast_a = breast_a.filter("mean_cn", "out_eq", (1.5, 3))
breast_a = breast_a.filter("mad", "lt_eq", 0.95, percentile=True)
#s_filter.umap(use_highly_variable=True)
#s_filter.plot_umap(outpath="data/10x_breastA/umap.png")
#s_filter.cluster("hdbscan", embeddings="umap")
#s_filter.plot_clusters("heatmap", "data/10x_breastA/hdbscan_clusters.png")
breast_b = phylics.Sample.from_file("/scratch/trcanmed/phylics/dataset/10x_breastB/SegCopy", sample_name="S_B")
df = pd.read_csv("/scratch/trcanmed/phylics/dataset/10x_breastB/results.txt", sep="\t", index_col=0).dropna()
ploidy = df['Copy_Number']
breast_b.add_annotation(ploidy, "mean_cn")
breast_b.mad()
breast_b.variable_features()
breast_b = breast_b.filter("mean_cn", "out_eq", (1.5, 3))
breast_b = breast_b.filter("mad", "lt_eq", 0.95, percentile=True)
#breast_b.umap(use_highly_variable=True)
#breast_b.cluster("hdbscan", embeddings="umap",  min_cluster_size=15, metric='manhattan')
#breast_b.plot_clusters("scatter", "data/10x_breastB/hdbscan_clusters_scatter.png")

breast_c = phylics.Sample.from_file("/scratch/trcanmed/phylics/dataset/10x_breastC/SegCopy", sample_name="S_C")
df = pd.read_csv("/scratch/trcanmed/phylics/dataset/10x_breastC/results.txt", sep="\t", index_col=0).dropna()
ploidy = df['Copy_Number']
breast_c.add_annotation(ploidy, "mean_cn")
breast_c.mad()
breast_c.variable_features()
breast_c = breast_c.filter("mean_cn", "out_eq", (1.5, 3))
breast_c = breast_c.filter("mad", "lt_eq", 0.95, percentile=True)

breast_d = phylics.Sample.from_file("/scratch/trcanmed/phylics/dataset/10x_breastD/SegCopy", sample_name="S_D")
df = pd.read_csv("/scratch/trcanmed/phylics/dataset/10x_breastD/results.txt", sep="\t", index_col=0).dropna()
ploidy = df['Copy_Number']
breast_d.add_annotation(ploidy, "mean_cn")
breast_d.mad()
breast_d.variable_features()
breast_d = breast_d.filter("mean_cn", "out_eq", (1.5, 3))
breast_d = breast_d.filter("mad", "lt_eq", 0.95, percentile=True)

breast_e = phylics.Sample.from_file("/scratch/trcanmed/phylics/dataset/10x_breastE/SegCopy", sample_name="S_E")
df = pd.read_csv("/scratch/trcanmed/phylics/dataset/10x_breastE/results.txt", sep="\t", index_col=0).dropna()
ploidy = df['Copy_Number']
breast_e.add_annotation(ploidy, "mean_cn")
breast_e.mad()
breast_e.variable_features()
breast_e = breast_e.filter("mean_cn", "out_eq", (1.5, 3))
breast_e = breast_e.filter("mad", "lt_eq", 0.95, percentile=True)

ss = phylics.MultiSample.from_list(breast_b, breast_c, breast_d, breast_e)
ss.SHscore(n_jobs=10)
#ss.plot_scores("/scratch/trcanmed/phylics/dataset/10x_multi_sample/scores.png")
ss.plot_dendrogram(outpath="/scratch/trcanmed/phylics/dataset/10x_multi_sample/agglomerative_heatmap.png")

navin_prim = phylics.Sample.from_file("/scratch/trcanmed/phylics/dataset/navin_primary/SegCopy", sample_name="primary")
df = pd.read_csv("/scratch/trcanmed/phylics/dataset/navin_primary/results.txt", sep="\t", index_col=0).dropna()
ploidy = df['Copy_Number']
navin_prim.add_annotation(ploidy, "mean_cn")
navin_prim.mad()
navin_prim.variable_features()
navin_prim = navin_prim.filter("mean_cn", "gt_eq", 3)
navin_prim = navin_prim.filter("mad", "lt_eq", 0.95, percentile=True)

navin_met = phylics.Sample.from_file("/scratch/trcanmed/phylics/dataset/navin_metastasis/SegCopy", sample_name="metastasis")
df = pd.read_csv("/scratch/trcanmed/phylics/dataset/navin_metastasis/results.txt", sep="\t", index_col=0).dropna()
ploidy = df['Copy_Number']
navin_met.add_annotation(ploidy, "mean_cn")
navin_met.mad()
navin_met.variable_features()
navin_met = navin_met.filter("mean_cn", "gt_eq", 3)
navin_met = navin_met.filter("mad", "lt_eq", 0.95, percentile=True)
#navin_met.umap(use_highly_variable=True)
#navin_met.cluster("hdbscan", embeddings="umap",  min_cluster_size=15, min_samples=1, cluster_selection_epsilon=0.5)
#navin_met.plot_clusters("scatter", "data/navin_metastasis/hdbscan_clusters_scatter.png")

ss = phylics.MultiSample.from_list(navin_prim, navin_met)
ss.SHscore(n_jobs=10)
#ss.plot_scores("data/navin_met_prim/scores.png")
ss.plot_dendrogram(outpath="/scratch/trcanmed/phylics/dataset/navin_met_prim/agglomerative_heatmap.png")

