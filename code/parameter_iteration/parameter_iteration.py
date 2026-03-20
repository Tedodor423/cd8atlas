import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import scipy.sparse as sp
import scrublet as scr
import scipy.stats as stats
from pathlib import Path


# SETTINGS ####
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.figsize"] = [5, 4]
plt.rcParams["savefig.bbox"] = "tight"

OUTPUT_NAME = "parameter_iteration_test"

root_dir = Path("/ceph/project/sharmalab/dnimrich/cd8atlas")
input_h5ad = root_dir / "data" / "huARdb_v2_GEX.CD8.all_genes.h5ad"
params_csv = root_dir / "code" / "parameter_iteration" / "parameters.csv"

figures_dir = root_dir / "output" / OUTPUT_NAME / "figures"
data_dir = root_dir / "output" / OUTPUT_NAME / "data"

figures_dir.mkdir(parents=True, exist_ok=True)
data_dir.mkdir(parents=True, exist_ok=True)

random_seed = 423
min_cells_per_gene = 20
##############


# 0. Load data
adata = sc.read_h5ad(str(input_h5ad))
print(f"Loaded: {input_h5ad} with shape {adata.shape}")

# subsampling
max_cells = 100_000
if max_cells and adata.n_obs > max_cells:
    idx = np.random.default_rng(random_seed).choice(adata.n_obs, size=max_cells, replace=False)
    adata = adata[idx, :].copy()
    print(f"After subsampling: {adata.shape}")

# preserve original counts
adata.layers["counts"] = adata.X.copy() 


# 1. QC
## 1.1 ambient RNA QC
### check if necessary

low_count_treshold = 1_000
frac_treshold = 0.05

low_count_cells = np.asarray(adata.X.sum(axis=1)) < low_count_treshold
if low_count_cells.mean() > frac_treshold:
    print("potential RNA contamination, clean up the data")
    exit()


## 1.2 count QC
### QC metrics

#### mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
#### ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
#### hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt","ribo", "hb"],  # what percentage of counts is mito, ribosomal, hemoglobin
    percent_top=[20, 50, 100],    # what percenteage of counts is top 20, 50, 100 genes
    inplace=True
)

### Get outliers with MAD

def outlier_MAD(adata, metric: str, nmads: int) -> pd.Series:
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * stats.median_abs_deviation(M)) | (
        np.median(M) + nmads * stats.median_abs_deviation(M) < M
    )
    return outlier

adata.obs["outlier"] = (
    outlier_MAD(adata, "log1p_total_counts", 5)
    | outlier_MAD(adata, "log1p_n_genes_by_counts", 5)
    | outlier_MAD(adata, "pct_counts_in_top_20_genes", 5)
)
adata.obs.outlier.value_counts()
adata.obs.outlier.value_counts(normalize=True)

adata.obs["mt_outlier"] = outlier_MAD(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)
adata.obs.mt_outlier.value_counts()
adata.obs.mt_outlier.value_counts(normalize=True)

### plot before QC
sns.displot(adata.obs["total_counts"], bins=100, kde=True)
plt.savefig(figures_dir / "qc_before_total_counts_hist.png")
plt.close()
sns.displot(adata.obs["pct_counts_mt"], bins=100, kde=True)
plt.savefig(figures_dir / "qc_before_pct_counts_mt_hist.png")
plt.close()

sc.pl.violin(adata, "pct_counts_mt", show=False)
plt.savefig(figures_dir / "qc_before_pct_counts_mt_violin.png", bbox_inches="tight")
plt.close()
qc_scatter_genes = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
plt.savefig(figures_dir / "qc_before_scatter_genes_total_mt.png")
plt.close()


### Filter
print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)]
adata = adata.copy()

print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

genes_before_filter = adata.n_vars
sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
print(
    f"Number of genes after min_cells>={min_cells_per_gene} filtering: "
    f"{adata.n_vars} (removed {genes_before_filter - adata.n_vars})"
)

### plot after count QC
sns.displot(adata.obs["total_counts"], bins=100, kde=True).savefig(figures_dir / "qc_after_total_counts_hist.png")
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
plt.savefig(figures_dir / "qc_after_scatter_genes_total_mt.png")
plt.close()



## 1.3 Doublet QC
n_obs_before_doublet_qc = adata.n_obs
### prepare input
counts_csr = adata.layers["counts"].tocsr() if not sp.issparse(adata.layers["counts"]) else adata.layers["counts"].tocsr()
scrublet = scr.Scrublet(counts_matrix=counts_csr, expected_doublet_rate=0.05)

### run scrublet
doublet_scores, predicted_doublets = scrublet.scrub_doublets()

### save scores, calls, histogram
adata.obs["doublet_score"] = doublet_scores
adata.obs["doublet_call"] = pd.Categorical(predicted_doublets)
adata.uns["scrublet_sim_scores"] = scrublet.doublet_scores_sim_

scrublet.plot_histogram()
plt.savefig(figures_dir / "qc_scrublet_histogram.png", dpi=300)
plt.close()

### Filter if necessary
if adata.obs["doublet_call"].nunique() > 1:
    print("There are predicted doublets, applying filter")
    adata = adata[~adata.obs["doublet_call"].astype(bool).to_numpy(), :]
else:
    print("No doublets predicted")

print(f"After doublet filtering: {adata.shape} (kept {adata.n_obs}/{n_obs_before_doublet_qc} cells)")


# 2. Normalisation
# ## retrieve count values - not necessary now
# adata.X = adata.layers["counts"].copy() 

sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
adata.layers["normalised"] = adata.X.copy()

sc.pp.log1p(adata)
adata.X = adata.X.astype(np.float32)  # for better storage
adata.layers["log1p"] = adata.X.copy()

log1p_library_size = np.asarray(adata.layers["log1p"].sum(axis=1)).ravel()
log1p_library_size


sns.histplot(adata.obs["total_counts"], bins=100, kde=False)
plt.savefig(figures_dir / "normalisation_total_counts_hist.png")
plt.close()
sns.histplot(log1p_library_size, bins=100, kde=False)
plt.savefig(figures_dir / "normalisation_log1p_hist.png")
plt.close()

#### Parameter Iteration ####

## load parameters from csv
params_df = pd.read_csv(params_csv)
param_names = ["hvg_n_top_genes", "n_pcs", "n_neighbors", "leiden_resolution", "min_dist"]
param_types = {"hvg_n_top_genes": int, "n_pcs": int, "n_neighbors": int, "leiden_resolution": float, "min_dist": float}
params = {}
for name in param_names:
    params[name] = params_df[name].dropna().astype(param_types[name]).tolist()

for name in param_names:
    print(f"{name}: {params[name]}")
print("total runs:", np.prod([len(params[name]) for name in param_names]))


## iteration
cluster_statistics = []
run_id = 0

for hvg_n_top_genes in params["hvg_n_top_genes"]:
    for n_pcs in params["n_pcs"]:
        for n_neighbors in params["n_neighbors"]:
            for leiden_resolution in params["leiden_resolution"]:
                for min_dist in params["min_dist"]:
                    run_id += 1
                    print(
                        f"Run {run_id}: hvg={hvg_n_top_genes}, n_pcs={n_pcs}, n_neighbors={n_neighbors}, resolution={leiden_resolution}, min_dist={min_dist}"
                    )

                    adata_i = adata.copy()

                    sc.pp.highly_variable_genes(
                        adata_i,
                        layer="log1p",
                        flavor="seurat",
                        n_top_genes=hvg_n_top_genes,
                        inplace=True,
                    )
                    adata_i = adata_i[:, adata_i.var["highly_variable"]].copy()

                    sc.pp.scale(adata_i, max_value=10)
                    sc.tl.pca(adata_i, n_comps=n_pcs, svd_solver="arpack", random_state=random_seed)
                    sc.pp.neighbors(adata_i, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_seed)
                    sc.tl.umap(adata_i, min_dist=min_dist, spread=1.0, random_state=random_seed)

                    leiden_key = f"leiden_res_{leiden_resolution}"
                    sc.tl.leiden(adata_i, resolution=leiden_resolution, key_added=leiden_key)

                    # only one UMAP plot per parameter combination
                    sc.pl.umap(adata_i, color=[leiden_key], show=False)
                    safe_tag = (
                        f"hvg{hvg_n_top_genes}_pcs{n_pcs}_nn{n_neighbors}_"
                        f"res{leiden_resolution:.2f}_md{min_dist:.2f}"
                    ).replace(".", "p")
                    plt.savefig(figures_dir / f"umap_{safe_tag}.png")
                    plt.close()

                    cluster_counts = adata_i.obs[leiden_key].value_counts()

                    mean_expr = np.asarray(adata_i.X.mean(axis=0)).ravel()
                    top_idx = np.argsort(mean_expr)[::-1][:20]
                    top_genes = ";".join(adata_i.var_names[top_idx].tolist())

                    cluster_statistics.append(
                        {
                            "run_id": run_id,
                            "hvg_n_top_genes": hvg_n_top_genes,
                            "n_pcs": n_pcs,
                            "n_neighbors": n_neighbors,
                            "leiden_resolution": leiden_resolution,
                            "min_dist": min_dist,
                            "n_cells": int(adata_i.n_obs),
                            "n_genes_after_hvg": int(adata_i.n_vars),
                            "n_clusters": int(cluster_counts.shape[0]),
                            "smallest_cluster_n": int(cluster_counts.min()),
                            "largest_cluster_fraction": float(cluster_counts.max() / adata_i.n_obs),
                            "top20_highly_expressed_genes": top_genes,
                        }
                    )


# 5. Save statistics
stats_df = pd.DataFrame(cluster_statistics)
stats_df.to_csv(data_dir / "cluster_parameter_stats.csv", index=False)
print(f"Saved: {data_dir / 'cluster_parameter_stats.csv'}")
