import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scrublet as scr
# import celltypist
# from celltypist import models
import scipy.stats as stats
from pathlib import Path
import bioreport


# region SETTINGS ####
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.figsize"] = [5, 4]
plt.rcParams["savefig.bbox"] = "tight"

OUTPUT_NAME = "test_pipeline"

root_dir = Path("/ceph/project/sharmalab/dnimrich/cd8atlas")

input_h5ad = root_dir / "data/" / "huARdb_v2_GEX.CD8.all_genes.h5ad" #"huARdb_v2_GEX.CD8.hvg4k.h5ad"

figures_dir = root_dir / "output" / OUTPUT_NAME / "figures"
data_dir = root_dir / "output" / OUTPUT_NAME / "data"

figures_dir.mkdir(parents=True, exist_ok=True)
data_dir.mkdir(parents=True, exist_ok=True)

random_seed = 423
min_cells_per_gene = 20
# endregion #############

# region 0. load/subsample data
adata = sc.read_h5ad(str(input_h5ad))
print(f"Loaded: {input_h5ad} with shape {adata.shape}")

## subsampling

subsample_amount = 100_000
if subsample_amount and adata.n_obs > subsample_amount:
    idx = np.random.default_rng(random_seed).choice(adata.n_obs, size=subsample_amount, replace=False)
    adata = adata[idx, :]
    print(f"After subsampling: {adata.shape}")

## preserve original counts
adata.layers["counts"] = adata.X.copy() 

# endregion

# ## save subsampled data
# adata.write_h5ad(data_dir / f"subsampled_{subsample_amount}.h5ad")

# ## retrieve subsampled data
# adata = sc.read_h5ad(data_dir / f"subsampled_100000.h5ad")
# print(f"Loaded data with shape {adata.shape}")

# region 1. QC

# region 1.1 ambient RNA QC
### check if necessary

low_count_treshold = 1_000
frac_treshold = 0.05

low_count_cells = np.asarray(adata.X.sum(axis=1)) < low_count_treshold
if low_count_cells.mean() > frac_treshold:
    print("potential RNA contamination, clean up the data")
    exit()

# endregion

# region 1.2 count QC
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

### Filter based on counts
print(f"Size of matrix before count filtering: {adata.shape}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)]

### Filter out mt, ribo, hb genes
adata = adata[:, ~(adata.var.mt | adata.var.ribo | adata.var.hb)]
print(f"Size of matrix after count filtering: {adata.shape}")

### Filter genes with too few detected cells before HVG selection
genes_before_filter = adata.n_vars
sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
print(
    f"Number of genes after min_cells >= {min_cells_per_gene} filtering: "
    f"{adata.n_vars} (removed {genes_before_filter - adata.n_vars})"
)

### plot after count QC
sns.displot(adata.obs["total_counts"], bins=100, kde=True).savefig(figures_dir / "qc_after_total_counts_hist.png")
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
plt.savefig(figures_dir / "qc_after_scatter_genes_total_mt.png")
plt.close()

# endregion

# region 1.3 Doublet QC
n_obs_before_doublet_qc = adata.n_obs  # for statistics
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
print(f"After doublet step: {adata.shape} (kept {adata.n_obs}/{n_obs_before_doublet_qc} cells)")

# endregion
# endregion

# region 2. Normalisation

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

# endregion

## save progress
adata.write_h5ad(data_dir / "qc+normalised.h5ad")

# ## retrieve progress
# adata = sc.read_h5ad(data_dir / "qc+normalised.h5ad")
# print(f"Loaded data with shape {adata.shape}")


# labeling, CC analysis

# region 3. Feature selection
## 3.1 Cell selection - none here

## 3.2 Gene exclusion

excluded_gene_sets = {
    "stress_genes": [
        "G0S2", "JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "CDKN1A"
    ],
    "cc_genes": [
        "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "TYMS",
        "FEN1", "MKI67", "TOP2A", "TK1", "RRM1", "RRM2", "HELLS", "UNG",
        "GINS2", "CDC6", "CDK1", "CDC20", "CCNB1", "CCNB2", "AURKA", "AURKB",
        "BUB1", "BUB1B", "UBE2C", "CENPF", "CENPE", "HMGB2",
    ],
    "bad_features": ["MALAT1", "NEAT1"],
}

excluded_gene_patterns = {
    "IFN_genes": r"^(IFI|IFIT|IFITM|ISG|MX|OAS|IRF|STAT1$|STAT2$|RSAD2$|BST2$)",
    "ccl_genes": r"^CCL",
    "MHC_genes": r"^(HLA-[A-Z]|B2M$|TAP1$|TAP2$|PSMB8$|PSMB9$)",
    "hist_genes": r"^HIST",
    "comp_genes": r"^(C1Q|C1R$|C1S$|C2$|C3$|CFB$|CFD$|C7$)",
    "ig_genes": r"^(IGH|IGK|IGL|JCHAIN$|MZB1$)",
    "hb_genes": r"^HB(?!P)",
}

gene_names_upper = adata.var_names.str.upper()
excluded_gene_mask = pd.Series(False, index=adata.var_names, dtype=bool)
excluded_gene_summary = {}

for gene_set_name, gene_list in excluded_gene_sets.items():
    current_mask = gene_names_upper.isin(gene_list)
    excluded_gene_mask |= current_mask
    excluded_gene_summary[gene_set_name] = int(current_mask.sum())

for gene_set_name, pattern in excluded_gene_patterns.items():
    current_mask = gene_names_upper.str.contains(pattern, regex=True)
    excluded_gene_mask |= current_mask
    excluded_gene_summary[gene_set_name] = int(current_mask.sum())

print("Genes matched by exclusion panel:")
for gene_set_name, n_genes in excluded_gene_summary.items():
    print(f"  {gene_set_name}: {n_genes}")

print(f"Removing {int(excluded_gene_mask.sum())} genes from exclusion panel before HVG selection")
adata.var["excluded_feature"] = excluded_gene_mask.to_numpy()
adata = adata[:, ~adata.var["excluded_feature"]].copy()
print(f"Shape after exclusion panel filtering: {adata.shape}")



## 3.3. HVG selection
### using seurat as seurat_v3 with original counts raises an environment error
sc.pp.highly_variable_genes(adata, layer='log1p', flavor='seurat', n_top_genes=4000, inplace=True)
print("HVG selected:", int(adata.var["highly_variable"].sum()))

#sc.pl.highly_variable_genes(adata, show=False)
sns.scatterplot(data=adata.var, x="means", y="dispersions", hue="highly_variable", s=5)
plt.savefig(figures_dir / "hvg_before.png")
plt.close()

### Filter
adata = adata[:, adata.var.highly_variable]
 
#sc.pl.highly_variable_genes(adata, show=False)
sns.scatterplot(data=adata.var, x="means", y="dispersions", hue="highly_variable", s=5)
plt.savefig(figures_dir / "hvg_after.png")
plt.close()

#TODO investigate shape

# ## save progress
# adata.write_h5ad(data_dir / "qc+normalised+hvg.h5ad")
# ## retrieve progress
# adata = sc.read_h5ad(data_dir / "qc+normalised+hvg.h5ad")
# print(f"Loaded data with shape {adata.shape}")

# plot top 20 genes for sanity check
sc.pl.highest_expr_genes(adata, n_top=20, show=False)
plt.savefig(figures_dir / "top20_expr_genes.png")
plt.close()

# endregion

# region 4. Dimension reduction (=Embedding) 

## PCA (linear reduction)
n_pcs = 30
sc.pp.scale(adata, max_value=10)  # scale to variance
sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", random_state=random_seed)

## examine PCA - plot variance ratio
sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs, log=True, show=False)
plt.savefig(figures_dir / "pca_variance_ratio.png")
plt.close()

## examine PCA - plot scatter with different colours
pca_color_cols = ["cell_subtype_3", "disease_type", "meta_tissue_type", "study_name", "total_counts"]
sc.pl.pca_scatter(adata, color=pca_color_cols, ncols=2, show=False)
plt.savefig(figures_dir / "pca_scatter.png", dpi=100)
plt.close()
 
# endregion

# region 5. Clustering
## create k neighbour graph
n_neighbors = 15
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_seed)
sc.tl.umap(adata, min_dist=0.3, spread=1.0, random_state=random_seed)

## detect clusters at one resolution (test pipeline default)
leiden_resolution = 0.5
leiden_key = "leiden_res_0.5"
sc.tl.leiden(adata, resolution=leiden_resolution, key_added=leiden_key)

## save umap with different colourings
umap_color_cols = [leiden_key, "cell_subtype_3", "disease_type", "meta_tissue_type"]
sc.pl.umap(adata, color=umap_color_cols, ncols=2, wspace=0.35, show=False)
plt.savefig(figures_dir / "umap_clustering_overview.png")
plt.close()

## save processed data
adata.write_h5ad(data_dir / "adata_processed.h5ad")

# endregion
