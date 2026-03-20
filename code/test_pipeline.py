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


OUTPUT_NAME = "test_pipeline"
root_dir = Path("/ceph/project/sharmalab/dnimrich/cd8atlas")
bioreport.start(OUTPUT_NAME, path=root_dir / "output", automatic_titles=True, git_image_links='Tedodor423/cd8atlas/refs/heads/main')

random_seed = 423

# region 0. load/subsample data
input_h5ad = root_dir / "data/" / "huARdb_v2_GEX.CD8.all_genes.h5ad" #"huARdb_v2_GEX.CD8.hvg4k.h5ad"
adata = sc.read_h5ad(bioreport.file_load(input_h5ad))
bioreport.log(f"Loaded adata with with shape {adata.shape}")

## subsampling

subsample_amount = 100_000
if subsample_amount and adata.n_obs > subsample_amount:
    idx = np.random.default_rng(random_seed).choice(adata.n_obs, size=subsample_amount, replace=False)
    adata = adata[idx, :]
    bioreport.log(f"Subsampled into: {adata.shape}")
# endregion

## save subsampled data
adata.write_h5ad(bioreport.file_save(f"subsampled_{subsample_amount}.h5ad"))

## retrieve subsampled data
adata = sc.read_h5ad(bioreport.file_load(root_dir / "output" / OUTPUT_NAME / "data" / f"subsampled_100000.h5ad"))
bioreport.log(f"Loaded subsampled data with shape {adata.shape}")

## preserve original counts
adata.layers["counts"] = adata.X.copy() 


# region 1. QC
bioreport.log("Quality Control", heading=2)

bioreport.log("Ambient RNA QC", heading=3)
# region 1.1 ambient RNA QC
### check if necessary

low_count_treshold = 1_000
frac_treshold = 0.05
bioreport.log(f"Checking if a cleanup is necessary, flagging cells with less than {low_count_treshold} counts as potential ambient RNA contamination")

low_count_cells = np.asarray(adata.X.sum(axis=1)) < low_count_treshold
if low_count_cells.mean() > frac_treshold:
    bioreport.log("Potential RNA contamination detected, cleaning up the data")
    bioreport.log("no cleanup code exists yet, exiting for now")
    exit()
else:
    bioreport.log("No significant RNA contamination detected, skipping cleanup")

# endregion

# region 1.2 count QC
bioreport.log("Count QC", heading=3)
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
bioreport.table(adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb", "pct_counts_in_top_20_genes"]].describe().T)

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
bioreport.log("count outliers with MAD:" + str(adata.obs["outlier"].value_counts()[True]))

adata.obs["mt_outlier"] = outlier_MAD(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)
bioreport.log("mitochondrial outliers with MAD:" + str(adata.obs["mt_outlier"].value_counts()[True]))

### plot before QC
sns.displot(adata.obs["total_counts"], bins=100, kde=True)
bioreport.figure("qc_before_total_counts_hist")

sns.displot(adata.obs["pct_counts_mt"], bins=100, kde=True)
bioreport.figure("qc_before_pct_counts_mt_hist", sameline=True)

sc.pl.violin(adata, "pct_counts_mt", show=False)
bioreport.figure("qc_before_pct_counts_mt_violin", sameline=True)

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
bioreport.figure("qc_before_scatter_genes_total_mt", sameline=True)

### Filter based on counts
bioreport.log(f"Size of matrix before count filtering: {adata.shape}")
before_filter_counts = adata.shape[0]  # statistics
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)]

### Filter out mt, ribo, hb genes
adata = adata[:, ~(adata.var.mt | adata.var.ribo | adata.var.hb)]
bioreport.log(f"Size of matrix after count filtering: {adata.shape} (removed {1 - adata.shape[0] / before_filter_counts:.2%} cells)")

### Filter genes with too few detected cells before HVG selection
genes_before_filter = adata.n_vars # statistics
min_cells_per_gene = 20
sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
bioreport.log(f"Filtering out genes expressed in less than {min_cells_per_gene} cells: {adata.n_vars} (removed {1 - adata.n_vars / genes_before_filter:.2%})")

### plot after count QC
bioreport.log("After count QC:", heading=4)
sns.displot(adata.obs["total_counts"], bins=100, kde=True)
bioreport.figure("qc_after_total_counts_hist", sameline=True)

sns.displot(adata.obs["pct_counts_mt"], bins=100, kde=True)
bioreport.figure("qc_after_pct_counts_mt_hist", sameline=True)

sc.pl.violin(adata, "pct_counts_mt", show=False)
bioreport.figure("qc_after_pct_counts_mt_violin", sameline=True)

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
bioreport.figure("qc_after_scatter_genes_total_mt", sameline=True)

# endregion

# region 1.3 Doublet QC
bioreport.log("Doublet QC", heading=3)
n_obs_before_doublet_qc = adata.n_obs  # statistics
### prepare input
counts_csr = adata.layers["counts"].tocsr() if not sp.issparse(adata.layers["counts"]) else adata.layers["counts"].tocsr()
expected_doublet_rate=0.05
scrublet = scr.Scrublet(counts_matrix=counts_csr, expected_doublet_rate=expected_doublet_rate, random_state=random_seed)

### run scrublet
bioreport.log(f"Running scrublet with expected doublet rate {expected_doublet_rate}")
doublet_scores, predicted_doublets = scrublet.scrub_doublets()

### save scores, calls, histogram
adata.obs["doublet_score"] = doublet_scores
adata.obs["doublet_call"] = pd.Categorical(predicted_doublets)
adata.uns["scrublet_sim_scores"] = scrublet.doublet_scores_sim_

scrublet.plot_histogram()
bioreport.figure("qc_scrublet_histogram")

### Filter if necessary
if adata.obs["doublet_call"].nunique() > 1:
    bioreport.log(f"There are {adata.obs['doublet_call'].nunique()} predicted doublets, applying filter")
    adata = adata[~adata.obs['doublet_call'].astype(bool).to_numpy(), :]
else:
    bioreport.log("No doublets predicted")
bioreport.log(f"After doublet step: {adata.shape} (removed {1 - adata.n_obs / n_obs_before_doublet_qc :.2%} cells)")

# endregion
# endregion

# region 2. Normalisation

# ## retrieve count values - not necessary now
# adata.X = adata.layers["counts"].copy() 
bioreport.log("Normalisation", heading=2)

sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
adata.layers["normalised"] = adata.X.copy()
bioreport.log("Normalised in reference to total counts")

sc.pp.log1p(adata)
adata.X = adata.X.astype(np.float32)  # for better storage
adata.layers["log1p"] = adata.X.copy()
bioreport.log("Transformed with log1p")

log1p_library_size = np.asarray(adata.layers["log1p"].sum(axis=1)).ravel()

sns.histplot(adata.obs["total_counts"], bins=100, kde=False)
bioreport.figure("normalisation_total_counts_hist")
sns.histplot(log1p_library_size, bins=100, kde=False)
bioreport.figure("normalisation_log1p_hist", sameline=True)


# endregion

## save progress
adata.write_h5ad(bioreport.file_save("qc+normalised.h5ad"))
bioreport.log(f"Saved adata of shape {adata.shape}")

# ## retrieve progress
# adata = sc.read_h5ad(data_dir / "qc+normalised.h5ad")
# print(f"Loaded data with shape {adata.shape}")


# labeling, CC analysis

# region 3. Feature selection
bioreport.log("Feature selection", heading=2)
## 3.1 Cell selection - none here

## 3.2 Gene exclusion
bioreport.log("Gene exclusion", heading=3)

excluded_gene_sets = {
    "stress_genes": [
        "G0S2", "JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "CDKN1A"
    ],
    "cc_genes": [
        "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "TYMS",
        "FEN1", "MKI67", "TOP2A", "TK1", "RRM1", "RRM2", "HELLS", "UNG",
        "GINS2", "CDC6", "CDK1", "CDC20", "CCNB1", "CCNB2", "AURKA", "AURKB",
        "BUB1", "BUB1B", "UBE2C", "CENPF", "CENPE", "HMGB2",
    ]
}
bioreport.log("Excluding the following genes:")
bioreport.table(excluded_gene_sets)

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

bioreport.log("Genes matched by exclusion panel:")
bioreport.table(excluded_gene_sets)

bioreport.log(f"Removing {int(excluded_gene_mask.sum())} genes from exclusion panel before HVG selection")
adata.var["excluded_feature"] = excluded_gene_mask.to_numpy()
adata = adata[:, ~adata.var["excluded_feature"]].copy()
bioreport.log(f"Shape after exclusion panel filtering: {adata.shape}")



## 3.3. HVG selection
bioreport.log("Highly Variable Gene selection", heading=3)
### using seurat as seurat_v3 with original counts raises an environment error
sc.pp.highly_variable_genes(adata, layer='log1p', flavor='seurat', n_top_genes=4000, inplace=True)
bioreport.log(f"HVG selected: {int(adata.var['highly_variable'].sum())}")

#sc.pl.highly_variable_genes(adata, show=False)
sns.scatterplot(data=adata.var, x="means", y="dispersions", hue="highly_variable", s=5)
bioreport.figure("hvg_before")

### Filter
adata = adata[:, adata.var.highly_variable]
 
#sc.pl.highly_variable_genes(adata, show=False)
sns.scatterplot(data=adata.var, x="means", y="dispersions", color="orange", s=5)
bioreport.figure("hvg_after", sameline=True)

#TODO investigate shape

# ## save progress
# adata.write_h5ad(data_dir / "qc+normalised+hvg.h5ad")
# ## retrieve progress
# adata = sc.read_h5ad(data_dir / "qc+normalised+hvg.h5ad")
# print(f"Loaded data with shape {adata.shape}")

# plot top 20 genes for sanity check
bioreport.log("Top 20 expressed genes after selection:", heading=3)
sc.pl.highest_expr_genes(adata, n_top=20, show=False)
bioreport.figure("top20_expr_genes")

# endregion

# region 4. Dimension reduction (=Embedding) 
bioreport.log("Dimension reduction", heading=2)

## PCA (linear reduction)
bioreport.log("Principal component analysis (PCA)", heading=3)
n_pcs = 30
var_scale_cutoff = 10
sc.pp.scale(adata, max_value=var_scale_cutoff)  # scale to variance
bioreport.log(f"Scaled data with max variance cutoff {var_scale_cutoff}")

sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", random_state=random_seed)
bioreport.log(f"Calculated PCA with {n_pcs} components")

## examine PCA - plot variance ratio
sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs, log=True, show=False)
bioreport.figure("pca_variance_ratio")

## examine PCA - plot scatter with different colours
pca_color_cols = ["cell_subtype_3", "disease_type", "meta_tissue_type", "study_name", "total_counts"]
sc.pl.pca_scatter(adata, color=pca_color_cols, ncols=2, show=False)
bioreport.figure("pca_scatter", dpi=100, sameline=True)
 
# endregion

# region 5. Clustering
bioreport.log("Clustering", heading=2)
## create k neighbour graph
n_neighbors = 15
min_dist = 0.3
spread = 1.0
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_seed)
bioreport.log(f"Calculated nearest {n_neighbors} neighbours")
sc.tl.umap(adata, min_dist=min_dist, spread=spread, random_state=random_seed)
bioreport.log(f"Calculated UMAP with min_dist {min_dist} and spread {spread}")

## detect clusters at one resolution (test pipeline default)
leiden_resolution = 0.5
leiden_key = "leiden_res_0.5"
sc.tl.leiden(adata, resolution=leiden_resolution, key_added=leiden_key)
bioreport.log(f"Detected clusters with leiden at resolution {leiden_resolution}")

## save umap with different colourings
umap_color_cols = [leiden_key, "cell_subtype_3", "disease_type", "meta_tissue_type"]
sc.pl.umap(adata, color=umap_color_cols, ncols=2, wspace=0.35, show=False)
bioreport.figure("umap_clustering_overview")

## save processed data
adata.write_h5ad(bioreport.file_save("adata_processed.h5ad"))
bioreport.log(f"Saved processed adata with shape {adata.shape}")

# endregion


bioreport.exportpdf()