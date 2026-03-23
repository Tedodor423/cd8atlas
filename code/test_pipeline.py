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

# # region 0. load/subsample data
# input_h5ad = root_dir / "data/" / "huARdb_v2_GEX.CD8.all_genes.h5ad" #"huARdb_v2_GEX.CD8.hvg4k.h5ad"
# adata = sc.read_h5ad(bioreport.file_load(input_h5ad))
# bioreport.log(f"Loaded adata with with shape {adata.shape}")

# # subsample into adata_subsampled

# subsample_amount = 100_000
# adata_subsampled = adata
# if subsample_amount and adata.n_obs > subsample_amount:
#     idx = np.random.default_rng(random_seed).choice(adata.n_obs, size=subsample_amount, replace=False)
#     adata_subsampled = adata[idx, :].copy()
#     bioreport.log(f"Subsampled into: {adata_subsampled.shape}")
# # endregion

# # save subsampled data
# adata_subsampled.write_h5ad(bioreport.file_save(root_dir / "data/" / f"subsampled_{subsample_amount}.h5ad"))

## retrieve subsampled data
adata = sc.read_h5ad(bioreport.file_load(root_dir / "data" / f"subsampled_100000.h5ad"))
bioreport.log(f"Loaded subsampled data with shape {adata.shape}")

## preserve original counts
adata.layers["counts"] = adata.X.copy() 


# region 1. QC
bioreport.log("1. Quality Control", heading=2)

bioreport.log("1.1 Ambient RNA QC", heading=3)
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
bioreport.log("1.2 Count QC", heading=3)
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

outlier_var = 5
adata.obs["outlier"] = (
    outlier_MAD(adata, "log1p_total_counts", outlier_var)
    | outlier_MAD(adata, "log1p_n_genes_by_counts", outlier_var)
    | outlier_MAD(adata, "pct_counts_in_top_20_genes", outlier_var)
)
bioreport.log(f"Number of count outliers with `MAD > {outlier_var}`: " + str(adata.obs["outlier"].value_counts()[True]))

adata.obs["mt_outlier"] = outlier_MAD(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)
bioreport.log(f"Number of mitochondrial outliers with `MAD > {outlier_var}`: " + str(adata.obs["mt_outlier"].value_counts()[True]))

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
bioreport.log(f"Filtering out genes expressed in less than {min_cells_per_gene} cells")
sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
bioreport.log(f"Removed {adata.n_vars} cells ({1 - adata.n_vars / genes_before_filter:.2%})")

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
bioreport.log("1.3 Doublet QC", heading=3)
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


# ## retrieve count values - not necessary now
# adata.X = adata.layers["counts"].copy() 

# region 2. Normalisation
bioreport.log("2. Normalisation", heading=2)

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

## retrieve progress
adata = sc.read_h5ad(bioreport.file_load("qc+normalised.h5ad"))
bioreport.log(f"Loaded data with shape {adata.shape}")


# labeling, CC analysis

# region 3. Feature selection
bioreport.log("3. Feature selection", heading=2)
## 3.1 Cell selection - none here

## 3.2 Gene exclusion
bioreport.log("3.2 Gene exclusion", heading=3)

excluded_gene_patterns = {
    "stress_genes": ["G0S2", r"^JUN", r"^FOS", "CDKN1A"],
    #"cc_genes": [r"^MCM", r"^CDC", r"^CDK"],
    "cc_genes": [
        # Source: https://www.nature.com/articles/s41592-025-02793-1#Sec30 
        # https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-025-02793-1/MediaObjects/41592_2025_2793_MOESM5_ESM.xlsx
        # - top 50 genes from the gene expression programmes CellCycle-G2M, CellCycle-S, CellCycle-Late-S
        "CDC20","CCNB2","DLGAP5","PLK1","UBE2C","CCNB1","HMMR","BIRC5","CENPF","ASPM","TPX2","TOP2A","TROAP","KIF14","CENPA","CDKN3","MKI67","KIF20A","CEP55","CENPE","HMGB2","GTSE1","PTTG1","NEK2","PIMREG","KIF2C","AURKB","FAM64A","CDCA8","CDCA3","KIF4A","KIF23","DEPDC1","NUF2","AURKA","CCNA2","BUB1","DEPDC1B","CDCA2","STMN1","CKAP2L", #"CDK1","NUSAP1","TUBA1B","HMGB3","TTK","KNL1","KIF15","HJURP","KIFC1","KPNA2","PRR11","NCAPG","CKS2","CKS1B","CENPW","TUBB","CDC25C","KIF11","HMGN2","TUBA1C","SGO2","SMC4","FOXM1","CIT","ARL6IP1","SHCBP1","SGO1","TUBB4B","PBK","PSRC1","SAPCD2","H2AFZ","KNSTRN","HN1","UBE2S","H4C3","PRC1","JPT1","KIF18B","MXD3","SPC25","SPC24","ANLN","RACGAP1","ECT2","ARHGAP11A","AB_CD71","BUB1B","MAD2L1","PIF1","SPAG5","KIAA1524","NCAPH","H2AZ1","ANP32E","H2AFV","CCNF","SMC2","PCLAF","KIF20B","CENPN","KIF18A","SKA3","NDC80","PARPBP","POC1A","FAM72B","RAD21","HMGB1","FAM72C","ESPL1","DIAPH3","PHF19","OIP5","APOBEC3B","SKA1","SPDL1","INCENP","ZWINT","NUCKS1","MTFR2","FAM72D","NCAPD2","CIP2A","CKAP5","NEIL3","FBXO43","CKAP2","KIF22","MELK","CENPM","TRIP13","TMPO","GPSM2","H2AFX","SKA2","ARHGEF39","CALM3","DTYMK","LMNB1","UBE2T","GAS2L3","BCL2L12","STIL","TYMS","TK1","H2AX","PLK4","RPL39L","SMTN","FAM83D","CDC25B","DCTN3","POLQ","LSM5","GGH","HNRNPA2B1","CALM2","HDGF","TICRR","AB_CD61","TRAIP","AB_CD278","TACC3","RAN","MND1","ANP32B","AB_CD9","H2AC13","SPTBN1","BORA","MZT1","RDM1","C14orf80","H3C8","CDCA5","NCAPG2","ASF1B","TBXT","GSG2","HP1BP3","LBR","GINS1","MYBL2","BRCA2","EMC9","ARHGAP11B","DEK","GEN1","DDX39A","TPM4","H2AC14","FANCI","UBALD2","BUB3","NUDT1","REEP4","LGALS1","ERCC6L",
        "GINS2","TK1","TYMS","DTL","UHRF1","PCLAF","CLSPN","CDC45","MCM10","FAM111B","CDT1","CDC6","E2F1","PCNA","RRM2","MCM4","ASF1B","MCM2","CENPU","STMN1","MCM7","DUT","PKMYT1","MYBL2","MCM5","ZWINT","CHEK1","DHFR","HELLS","FEN1","MCM3","BRCA1","CDCA7","EXO1","MCM6","ORC6","PAQR4","TUBA1B","C16orf59","RAD51AP1","DSCC1", #"ORC1","CDCA5","RAD51","CENPM","ESCO2","MND1","CHAF1A","TCF19","RAD54L","AB_CD71","CCNE2","RMI2","E2F8","TEDC2","WDR34","NASP","CDC25A","GMNN","ATAD2","E2F2","SMC2","DMC1","CENPW","RBBP8","DEK","ATAD5","TMEM106C","TIMELESS","H2AFZ","CTNNAL1","POLE2","WDR76","H2AZ1","FANCI","SLBP","BRCA2","HMGB2","GGH","RFC4","RRM1","DIAPH3","WDHD1","LIG1","CENPH","DNAJC9","TUBB","DNMT1","HMGN2","GINS1","E2F7","UBE2T","PTMA","UNG","CHAF1B","SGO1","BRIP1","AUNIP","TYMSOS","RPA3","MELK","MKI67","FBXO5","CARHSP1","NUDT1","CCNE1","SPC25","SHCBP1","POLQ","LMNB1","XRCC2","RFC2","MSH6","SKA3","CENPN","HMGB1","RIBC2","PSMC3IP","EZH2","CDCA7L","SNRNP25","C19orf48","POC1A","CDK1","POLD3","RANBP1","CENPK","MAD2L1","RECQL4","SIVA1","USP1","ACOT7","RFC5","ZNF367","DNA2","MTHFD1","GINS3","IFI27L1","H2AFV","GAPDH","IDH2","H4C3","HNRNPAB","SMC3","DSN1","NCAPH","NCAPG2","RFC3","CENPP","CDK2","CENPX","RNASEH2A","BIRC5","SPC24","NEIL3","H2AFY","ZGRF1","KIAA1524","OIP5","ANP32B","CEP55","CDK4","FANCG","SLC29A1","GINS4","KNTC1","SLC43A3","SLC25A5","PXMP2","GGCT","TFDP1","POLA1","HNRNPF","FABP5","CBX5","POLD1","MTFR2","KIFC1","TP73","SMC1A","FANCA","PAICS","PPIA","NSD2","C21orf58","PHGDH","STIL","RAN","SNRPB","TOPBP1","PLK4","DTYMK","STRIP2","FOXM1","BARD1","HMGA1","EME1","MSH2","SYCE2","MMS22L",
        "HIST1H1B","HIST1H4C","HIST1H3G","HIST1H3B","RRM2","SPC25","CDCA5","NUSAP1","HIST1H2AJ","HJURP","TOP2A","HIST1H3F","TYMS","KIF15","PKMYT1","ESCO2","CDK1","E2F8","KIFC1","UBE2C","DIAPH3","CKAP2L","MKI67","KIF23","TUBB","NCAPG","HIST1H2AH","PBK","GTSE1","AURKB","H4C3","ANLN","MELK","KNL1","CDCA2","C2orf48","STMN1","KIF11","HIST1H2AL","CDCA3","SPC24", #"TUBA1B","MND1","CCNA2","RAD51AP1","KIF2C","PCLAF","SKA3","FBXO5","TK1","ZWINT","HIST1H3D","CDCA8","SKA1","BUB1B","POLQ","KIF4A","H2AFX","SHCBP1","KIF18B","ASF1B","ASPM","HIST1H2AI","MYBL2","TPX2","HIST1H1C","ESPL1","HIST1H1D","MXD3","HIST1H3C","E2F7","CDC25C","HIST1H2AM","NCAPH","FOXM1","FANCI","FAM111B","BIRC5","HIST2H2AC","RAD54L","UBE2T","FBXO43","HIST1H2BH","STIL","NEIL3","CENPM","DLGAP5","EZH2","PCNA","HIST1H1E","HMMR","CLSPN","CENPF","NDC80","CKS1B","SAC3D1","CDT1","CENPU","HMGB1","CKS2","WDR34","TTK","ATAD2","TICRR","SMC4","CENPA","HMGB2","CENPE","PRC1","KIF14","C21orf58","BUB1","AB_CD38-2","DEPDC1","HMGN2","HASPIN","SGO1","RACGAP1","DSCC1","TMPO","GGH","HIST1H4H","TRIP13","ARHGAP11A","HIST1H2BM","HIST1H2AG","HIST2H4B","HIST1H2BO","CDC45","NCAPG2","ECT2","KIF22","H2AFZ","DUT","UHRF1","KIAA1524","EXO1","CENPW","SMC2","HIST1H2BL","SPAG5","POC1A","CIT","CIP2A","TUBB4B","RAD51","E2F2","CCNF","SGO2","RRM1","HIST1H4D","KIF24","BRCA1","TMEM106C","HIST1H2BF","LMNB1","KIF20A","DEPDC1B","CDKN2C","H2AFV","CENPN","MTFR2","APOBEC3B","ERCC6L","BRCA2","CEP55","TACC3","AURKA","HIST1H2AB","CCDC150","FEN1","EME1","GSG2","HIST2H2AB","PLK4","WDR62","NUF2","TCF19","INCENP","CCDC34","HIST1H1A","C1QL1","HIST1H4I","C16orf59","ORC6","H2AZ1","MCM10","DHFR","CDCA4","CARHSP1","TYMSOS","H1-1","HIST1H2BI","H1-5","CDC25A","RECQL4","C17orf53","HIST1H2BB","ARHGEF39","DNAJC9",
    ],
    "IFN_genes": [r"^IRF", r"^IFN"],
    "ccl_genes": [r"CCL"],
    "MHC_genes": [r"^HLA-[A-Z]"],
    "hist_genes": [r"HIST"],
    "comp_genes": [r"^C1Q"],
    "ig_genes": [r"^IGH", r"^IGK", r"^IGL"],
    "hb_genes": [r"^HBA|^HBB"],
    "miscelaneous_bad_features": []
}
bioreport.log("Excluding genes matching the following patterns:")
bioreport.table(excluded_gene_patterns)

### join all the lists together
all_excluded_gene_patterns = [pattern for patterns in excluded_gene_patterns.values() for pattern in patterns]

### create the mask with regex
excluded_gene_mask = pd.Series(False, index=adata.var_names, dtype=bool)
for pattern in all_excluded_gene_patterns:
    excluded_gene_mask |= adata.var_names.str.upper().str.contains(pattern.upper(), regex=True)

### Filter
bioreport.log(f"Removing {int(excluded_gene_mask.sum())} genes ({excluded_gene_mask.sum() / len(adata.var_names) * 100:.2f}%)")
adata = adata[:, ~excluded_gene_mask.to_numpy()].copy()
bioreport.log(f"Shape after exclusion panel filtering: {adata.shape}")



# ## 3.3. HVG selection
# bioreport.log("3.3 Highly Variable Gene selection", heading=3)
# ### using seurat as seurat_v3 with original counts raises an environment error
# # sc.pp.highly_variable_genes(adata, layer='log1p', flavor='seurat', n_top_genes=4000, inplace=True)
# sc.pp.highly_variable_genes(adata, layer='counts', flavor='seurat_v3', n_top_genes=4000, inplace=True)
# bioreport.log(f"HVGs selected: {int(adata.var['highly_variable'].sum())}")

# #sc.pl.highly_variable_genes(adata, show=False)
# sns.scatterplot(data=adata.var, x="means", y="variances_norm", hue="highly_variable", s=5)
# bioreport.figure("hvg_before")

# ### Filter
# adata = adata[:, adata.var.highly_variable]
 
# #sc.pl.highly_variable_genes(adata, show=False)
# sns.scatterplot(data=adata.var, x="means", y="variances_norm", color="orange", s=5)
# bioreport.figure("hvg_after", sameline=True)


## 3.3. HVG selection
bioreport.log("3.3 Highly Variable Gene selection", heading=3)
### using seurat as seurat_v3 with original counts raises an environment error
# sc.pp.highly_variable_genes(adata, layer='log1p', flavor='seurat', n_top_genes=4000, inplace=True)
sc.pp.highly_variable_genes(adata, layer='counts', flavor='seurat_v3', n_top_genes=2000, inplace=True)
bioreport.log(f"HVGs selected: {int(adata.var['highly_variable'].sum())}")

#sc.pl.highly_variable_genes(adata, show=False)
sns.scatterplot(data=adata.var, x="means", y="variances_norm", hue="highly_variable", s=5)
bioreport.figure("hvg_before2")

### Filter
adata = adata[:, adata.var.highly_variable]
 
sc.pl.highly_variable_genes(adata, show=False)
# sns.scatterplot(data=adata.var, x="means", y="variances_norm", color="orange", s=5)
bioreport.figure("hvg_after2", sameline=True)

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

# region 4. Dimensional reduction (=Embedding) 
bioreport.log("4. Dimensional reduction", heading=2)

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
bioreport.log("5. Clustering", heading=2)
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
