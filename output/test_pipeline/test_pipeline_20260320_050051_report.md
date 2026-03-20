# test_pipeline

- Timestamp: `2026-03-20 05:00:51`
- Source file: `/ceph/project/sharmalab/dnimrich/cd8atlas/code/test_pipeline.py`


*Loading from [../../data/huARdb_v2_GEX.CD8.all_genes.h5ad](../../data/huARdb_v2_GEX.CD8.all_genes.h5ad)*

Loaded adata with with shape (1125712, 19957)

Subsampled into: (100000, 19957)

*Saving into [../../data/subsampled_100000.h5ad](../../data/subsampled_100000.h5ad)*

---
## 1. Quality Control


### 1.1 Ambient RNA QC


Checking if a cleanup is necessary, flagging cells with less than 1000 counts as potential ambient RNA contamination

No significant RNA contamination detected, skipping cleanup

### 1.2 Count QC


| index | count | mean | std | min | 25% | 50% | 75% | max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| total_counts | 1125712.00 | 3925.53 | 2724.44 | 214.00 | 2348.00 | 3420.00 | 4764.00 | 90450.00 |
| n_genes_by_counts | 1125712.00 | 1433.83 | 641.46 | 163.00 | 1010.00 | 1344.00 | 1727.00 | 8466.00 |
| pct_counts_mt | 1125712.00 | 4.61 | 2.42 | 0.00 | 3.00 | 4.26 | 5.74 | 61.19 |
| pct_counts_ribo | 1125712.00 | 28.52 | 9.63 | 0.46 | 21.52 | 27.63 | 35.04 | 69.67 |
| pct_counts_hb | 1125712.00 | 0.01 | 0.16 | 0.00 | 0.00 | 0.00 | 0.00 | 49.51 |
| pct_counts_in_top_20_genes | 1125712.00 | 20.88 | 3.43 | 5.70 | 18.48 | 20.55 | 22.94 | 65.92 |


count outliers with `MAD > 5`:17798

mitochondrial outliers with MAD > 5:89879

[<img src="test_pipeline_20260320_050051_figures/qc_before_total_counts_hist.png" width="23%">](test_pipeline_20260320_050051_figures/qc_before_total_counts_hist.png)
[<img src="test_pipeline_20260320_050051_figures/qc_before_pct_counts_mt_hist.png" width="23%">](test_pipeline_20260320_050051_figures/qc_before_pct_counts_mt_hist.png)
[<img src="test_pipeline_20260320_050051_figures/qc_before_pct_counts_mt_violin.png" width="23%">](test_pipeline_20260320_050051_figures/qc_before_pct_counts_mt_violin.png)
[<img src="test_pipeline_20260320_050051_figures/qc_before_scatter_genes_total_mt.png" width="23%">](test_pipeline_20260320_050051_figures/qc_before_scatter_genes_total_mt.png)


Size of matrix before count filtering: (1125712, 19957)

Size of matrix after count filtering: (1022355, 19832) (removed 9.18% cells)

Filtering out genes expressed in less than 20 cells

Removed 16849 cells (15.04%)

#### After count QC:


[<img src="test_pipeline_20260320_050051_figures/qc_after_total_counts_hist.png" width="23%">](test_pipeline_20260320_050051_figures/qc_after_total_counts_hist.png)
[<img src="test_pipeline_20260320_050051_figures/qc_after_pct_counts_mt_hist.png" width="23%">](test_pipeline_20260320_050051_figures/qc_after_pct_counts_mt_hist.png)
[<img src="test_pipeline_20260320_050051_figures/qc_after_pct_counts_mt_violin.png" width="23%">](test_pipeline_20260320_050051_figures/qc_after_pct_counts_mt_violin.png)
[<img src="test_pipeline_20260320_050051_figures/qc_after_scatter_genes_total_mt.png" width="23%">](test_pipeline_20260320_050051_figures/qc_after_scatter_genes_total_mt.png)


---
## 2. Normalisation


Normalised in reference to total counts

Transformed with log1p

[<img src="test_pipeline_20260320_050051_figures/normalisation_total_counts_hist.png" width="48%">](test_pipeline_20260320_050051_figures/normalisation_total_counts_hist.png)
[<img src="test_pipeline_20260320_050051_figures/normalisation_log1p_hist.png" width="48%">](test_pipeline_20260320_050051_figures/normalisation_log1p_hist.png)


*Saving into [test_pipeline_20260320_050051_data/qc+normalised.h5ad](test_pipeline_20260320_050051_data/qc+normalised.h5ad)*

Saved adata of shape (1022355, 16849)

---
## 3. Feature selection


### 3.2 Gene exclusion


Excluding genes matching the following patterns:

|  |  |
| --- | --- |
| stress_genes | ['G0S2', '^JUN', '^FOS', 'CDKN1A'] |
| cc_genes | ['CDC20', 'CCNB2', 'DLGAP5', 'PLK... |
| IFN_genes | ['^IRF', '^IFN'] |
| ccl_genes | ['CCL'] |
| MHC_genes | ['^HLA-[A-Z]'] |
| hist_genes | ['HIST'] |
| comp_genes | ['^C1Q'] |
| ig_genes | ['^IGH', '^IGK', '^IGL'] |
| hb_genes | ['^HBA\|^HBB'] |
| miscelaneous_bad_features | [] |


Removing 312 genes (1.85%)

Shape after exclusion panel filtering: (1022355, 16537)

### 3.3 Highly Variable Gene selection


HVG selected: 4000

[<img src="test_pipeline_20260320_050051_figures/hvg_before.png" width="48%">](test_pipeline_20260320_050051_figures/hvg_before.png)
[<img src="test_pipeline_20260320_050051_figures/hvg_after.png" width="48%">](test_pipeline_20260320_050051_figures/hvg_after.png)


### Top 20 expressed genes after selection:


[<img src="test_pipeline_20260320_050051_figures/top20_expr_genes.png" width="50%">](test_pipeline_20260320_050051_figures/top20_expr_genes.png)


---
## 4. Dimension reduction


### Principal component analysis (PCA)


Scaled data with max variance cutoff 10

Calculated PCA with 30 components

[<img src="test_pipeline_20260320_050051_figures/pca_variance_ratio.png" width="48%">](test_pipeline_20260320_050051_figures/pca_variance_ratio.png)
[<img src="test_pipeline_20260320_050051_figures/pca_scatter.png" width="48%">](test_pipeline_20260320_050051_figures/pca_scatter.png)


---
## 5. Clustering


Calculated nearest 15 neighbours

Calculated UMAP with min_dist 0.3 and spread 1.0

Detected clusters with leiden at resolution 0.5

[<img src="test_pipeline_20260320_050051_figures/umap_clustering_overview.png" width="50%">](test_pipeline_20260320_050051_figures/umap_clustering_overview.png)


*Saving into [test_pipeline_20260320_050051_data/adata_processed.h5ad](test_pipeline_20260320_050051_data/adata_processed.h5ad)*

Saved processed adata with shape (1022355, 4000)
