# test_cc_genes_clustering

- Timestamp: `2026-03-24 13:41:50`
- Source file: `/ceph/project/sharmalab/dnimrich/cd8atlas/code/pipeline_elements.py`


*Loading from [../../data/subsampled_100000.h5ad](../../data/subsampled_100000.h5ad)*

Loaded adata with with shape (100000, 19957)

---
## 1. Quality Control


### 1.1 Ambient RNA QC


Checking if a cleanup is necessary, flagging cells with less than 1000 counts as potential ambient RNA contamination

No significant RNA contamination detected, skipping cleanup

### 1.2 Count QC


| index | count | mean | std | min | 25% | 50% | 75% | max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| total_counts | 100000.00 | 3924.63 | 2698.54 | 248.00 | 2350.00 | 3424.00 | 4779.00 | 65975.00 |
| n_genes_by_counts | 100000.00 | 1434.59 | 638.66 | 201.00 | 1009.00 | 1346.00 | 1730.00 | 6872.00 |
| pct_counts_mt | 100000.00 | 4.61 | 2.41 | 0.00 | 3.02 | 4.26 | 5.73 | 47.37 |
| pct_counts_ribo | 100000.00 | 28.53 | 9.60 | 0.68 | 21.57 | 27.68 | 35.02 | 64.87 |
| pct_counts_hb | 100000.00 | 0.01 | 0.23 | 0.00 | 0.00 | 0.00 | 0.00 | 49.51 |
| pct_counts_in_top_20_genes | 100000.00 | 20.88 | 3.43 | 8.44 | 18.47 | 20.55 | 22.95 | 59.58 |


[<img src="test_cc_genes_clustering_20260324_134150_figures/qc_before_total_counts_hist.png" width="23%">](test_cc_genes_clustering_20260324_134150_figures/qc_before_total_counts_hist.png)
[<img src="test_cc_genes_clustering_20260324_134150_figures/qc_before_pct_counts_mt_hist.png" width="23%">](test_cc_genes_clustering_20260324_134150_figures/qc_before_pct_counts_mt_hist.png)
[<img src="test_cc_genes_clustering_20260324_134150_figures/qc_before_pct_counts_mt_violin.png" width="23%">](test_cc_genes_clustering_20260324_134150_figures/qc_before_pct_counts_mt_violin.png)
[<img src="test_cc_genes_clustering_20260324_134150_figures/qc_before_scatter_genes_total_mt.png" width="23%">](test_cc_genes_clustering_20260324_134150_figures/qc_before_scatter_genes_total_mt.png)


Number of count outliers with `MAD > 5`: 1548

Number of mitochondrial outliers with `MAD > 5` or `pct_counts_mt > 8`: 7739

Filtered out outliers: 8902 (8.90%)

Filtering out genes expressed in less than 20 cells

Removed 14368 genes (28.01%)

#### After count QC:


[<img src="test_cc_genes_clustering_20260324_134150_figures/qc_after_total_counts_MAD=5.png" width="23%">](test_cc_genes_clustering_20260324_134150_figures/qc_after_total_counts_MAD=5.png)
[<img src="test_cc_genes_clustering_20260324_134150_figures/qc_after_pct_counts_mt_MAD=5.png" width="23%">](test_cc_genes_clustering_20260324_134150_figures/qc_after_pct_counts_mt_MAD=5.png)
[<img src="test_cc_genes_clustering_20260324_134150_figures/qc_after_pct_counts_mt_violin_MAD=5.png" width="23%">](test_cc_genes_clustering_20260324_134150_figures/qc_after_pct_counts_mt_violin_MAD=5.png)
[<img src="test_cc_genes_clustering_20260324_134150_figures/qc_after_scatter_genes_total_mt_MAD=5.png" width="23%">](test_cc_genes_clustering_20260324_134150_figures/qc_after_scatter_genes_total_mt_MAD=5.png)


### 1.4 Gene exclusion


Excluding genes matching the following patterns:

_Note: Showing at most 10 rows and 10 columns from data with 11 rows and 2 columns._

|  |  |
| --- | --- |
| mt_genes | ['MT-'] |
| ribo_genes | ['RPS', 'RPL'] |
| hb_genes | ['^HBA\|^HBB'] |
| stress_genes | ['G0S2', '^JUN', '^FOS', 'CDKN1A'] |
| IFN_genes | ['^IRF', '^IFN'] |
| ccl_genes | ['CCL'] |
| MHC_genes | ['^HLA-[A-Z]'] |
| hist_genes | ['HIST'] |
| comp_genes | ['^C1Q'] |
| ig_genes | ['^IGH', '^IGK', '^IGL'] |


Removing 340 genes (2.37%)

Shape after exclusion panel filtering: (91098, 14028)

---
## 2. Normalisation


Normalised in reference to total counts

Transformed with log1p

[<img src="test_cc_genes_clustering_20260324_134150_figures/normalisation_total_counts_hist.png" width="48%">](test_cc_genes_clustering_20260324_134150_figures/normalisation_total_counts_hist.png)
[<img src="test_cc_genes_clustering_20260324_134150_figures/normalisation_log1p_hist.png" width="48%">](test_cc_genes_clustering_20260324_134150_figures/normalisation_log1p_hist.png)


---
## 3. Feature selection


### 3.3 Highly Variable Gene selection


HVGs selected: 4000

[<img src="test_cc_genes_clustering_20260324_134150_figures/hvg_before_seurat_4000.png" width="48%">](test_cc_genes_clustering_20260324_134150_figures/hvg_before_seurat_4000.png)
[<img src="test_cc_genes_clustering_20260324_134150_figures/hvg_after_seurat_4000.png" width="48%">](test_cc_genes_clustering_20260324_134150_figures/hvg_after_seurat_4000.png)


### Top 20 expressed genes after selection:


[<img src="test_cc_genes_clustering_20260324_134150_figures/top20_expr_genes.png" width="50%">](test_cc_genes_clustering_20260324_134150_figures/top20_expr_genes.png)


CC genes coloured:

[<img src="test_cc_genes_clustering_20260324_134150_figures/cc_genes_hvg.png" width="50%">](test_cc_genes_clustering_20260324_134150_figures/cc_genes_hvg.png)


---
## 4. Dimensional reduction


### Principal component analysis (PCA)


Scaled data with max variance cutoff 10

Calculated PCA with 30 components

[<img src="test_cc_genes_clustering_20260324_134150_figures/pca_variance_ratio_n=30.png" width="48%">](test_cc_genes_clustering_20260324_134150_figures/pca_variance_ratio_n=30.png)
[<img src="test_cc_genes_clustering_20260324_134150_figures/pca_scatter_n=30.png" width="48%">](test_cc_genes_clustering_20260324_134150_figures/pca_scatter_n=30.png)


---
## 5. Clustering


Calculated nearest 15 neighbours using 30 PCs

Calculated UMAP with min_dist 0.3 and spread 1.0

Detected clusters with leiden at resolution 0.5

[<img src="test_cc_genes_clustering_20260324_134150_figures/umap_clustering_overview_n=15_mindist=0.3_res=0.5.png" width="50%">](test_cc_genes_clustering_20260324_134150_figures/umap_clustering_overview_n=15_mindist=0.3_res=0.5.png)


*Saving into [test_cc_genes_clustering_20260324_134150_data/adata_umap_clustering_n=15_mindist=0.3_res=0.5.h5ad](test_cc_genes_clustering_20260324_134150_data/adata_umap_clustering_n=15_mindist=0.3_res=0.5.h5ad)*


> Saved adata with shape (91098, 4000)
