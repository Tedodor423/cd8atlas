# test_pipeline

- Timestamp: `2026-03-20 04:49:42`
- Source file: `/ceph/project/sharmalab/dnimrich/cd8atlas/code/test_pipeline.py`


*Loading from [../../data/subsampled_100000.h5ad](../../data/subsampled_100000.h5ad)*

Loaded subsampled data with shape (100000, 19957)

---
## Quality Control


### Ambient RNA QC


Checking if a cleanup is necessary, flagging cells with less than 1000 counts as potential ambient RNA contamination

No significant RNA contamination detected, skipping cleanup

### Count QC


| index | count | mean | std | min | 25% | 50% | 75% | max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| total_counts | 100000.00 | 3924.63 | 2698.54 | 248.00 | 2350.00 | 3424.00 | 4779.00 | 65975.00 |
| n_genes_by_counts | 100000.00 | 1434.59 | 638.66 | 201.00 | 1009.00 | 1346.00 | 1730.00 | 6872.00 |
| pct_counts_mt | 100000.00 | 4.61 | 2.41 | 0.00 | 3.02 | 4.26 | 5.73 | 47.37 |
| pct_counts_ribo | 100000.00 | 28.53 | 9.60 | 0.68 | 21.57 | 27.68 | 35.02 | 64.87 |
| pct_counts_hb | 100000.00 | 0.01 | 0.23 | 0.00 | 0.00 | 0.00 | 0.00 | 49.51 |
| pct_counts_in_t... | 100000.00 | 20.88 | 3.43 | 8.44 | 18.47 | 20.55 | 22.95 | 59.58 |


count outliers with MAD:1548

mitochondrial outliers with MAD:8020

[<img src="test_pipeline_20260320_044942_figures/qc_before_total_counts_hist.png" width="23%">](test_pipeline_20260320_044942_figures/qc_before_total_counts_hist.png)
[<img src="test_pipeline_20260320_044942_figures/qc_before_pct_counts_mt_hist.png" width="23%">](test_pipeline_20260320_044942_figures/qc_before_pct_counts_mt_hist.png)
[<img src="test_pipeline_20260320_044942_figures/qc_before_pct_counts_mt_violin.png" width="23%">](test_pipeline_20260320_044942_figures/qc_before_pct_counts_mt_violin.png)
[<img src="test_pipeline_20260320_044942_figures/qc_before_scatter_genes_total_mt.png" width="23%">](test_pipeline_20260320_044942_figures/qc_before_scatter_genes_total_mt.png)


Size of matrix before count filtering: (100000, 19957)

Size of matrix after count filtering: (90828, 19832) (removed 9.17% cells)

Filtering out genes expressed in less than 20 cells

Removed 14250 cells (28.15%)

#### After count QC:


[<img src="test_pipeline_20260320_044942_figures/qc_after_total_counts_hist.png" width="23%">](test_pipeline_20260320_044942_figures/qc_after_total_counts_hist.png)
[<img src="test_pipeline_20260320_044942_figures/qc_after_pct_counts_mt_hist.png" width="23%">](test_pipeline_20260320_044942_figures/qc_after_pct_counts_mt_hist.png)
[<img src="test_pipeline_20260320_044942_figures/qc_after_pct_counts_mt_violin.png" width="23%">](test_pipeline_20260320_044942_figures/qc_after_pct_counts_mt_violin.png)
[<img src="test_pipeline_20260320_044942_figures/qc_after_scatter_genes_total_mt.png" width="23%">](test_pipeline_20260320_044942_figures/qc_after_scatter_genes_total_mt.png)


---
## Normalisation


Normalised in reference to total counts

Transformed with log1p

[<img src="test_pipeline_20260320_044942_figures/normalisation_total_counts_hist.png" width="48%">](test_pipeline_20260320_044942_figures/normalisation_total_counts_hist.png)
[<img src="test_pipeline_20260320_044942_figures/normalisation_log1p_hist.png" width="48%">](test_pipeline_20260320_044942_figures/normalisation_log1p_hist.png)


*Saving into [test_pipeline_20260320_044942_data/qc+normalised.h5ad](test_pipeline_20260320_044942_data/qc+normalised.h5ad)*

Saved adata of shape (90828, 14250)

---
## Feature selection


### Gene exclusion


Excluding genes matching the following patterns:

|  |  |
| --- | --- |
| stress_genes | ['G0S2', '^JUN'... |
| cc_genes | ['CDC20', 'CCNB... |
| IFN_genes | ['^IRF', '^IFN'] |
| ccl_genes | ['CCL'] |
| MHC_genes | ['^HLA-[A-Z]'] |
| hist_genes | ['HIST'] |
| comp_genes | ['^C1Q'] |
| ig_genes | ['^IGH', '^IGK'... |
| hb_genes | ['^HBA\|^HBB'] |
| miscelaneous_ba... | [] |


Removing 285 genes (2.00%)

Shape after exclusion panel filtering: (90828, 13965)

### Highly Variable Gene selection


HVG selected: 4000

[<img src="test_pipeline_20260320_044942_figures/hvg_before.png" width="48%">](test_pipeline_20260320_044942_figures/hvg_before.png)
[<img src="test_pipeline_20260320_044942_figures/hvg_after.png" width="48%">](test_pipeline_20260320_044942_figures/hvg_after.png)


### Top 20 expressed genes after selection:


[<img src="test_pipeline_20260320_044942_figures/top20_expr_genes.png" width="50%">](test_pipeline_20260320_044942_figures/top20_expr_genes.png)


---
## Dimension reduction


### Principal component analysis (PCA)


Scaled data with max variance cutoff 10

Calculated PCA with 30 components

[<img src="test_pipeline_20260320_044942_figures/pca_variance_ratio.png" width="48%">](test_pipeline_20260320_044942_figures/pca_variance_ratio.png)
[<img src="test_pipeline_20260320_044942_figures/pca_scatter.png" width="48%">](test_pipeline_20260320_044942_figures/pca_scatter.png)


---
## Clustering


Calculated nearest 15 neighbours

Calculated UMAP with min_dist 0.3 and spread 1.0

Detected clusters with leiden at resolution 0.5

[<img src="test_pipeline_20260320_044942_figures/umap_clustering_overview.png" width="50%">](test_pipeline_20260320_044942_figures/umap_clustering_overview.png)


*Saving into [test_pipeline_20260320_044942_data/adata_processed.h5ad](test_pipeline_20260320_044942_data/adata_processed.h5ad)*

Saved processed adata with shape (90828, 4000)
