# parameter_iteration

- Timestamp: `2026-04-21 03:53:52`
- Source file: `/ceph/project/sharmalab/dnimrich/cd8atlas/code/pipeline_elements.py`


---
# Loading and preprocessing data


*Loading from [../../data/subsampled_100000.h5ad](../../data/subsampled_100000.h5ad)*

Loaded adata with with shape (100000, 19957)

Initialized `counts` layer from adata.X

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


[<img src="parameter_iteration_20260421_035352_figures/qc_before_total_counts_hist.png" width="23%">](parameter_iteration_20260421_035352_figures/qc_before_total_counts_hist.png)
[<img src="parameter_iteration_20260421_035352_figures/qc_before_pct_counts_mt_hist.png" width="23%">](parameter_iteration_20260421_035352_figures/qc_before_pct_counts_mt_hist.png)
[<img src="parameter_iteration_20260421_035352_figures/qc_before_pct_counts_mt_violin.png" width="23%">](parameter_iteration_20260421_035352_figures/qc_before_pct_counts_mt_violin.png)
[<img src="parameter_iteration_20260421_035352_figures/qc_before_scatter_genes_total_mt.png" width="23%">](parameter_iteration_20260421_035352_figures/qc_before_scatter_genes_total_mt.png)


Number of count outliers with `MAD > 5`: 1548

Number of mitochondrial outliers with `MAD > 5` or `pct_counts_mt > 8`: 7739

Filtered out outliers: 8902 (8.90%)

Filtering out genes expressed in less than 20 cells

Removed 14368 genes (28.01%)

#### After count QC:


[<img src="parameter_iteration_20260421_035352_figures/qc_after_total_counts_MAD=5.png" width="23%">](parameter_iteration_20260421_035352_figures/qc_after_total_counts_MAD=5.png)
[<img src="parameter_iteration_20260421_035352_figures/qc_after_pct_counts_mt_MAD=5.png" width="23%">](parameter_iteration_20260421_035352_figures/qc_after_pct_counts_mt_MAD=5.png)
[<img src="parameter_iteration_20260421_035352_figures/qc_after_pct_counts_mt_violin_MAD=5.png" width="23%">](parameter_iteration_20260421_035352_figures/qc_after_pct_counts_mt_violin_MAD=5.png)
[<img src="parameter_iteration_20260421_035352_figures/qc_after_scatter_genes_total_mt_MAD=5.png" width="23%">](parameter_iteration_20260421_035352_figures/qc_after_scatter_genes_total_mt_MAD=5.png)


### 1.3 Doublet QC


Running scrublet with expected doublet rate 0.05

[<img src="parameter_iteration_20260421_035352_figures/qc_scrublet_histogram.png" width="50%">](parameter_iteration_20260421_035352_figures/qc_scrublet_histogram.png)


There are 2 predicted doublets, applying filter

After doublet step: (91096, 14368) (removed 0.00% cells)

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


Removing 343 genes (2.39%)

Shape after exclusion panel filtering: (91096, 14025)

---
## 2. Normalisation


Normalised in reference to total counts

Transformed with log1p

[<img src="parameter_iteration_20260421_035352_figures/normalisation_total_counts_hist.png" width="48%">](parameter_iteration_20260421_035352_figures/normalisation_total_counts_hist.png)
[<img src="parameter_iteration_20260421_035352_figures/normalisation_log1p_hist.png" width="48%">](parameter_iteration_20260421_035352_figures/normalisation_log1p_hist.png)


---
# Hardcoded parameter grid


| parameter | values |
| --- | --- |
| hvg_n_top_genes | [4000] |
| n_pcs | [20, 30] |
| n_neighbors | [10, 15, 20] |
| min_dist | [0.1, 0.3, 0.5] |
| leiden_resolution | [0.3, 0.5, 1.0] |


Running 18 clustering evaluations (1 HVG settings x 2 PC settings x 3 neighbor settings x 3 Leiden resolutions)

Rendering 18 UMAP panels for visualization (3 min_dist settings for each clustering graph)

---
# UMAP summary grid


---
## HVG 4000, PCs 20


Columns are n_neighbors = 10, 15, 20. Each panel contains all Leiden resolutions: 0.3, 0.5, 1.0.

### min_dist = 0.1


[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn10_md0.1_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn10_md0.1_leiden_resolutions.png)
[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn15_md0.1_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn15_md0.1_leiden_resolutions.png)
[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn20_md0.1_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn20_md0.1_leiden_resolutions.png)


### min_dist = 0.3


[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn10_md0.3_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn10_md0.3_leiden_resolutions.png)
[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn15_md0.3_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn15_md0.3_leiden_resolutions.png)
[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn20_md0.3_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn20_md0.3_leiden_resolutions.png)


### min_dist = 0.5


[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn10_md0.5_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn10_md0.5_leiden_resolutions.png)
[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn15_md0.5_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn15_md0.5_leiden_resolutions.png)
[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn20_md0.5_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs20_nn20_md0.5_leiden_resolutions.png)


---
## HVG 4000, PCs 30


Columns are n_neighbors = 10, 15, 20. Each panel contains all Leiden resolutions: 0.3, 0.5, 1.0.

### min_dist = 0.1


[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn10_md0.1_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn10_md0.1_leiden_resolutions.png)
[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn15_md0.1_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn15_md0.1_leiden_resolutions.png)
[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn20_md0.1_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn20_md0.1_leiden_resolutions.png)


### min_dist = 0.3


[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn10_md0.3_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn10_md0.3_leiden_resolutions.png)
[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn15_md0.3_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn15_md0.3_leiden_resolutions.png)
[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn20_md0.3_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn20_md0.3_leiden_resolutions.png)


### min_dist = 0.5


[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn10_md0.5_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn10_md0.5_leiden_resolutions.png)
[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn15_md0.5_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn15_md0.5_leiden_resolutions.png)
[<img src="parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn20_md0.5_leiden_resolutions.png" width="31%">](parameter_iteration_20260421_035352_figures/hvg4000_pcs30_nn20_md0.5_leiden_resolutions.png)


*Saving into [parameter_iteration_20260421_035352_data/parameter_iteration_summary_statistics.csv](parameter_iteration_20260421_035352_data/parameter_iteration_summary_statistics.csv)*

---
# Summary statistics


These statistics summarize clustering runs only. UMAP `min_dist` is treated as a visualization parameter and is not included in the CSV or summary plots.

Full summary statistics CSV: [parameter_iteration_20260421_035352_data/parameter_iteration_summary_statistics.csv](parameter_iteration_20260421_035352_data/parameter_iteration_summary_statistics.csv)

_Note: Showing at most 50 rows and 10 columns from data with 18 rows and 12 columns._

| index | run_no | name | base_name | n_cells | n_genes | hvg_n_top_genes | n_pcs | n_neighbors | leiden_resolution | n_clusters |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.00 | 1.00 | hvg4000_pcs20_nn10_res0.3 | hvg4000_pcs20_nn10 | 91096.00 | 4000.00 | 4000.00 | 20.00 | 10.00 | 0.30 | 16.00 |
| 1.00 | 2.00 | hvg4000_pcs20_nn10_res0.5 | hvg4000_pcs20_nn10 | 91096.00 | 4000.00 | 4000.00 | 20.00 | 10.00 | 0.50 | 20.00 |
| 2.00 | 3.00 | hvg4000_pcs20_nn10_res1.0 | hvg4000_pcs20_nn10 | 91096.00 | 4000.00 | 4000.00 | 20.00 | 10.00 | 1.00 | 32.00 |
| 3.00 | 4.00 | hvg4000_pcs20_nn15_res0.3 | hvg4000_pcs20_nn15 | 91096.00 | 4000.00 | 4000.00 | 20.00 | 15.00 | 0.30 | 16.00 |
| 4.00 | 5.00 | hvg4000_pcs20_nn15_res0.5 | hvg4000_pcs20_nn15 | 91096.00 | 4000.00 | 4000.00 | 20.00 | 15.00 | 0.50 | 19.00 |
| 5.00 | 6.00 | hvg4000_pcs20_nn15_res1.0 | hvg4000_pcs20_nn15 | 91096.00 | 4000.00 | 4000.00 | 20.00 | 15.00 | 1.00 | 30.00 |
| 6.00 | 7.00 | hvg4000_pcs20_nn20_res0.3 | hvg4000_pcs20_nn20 | 91096.00 | 4000.00 | 4000.00 | 20.00 | 20.00 | 0.30 | 14.00 |
| 7.00 | 8.00 | hvg4000_pcs20_nn20_res0.5 | hvg4000_pcs20_nn20 | 91096.00 | 4000.00 | 4000.00 | 20.00 | 20.00 | 0.50 | 19.00 |
| 8.00 | 9.00 | hvg4000_pcs20_nn20_res1.0 | hvg4000_pcs20_nn20 | 91096.00 | 4000.00 | 4000.00 | 20.00 | 20.00 | 1.00 | 28.00 |
| 9.00 | 10.00 | hvg4000_pcs30_nn10_res0.3 | hvg4000_pcs30_nn10 | 91096.00 | 4000.00 | 4000.00 | 30.00 | 10.00 | 0.30 | 19.00 |
| 10.00 | 11.00 | hvg4000_pcs30_nn10_res0.5 | hvg4000_pcs30_nn10 | 91096.00 | 4000.00 | 4000.00 | 30.00 | 10.00 | 0.50 | 26.00 |
| 11.00 | 12.00 | hvg4000_pcs30_nn10_res1.0 | hvg4000_pcs30_nn10 | 91096.00 | 4000.00 | 4000.00 | 30.00 | 10.00 | 1.00 | 35.00 |
| 12.00 | 13.00 | hvg4000_pcs30_nn15_res0.3 | hvg4000_pcs30_nn15 | 91096.00 | 4000.00 | 4000.00 | 30.00 | 15.00 | 0.30 | 18.00 |
| 13.00 | 14.00 | hvg4000_pcs30_nn15_res0.5 | hvg4000_pcs30_nn15 | 91096.00 | 4000.00 | 4000.00 | 30.00 | 15.00 | 0.50 | 23.00 |
| 14.00 | 15.00 | hvg4000_pcs30_nn15_res1.0 | hvg4000_pcs30_nn15 | 91096.00 | 4000.00 | 4000.00 | 30.00 | 15.00 | 1.00 | 33.00 |
| 15.00 | 16.00 | hvg4000_pcs30_nn20_res0.3 | hvg4000_pcs30_nn20 | 91096.00 | 4000.00 | 4000.00 | 30.00 | 20.00 | 0.30 | 16.00 |
| 16.00 | 17.00 | hvg4000_pcs30_nn20_res0.5 | hvg4000_pcs30_nn20 | 91096.00 | 4000.00 | 4000.00 | 30.00 | 20.00 | 0.50 | 23.00 |
| 17.00 | 18.00 | hvg4000_pcs30_nn20_res1.0 | hvg4000_pcs30_nn20 | 91096.00 | 4000.00 | 4000.00 | 30.00 | 20.00 | 1.00 | 31.00 |


---
## Parameter iteration summary plots


[<img src="parameter_iteration_20260421_035352_figures/parameter_iteration_clusters_vs_dominance.png" width="50%">](parameter_iteration_20260421_035352_figures/parameter_iteration_clusters_vs_dominance.png)


[<img src="parameter_iteration_20260421_035352_figures/parameter_iteration_heatmap_clusters.png" width="50%">](parameter_iteration_20260421_035352_figures/parameter_iteration_heatmap_clusters.png)


[<img src="parameter_iteration_20260421_035352_figures/parameter_iteration_heatmap_smallest_cluster.png" width="50%">](parameter_iteration_20260421_035352_figures/parameter_iteration_heatmap_smallest_cluster.png)


[<img src="parameter_iteration_20260421_035352_figures/parameter_iteration_top_gene_frequency.png" width="50%">](parameter_iteration_20260421_035352_figures/parameter_iteration_top_gene_frequency.png)

