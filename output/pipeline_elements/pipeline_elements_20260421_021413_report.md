# pipeline_elements

- Timestamp: `2026-04-21 02:14:13`
- Source file: `/ceph/project/sharmalab/dnimrich/cd8atlas/code/pipeline_elements.py`


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


[<img src="pipeline_elements_20260421_021413_figures/qc_before_total_counts_hist.png" width="23%">](pipeline_elements_20260421_021413_figures/qc_before_total_counts_hist.png)
[<img src="pipeline_elements_20260421_021413_figures/qc_before_pct_counts_mt_hist.png" width="23%">](pipeline_elements_20260421_021413_figures/qc_before_pct_counts_mt_hist.png)
[<img src="pipeline_elements_20260421_021413_figures/qc_before_pct_counts_mt_violin.png" width="23%">](pipeline_elements_20260421_021413_figures/qc_before_pct_counts_mt_violin.png)
[<img src="pipeline_elements_20260421_021413_figures/qc_before_scatter_genes_total_mt.png" width="23%">](pipeline_elements_20260421_021413_figures/qc_before_scatter_genes_total_mt.png)


Number of count outliers with `MAD > 5`: 1548

Number of mitochondrial outliers with `MAD > 5` or `pct_counts_mt > 8`: 7739

Filtered out outliers: 8902 (8.90%)

Filtering out genes expressed in less than 20 cells

Removed 14368 genes (28.01%)

#### After count QC:


[<img src="pipeline_elements_20260421_021413_figures/qc_after_total_counts_MAD=5.png" width="23%">](pipeline_elements_20260421_021413_figures/qc_after_total_counts_MAD=5.png)
[<img src="pipeline_elements_20260421_021413_figures/qc_after_pct_counts_mt_MAD=5.png" width="23%">](pipeline_elements_20260421_021413_figures/qc_after_pct_counts_mt_MAD=5.png)
[<img src="pipeline_elements_20260421_021413_figures/qc_after_pct_counts_mt_violin_MAD=5.png" width="23%">](pipeline_elements_20260421_021413_figures/qc_after_pct_counts_mt_violin_MAD=5.png)
[<img src="pipeline_elements_20260421_021413_figures/qc_after_scatter_genes_total_mt_MAD=5.png" width="23%">](pipeline_elements_20260421_021413_figures/qc_after_scatter_genes_total_mt_MAD=5.png)


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

Shape after exclusion panel filtering: (91098, 14025)

---
## 2. Normalisation


Normalised in reference to total counts

Transformed with log1p

[<img src="pipeline_elements_20260421_021413_figures/normalisation_total_counts_hist.png" width="48%">](pipeline_elements_20260421_021413_figures/normalisation_total_counts_hist.png)
[<img src="pipeline_elements_20260421_021413_figures/normalisation_log1p_hist.png" width="48%">](pipeline_elements_20260421_021413_figures/normalisation_log1p_hist.png)


*Saving into [../qc+subsampled_100000.h5ad](../qc+subsampled_100000.h5ad)*


> Saved adata with shape (91098, 14025)
