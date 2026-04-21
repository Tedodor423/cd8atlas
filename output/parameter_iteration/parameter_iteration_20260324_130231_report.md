# parameter_iteration

- Timestamp: `2026-03-24 13:02:31`
- Source file: `/ceph/project/sharmalab/dnimrich/cd8atlas/code/pipeline_elements.py`


---
# Loading and preprocessing data


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


[<img src="parameter_iteration_20260324_130231_figures/qc_before_total_counts_hist.png" width="23%">](parameter_iteration_20260324_130231_figures/qc_before_total_counts_hist.png)
[<img src="parameter_iteration_20260324_130231_figures/qc_before_pct_counts_mt_hist.png" width="23%">](parameter_iteration_20260324_130231_figures/qc_before_pct_counts_mt_hist.png)
[<img src="parameter_iteration_20260324_130231_figures/qc_before_pct_counts_mt_violin.png" width="23%">](parameter_iteration_20260324_130231_figures/qc_before_pct_counts_mt_violin.png)
[<img src="parameter_iteration_20260324_130231_figures/qc_before_scatter_genes_total_mt.png" width="23%">](parameter_iteration_20260324_130231_figures/qc_before_scatter_genes_total_mt.png)


Number of count outliers with `MAD > 5`: 1548

Number of mitochondrial outliers with `MAD > 5` or `pct_counts_mt > 8`: 7739

Filtered out outliers: 8902 (8.90%)

Filtering out genes expressed in less than 20 cells

Removed 14368 genes (28.01%)

#### After count QC:


[<img src="parameter_iteration_20260324_130231_figures/qc_after_total_counts_MAD=5.png" width="23%">](parameter_iteration_20260324_130231_figures/qc_after_total_counts_MAD=5.png)
[<img src="parameter_iteration_20260324_130231_figures/qc_after_pct_counts_mt_MAD=5.png" width="23%">](parameter_iteration_20260324_130231_figures/qc_after_pct_counts_mt_MAD=5.png)
[<img src="parameter_iteration_20260324_130231_figures/qc_after_pct_counts_mt_violin_MAD=5.png" width="23%">](parameter_iteration_20260324_130231_figures/qc_after_pct_counts_mt_violin_MAD=5.png)
[<img src="parameter_iteration_20260324_130231_figures/qc_after_scatter_genes_total_mt_MAD=5.png" width="23%">](parameter_iteration_20260324_130231_figures/qc_after_scatter_genes_total_mt_MAD=5.png)


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

[<img src="parameter_iteration_20260324_130231_figures/normalisation_total_counts_hist.png" width="48%">](parameter_iteration_20260324_130231_figures/normalisation_total_counts_hist.png)
[<img src="parameter_iteration_20260324_130231_figures/normalisation_log1p_hist.png" width="48%">](parameter_iteration_20260324_130231_figures/normalisation_log1p_hist.png)


*Saving into [../parameter-iteration_preprocessed.h5ad](../parameter-iteration_preprocessed.h5ad)*


> Saved adata with shape (91098, 14028)

---
# Loaded parameters from /ceph/project/sharmalab/dnimrich/cd8atlas/code/parameter_iteration/parameters.csv:


| index | hvg_n_top_genes | n_pcs | n_neighbors | leiden_resolution | min_dist |
| --- | --- | --- | --- | --- | --- |
| 0.00 | 2000.00 | 20.00 | 10.00 | 0.30 | 0.10 |
| 1.00 | 3000.00 | 30.00 | 15.00 | 0.50 | 0.30 |
| 2.00 | 4000.00 | nan | 20.00 | 1.00 | 0.50 |


Running with 162 parameter combinations

---
## Run 1/162: hvg2000_pcs(20,)_nn(10,)_res0.3_md0.1


---
## Run 2/162: hvg2000_pcs(20,)_nn(10,)_res0.3_md0.3


---
## Run 3/162: hvg2000_pcs(20,)_nn(10,)_res0.3_md0.5


---
## Run 4/162: hvg2000_pcs(20,)_nn(10,)_res0.5_md0.1


---
## Run 5/162: hvg2000_pcs(20,)_nn(10,)_res0.5_md0.3


---
## Run 6/162: hvg2000_pcs(20,)_nn(10,)_res0.5_md0.5


---
## Run 7/162: hvg2000_pcs(20,)_nn(10,)_res1.0_md0.1


---
## Run 8/162: hvg2000_pcs(20,)_nn(10,)_res1.0_md0.3


---
## Run 9/162: hvg2000_pcs(20,)_nn(10,)_res1.0_md0.5


---
## Run 10/162: hvg2000_pcs(20,)_nn(15,)_res0.3_md0.1


---
## Run 11/162: hvg2000_pcs(20,)_nn(15,)_res0.3_md0.3


---
## Run 12/162: hvg2000_pcs(20,)_nn(15,)_res0.3_md0.5


---
## Run 13/162: hvg2000_pcs(20,)_nn(15,)_res0.5_md0.1


---
## Run 14/162: hvg2000_pcs(20,)_nn(15,)_res0.5_md0.3


---
## Run 15/162: hvg2000_pcs(20,)_nn(15,)_res0.5_md0.5


---
## Run 16/162: hvg2000_pcs(20,)_nn(15,)_res1.0_md0.1


---
## Run 17/162: hvg2000_pcs(20,)_nn(15,)_res1.0_md0.3


---
## Run 18/162: hvg2000_pcs(20,)_nn(15,)_res1.0_md0.5


---
## Run 19/162: hvg2000_pcs(20,)_nn(20,)_res0.3_md0.1


---
## Run 20/162: hvg2000_pcs(20,)_nn(20,)_res0.3_md0.3


---
## Run 21/162: hvg2000_pcs(20,)_nn(20,)_res0.3_md0.5


---
## Run 22/162: hvg2000_pcs(20,)_nn(20,)_res0.5_md0.1


---
## Run 23/162: hvg2000_pcs(20,)_nn(20,)_res0.5_md0.3


---
## Run 24/162: hvg2000_pcs(20,)_nn(20,)_res0.5_md0.5


---
## Run 25/162: hvg2000_pcs(20,)_nn(20,)_res1.0_md0.1


---
## Run 26/162: hvg2000_pcs(20,)_nn(20,)_res1.0_md0.3


---
## Run 27/162: hvg2000_pcs(20,)_nn(20,)_res1.0_md0.5


---
## Run 28/162: hvg2000_pcs(30,)_nn(10,)_res0.3_md0.1


---
## Run 29/162: hvg2000_pcs(30,)_nn(10,)_res0.3_md0.3


---
## Run 30/162: hvg2000_pcs(30,)_nn(10,)_res0.3_md0.5


---
## Run 31/162: hvg2000_pcs(30,)_nn(10,)_res0.5_md0.1


---
## Run 32/162: hvg2000_pcs(30,)_nn(10,)_res0.5_md0.3


---
## Run 33/162: hvg2000_pcs(30,)_nn(10,)_res0.5_md0.5


---
## Run 34/162: hvg2000_pcs(30,)_nn(10,)_res1.0_md0.1


---
## Run 35/162: hvg2000_pcs(30,)_nn(10,)_res1.0_md0.3


---
## Run 36/162: hvg2000_pcs(30,)_nn(10,)_res1.0_md0.5


---
## Run 37/162: hvg2000_pcs(30,)_nn(15,)_res0.3_md0.1


---
## Run 38/162: hvg2000_pcs(30,)_nn(15,)_res0.3_md0.3


---
## Run 39/162: hvg2000_pcs(30,)_nn(15,)_res0.3_md0.5


---
## Run 40/162: hvg2000_pcs(30,)_nn(15,)_res0.5_md0.1


---
## Run 41/162: hvg2000_pcs(30,)_nn(15,)_res0.5_md0.3


---
## Run 42/162: hvg2000_pcs(30,)_nn(15,)_res0.5_md0.5


---
## Run 43/162: hvg2000_pcs(30,)_nn(15,)_res1.0_md0.1


---
## Run 44/162: hvg2000_pcs(30,)_nn(15,)_res1.0_md0.3


---
## Run 45/162: hvg2000_pcs(30,)_nn(15,)_res1.0_md0.5


---
## Run 46/162: hvg2000_pcs(30,)_nn(20,)_res0.3_md0.1


---
## Run 47/162: hvg2000_pcs(30,)_nn(20,)_res0.3_md0.3


---
## Run 48/162: hvg2000_pcs(30,)_nn(20,)_res0.3_md0.5


---
## Run 49/162: hvg2000_pcs(30,)_nn(20,)_res0.5_md0.1


---
## Run 50/162: hvg2000_pcs(30,)_nn(20,)_res0.5_md0.3


---
## Run 51/162: hvg2000_pcs(30,)_nn(20,)_res0.5_md0.5


---
## Run 52/162: hvg2000_pcs(30,)_nn(20,)_res1.0_md0.1


---
## Run 53/162: hvg2000_pcs(30,)_nn(20,)_res1.0_md0.3


---
## Run 54/162: hvg2000_pcs(30,)_nn(20,)_res1.0_md0.5


---
## Run 55/162: hvg3000_pcs(20,)_nn(10,)_res0.3_md0.1


---
## Run 56/162: hvg3000_pcs(20,)_nn(10,)_res0.3_md0.3


---
## Run 57/162: hvg3000_pcs(20,)_nn(10,)_res0.3_md0.5


---
## Run 58/162: hvg3000_pcs(20,)_nn(10,)_res0.5_md0.1


---
## Run 59/162: hvg3000_pcs(20,)_nn(10,)_res0.5_md0.3


---
## Run 60/162: hvg3000_pcs(20,)_nn(10,)_res0.5_md0.5


---
## Run 61/162: hvg3000_pcs(20,)_nn(10,)_res1.0_md0.1


---
## Run 62/162: hvg3000_pcs(20,)_nn(10,)_res1.0_md0.3


---
## Run 63/162: hvg3000_pcs(20,)_nn(10,)_res1.0_md0.5


---
## Run 64/162: hvg3000_pcs(20,)_nn(15,)_res0.3_md0.1


---
## Run 65/162: hvg3000_pcs(20,)_nn(15,)_res0.3_md0.3


---
## Run 66/162: hvg3000_pcs(20,)_nn(15,)_res0.3_md0.5


---
## Run 67/162: hvg3000_pcs(20,)_nn(15,)_res0.5_md0.1


---
## Run 68/162: hvg3000_pcs(20,)_nn(15,)_res0.5_md0.3


---
## Run 69/162: hvg3000_pcs(20,)_nn(15,)_res0.5_md0.5


---
## Run 70/162: hvg3000_pcs(20,)_nn(15,)_res1.0_md0.1


---
## Run 71/162: hvg3000_pcs(20,)_nn(15,)_res1.0_md0.3


---
## Run 72/162: hvg3000_pcs(20,)_nn(15,)_res1.0_md0.5


---
## Run 73/162: hvg3000_pcs(20,)_nn(20,)_res0.3_md0.1


---
## Run 74/162: hvg3000_pcs(20,)_nn(20,)_res0.3_md0.3


---
## Run 75/162: hvg3000_pcs(20,)_nn(20,)_res0.3_md0.5


---
## Run 76/162: hvg3000_pcs(20,)_nn(20,)_res0.5_md0.1


---
## Run 77/162: hvg3000_pcs(20,)_nn(20,)_res0.5_md0.3


---
## Run 78/162: hvg3000_pcs(20,)_nn(20,)_res0.5_md0.5


---
## Run 79/162: hvg3000_pcs(20,)_nn(20,)_res1.0_md0.1


---
## Run 80/162: hvg3000_pcs(20,)_nn(20,)_res1.0_md0.3


---
## Run 81/162: hvg3000_pcs(20,)_nn(20,)_res1.0_md0.5


---
## Run 82/162: hvg3000_pcs(30,)_nn(10,)_res0.3_md0.1


---
## Run 83/162: hvg3000_pcs(30,)_nn(10,)_res0.3_md0.3


---
## Run 84/162: hvg3000_pcs(30,)_nn(10,)_res0.3_md0.5


---
## Run 85/162: hvg3000_pcs(30,)_nn(10,)_res0.5_md0.1


---
## Run 86/162: hvg3000_pcs(30,)_nn(10,)_res0.5_md0.3


---
## Run 87/162: hvg3000_pcs(30,)_nn(10,)_res0.5_md0.5


---
## Run 88/162: hvg3000_pcs(30,)_nn(10,)_res1.0_md0.1


---
## Run 89/162: hvg3000_pcs(30,)_nn(10,)_res1.0_md0.3


---
## Run 90/162: hvg3000_pcs(30,)_nn(10,)_res1.0_md0.5


---
## Run 91/162: hvg3000_pcs(30,)_nn(15,)_res0.3_md0.1


---
## Run 92/162: hvg3000_pcs(30,)_nn(15,)_res0.3_md0.3


---
## Run 93/162: hvg3000_pcs(30,)_nn(15,)_res0.3_md0.5


---
## Run 94/162: hvg3000_pcs(30,)_nn(15,)_res0.5_md0.1


---
## Run 95/162: hvg3000_pcs(30,)_nn(15,)_res0.5_md0.3


---
## Run 96/162: hvg3000_pcs(30,)_nn(15,)_res0.5_md0.5


---
## Run 97/162: hvg3000_pcs(30,)_nn(15,)_res1.0_md0.1


---
## Run 98/162: hvg3000_pcs(30,)_nn(15,)_res1.0_md0.3


---
## Run 99/162: hvg3000_pcs(30,)_nn(15,)_res1.0_md0.5


---
## Run 100/162: hvg3000_pcs(30,)_nn(20,)_res0.3_md0.1


---
## Run 101/162: hvg3000_pcs(30,)_nn(20,)_res0.3_md0.3


---
## Run 102/162: hvg3000_pcs(30,)_nn(20,)_res0.3_md0.5


---
## Run 103/162: hvg3000_pcs(30,)_nn(20,)_res0.5_md0.1


---
## Run 104/162: hvg3000_pcs(30,)_nn(20,)_res0.5_md0.3


---
## Run 105/162: hvg3000_pcs(30,)_nn(20,)_res0.5_md0.5


---
## Run 106/162: hvg3000_pcs(30,)_nn(20,)_res1.0_md0.1


---
## Run 107/162: hvg3000_pcs(30,)_nn(20,)_res1.0_md0.3


---
## Run 108/162: hvg3000_pcs(30,)_nn(20,)_res1.0_md0.5


---
## Run 109/162: hvg4000_pcs(20,)_nn(10,)_res0.3_md0.1


---
## Run 110/162: hvg4000_pcs(20,)_nn(10,)_res0.3_md0.3


---
## Run 111/162: hvg4000_pcs(20,)_nn(10,)_res0.3_md0.5


---
## Run 112/162: hvg4000_pcs(20,)_nn(10,)_res0.5_md0.1


---
## Run 113/162: hvg4000_pcs(20,)_nn(10,)_res0.5_md0.3


---
## Run 114/162: hvg4000_pcs(20,)_nn(10,)_res0.5_md0.5


---
## Run 115/162: hvg4000_pcs(20,)_nn(10,)_res1.0_md0.1


---
## Run 116/162: hvg4000_pcs(20,)_nn(10,)_res1.0_md0.3


---
## Run 117/162: hvg4000_pcs(20,)_nn(10,)_res1.0_md0.5


---
## Run 118/162: hvg4000_pcs(20,)_nn(15,)_res0.3_md0.1


---
## Run 119/162: hvg4000_pcs(20,)_nn(15,)_res0.3_md0.3


---
## Run 120/162: hvg4000_pcs(20,)_nn(15,)_res0.3_md0.5


---
## Run 121/162: hvg4000_pcs(20,)_nn(15,)_res0.5_md0.1


---
## Run 122/162: hvg4000_pcs(20,)_nn(15,)_res0.5_md0.3


---
## Run 123/162: hvg4000_pcs(20,)_nn(15,)_res0.5_md0.5


---
## Run 124/162: hvg4000_pcs(20,)_nn(15,)_res1.0_md0.1


---
## Run 125/162: hvg4000_pcs(20,)_nn(15,)_res1.0_md0.3


---
## Run 126/162: hvg4000_pcs(20,)_nn(15,)_res1.0_md0.5


---
## Run 127/162: hvg4000_pcs(20,)_nn(20,)_res0.3_md0.1


---
## Run 128/162: hvg4000_pcs(20,)_nn(20,)_res0.3_md0.3


---
## Run 129/162: hvg4000_pcs(20,)_nn(20,)_res0.3_md0.5


---
## Run 130/162: hvg4000_pcs(20,)_nn(20,)_res0.5_md0.1


---
## Run 131/162: hvg4000_pcs(20,)_nn(20,)_res0.5_md0.3


---
## Run 132/162: hvg4000_pcs(20,)_nn(20,)_res0.5_md0.5


---
## Run 133/162: hvg4000_pcs(20,)_nn(20,)_res1.0_md0.1


---
## Run 134/162: hvg4000_pcs(20,)_nn(20,)_res1.0_md0.3


---
## Run 135/162: hvg4000_pcs(20,)_nn(20,)_res1.0_md0.5


---
## Run 136/162: hvg4000_pcs(30,)_nn(10,)_res0.3_md0.1


---
## Run 137/162: hvg4000_pcs(30,)_nn(10,)_res0.3_md0.3


---
## Run 138/162: hvg4000_pcs(30,)_nn(10,)_res0.3_md0.5


---
## Run 139/162: hvg4000_pcs(30,)_nn(10,)_res0.5_md0.1


---
## Run 140/162: hvg4000_pcs(30,)_nn(10,)_res0.5_md0.3


---
## Run 141/162: hvg4000_pcs(30,)_nn(10,)_res0.5_md0.5


---
## Run 142/162: hvg4000_pcs(30,)_nn(10,)_res1.0_md0.1


---
## Run 143/162: hvg4000_pcs(30,)_nn(10,)_res1.0_md0.3


---
## Run 144/162: hvg4000_pcs(30,)_nn(10,)_res1.0_md0.5


---
## Run 145/162: hvg4000_pcs(30,)_nn(15,)_res0.3_md0.1


---
## Run 146/162: hvg4000_pcs(30,)_nn(15,)_res0.3_md0.3


---
## Run 147/162: hvg4000_pcs(30,)_nn(15,)_res0.3_md0.5


---
## Run 148/162: hvg4000_pcs(30,)_nn(15,)_res0.5_md0.1


---
## Run 149/162: hvg4000_pcs(30,)_nn(15,)_res0.5_md0.3


---
## Run 150/162: hvg4000_pcs(30,)_nn(15,)_res0.5_md0.5


---
## Run 151/162: hvg4000_pcs(30,)_nn(15,)_res1.0_md0.1


---
## Run 152/162: hvg4000_pcs(30,)_nn(15,)_res1.0_md0.3


---
## Run 153/162: hvg4000_pcs(30,)_nn(15,)_res1.0_md0.5


---
## Run 154/162: hvg4000_pcs(30,)_nn(20,)_res0.3_md0.1


---
## Run 155/162: hvg4000_pcs(30,)_nn(20,)_res0.3_md0.3


---
## Run 156/162: hvg4000_pcs(30,)_nn(20,)_res0.3_md0.5


---
## Run 157/162: hvg4000_pcs(30,)_nn(20,)_res0.5_md0.1


---
## Run 158/162: hvg4000_pcs(30,)_nn(20,)_res0.5_md0.3


---
## Run 159/162: hvg4000_pcs(30,)_nn(20,)_res0.5_md0.5


---
## Run 160/162: hvg4000_pcs(30,)_nn(20,)_res1.0_md0.1


---
## Run 161/162: hvg4000_pcs(30,)_nn(20,)_res1.0_md0.3


---
## Run 162/162: hvg4000_pcs(30,)_nn(20,)_res1.0_md0.5


_Note: Showing at most 10 rows and 10 columns from data with 162 rows and 13 columns._

| run_no | name | n_cells | n_genes | hvg_n_top_genes | n_pcs | n_neighbors | leiden_resolution | min_dist | n_clusters |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.00 | hvg2000_pcs(20,)_nn(10,)_res0.3_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 0.30 | 0.10 | 17.00 |
| 2.00 | hvg2000_pcs(20,)_nn(10,)_res0.3_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 0.30 | 0.30 | 17.00 |
| 3.00 | hvg2000_pcs(20,)_nn(10,)_res0.3_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 0.30 | 0.50 | 17.00 |
| 4.00 | hvg2000_pcs(20,)_nn(10,)_res0.5_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 0.50 | 0.10 | 22.00 |
| 5.00 | hvg2000_pcs(20,)_nn(10,)_res0.5_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 0.50 | 0.30 | 22.00 |
| 6.00 | hvg2000_pcs(20,)_nn(10,)_res0.5_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 0.50 | 0.50 | 22.00 |
| 7.00 | hvg2000_pcs(20,)_nn(10,)_res1.0_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 1.00 | 0.10 | 35.00 |
| 8.00 | hvg2000_pcs(20,)_nn(10,)_res1.0_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 1.00 | 0.30 | 35.00 |
| 9.00 | hvg2000_pcs(20,)_nn(10,)_res1.0_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 1.00 | 0.50 | 35.00 |
| 10.00 | hvg2000_pcs(20,)_nn(15,)_res0.3_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 15.00 | 0.30 | 0.10 | 14.00 |


---
# Parameter iteration summary plots


_Note: Showing at most 10 rows and 10 columns from data with 162 rows and 12 columns._

| index | run_no | name | n_cells | n_genes | hvg_n_top_genes | n_pcs | n_neighbors | leiden_resolution | min_dist | n_clusters |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.00 | 1.00 | hvg2000_pcs(20,)_nn(10,)_res0.3_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 0.30 | 0.10 | 17.00 |
| 1.00 | 2.00 | hvg2000_pcs(20,)_nn(10,)_res0.3_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 0.30 | 0.30 | 17.00 |
| 2.00 | 3.00 | hvg2000_pcs(20,)_nn(10,)_res0.3_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 0.30 | 0.50 | 17.00 |
| 3.00 | 4.00 | hvg2000_pcs(20,)_nn(10,)_res0.5_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 0.50 | 0.10 | 22.00 |
| 4.00 | 5.00 | hvg2000_pcs(20,)_nn(10,)_res0.5_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 0.50 | 0.30 | 22.00 |
| 5.00 | 6.00 | hvg2000_pcs(20,)_nn(10,)_res0.5_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 0.50 | 0.50 | 22.00 |
| 6.00 | 7.00 | hvg2000_pcs(20,)_nn(10,)_res1.0_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 1.00 | 0.10 | 35.00 |
| 7.00 | 8.00 | hvg2000_pcs(20,)_nn(10,)_res1.0_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 1.00 | 0.30 | 35.00 |
| 8.00 | 9.00 | hvg2000_pcs(20,)_nn(10,)_res1.0_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 10.00 | 1.00 | 0.50 | 35.00 |
| 9.00 | 10.00 | hvg2000_pcs(20,)_nn(15,)_res0.3_m... | 91098.00 | 4000.00 | 2000.00 | 20.00 | 15.00 | 0.30 | 0.10 | 14.00 |


[<img src="parameter_iteration_20260324_130231_figures/parameter_iteration_clusters_vs_dominance.png" width="50%">](parameter_iteration_20260324_130231_figures/parameter_iteration_clusters_vs_dominance.png)


[<img src="parameter_iteration_20260324_130231_figures/parameter_iteration_heatmap_clusters.png" width="50%">](parameter_iteration_20260324_130231_figures/parameter_iteration_heatmap_clusters.png)


[<img src="parameter_iteration_20260324_130231_figures/parameter_iteration_heatmap_smallest_cluster.png" width="50%">](parameter_iteration_20260324_130231_figures/parameter_iteration_heatmap_smallest_cluster.png)


[<img src="parameter_iteration_20260324_130231_figures/parameter_iteration_top_gene_frequency.png" width="50%">](parameter_iteration_20260324_130231_figures/parameter_iteration_top_gene_frequency.png)

