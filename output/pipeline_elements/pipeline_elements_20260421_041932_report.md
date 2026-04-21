# pipeline_elements

- Timestamp: `2026-04-21 04:19:32`
- Source file: `/ceph/project/sharmalab/dnimrich/cd8atlas/code/pipeline_elements.py`


*Loading from [../../data/qc+subsampled_100000.h5ad](../../data/qc+subsampled_100000.h5ad)*

Loaded adata with with shape (91098, 14025)

Preserved existing `counts` layer from loaded adata

Labeled 13256 genes from protein-coding_gene.txt

---
## 3. Feature selection


### 3.3 Highly Variable Gene selection


HVGs selected: 4000

[<img src="pipeline_elements_20260421_041932_figures/hvg_before_seurat_4000.png" width="48%">](pipeline_elements_20260421_041932_figures/hvg_before_seurat_4000.png)
[<img src="pipeline_elements_20260421_041932_figures/hvg_after_seurat_4000.png" width="48%">](pipeline_elements_20260421_041932_figures/hvg_after_seurat_4000.png)


### Top 20 expressed genes after selection:


[<img src="pipeline_elements_20260421_041932_figures/top20_expr_genes.png" width="50%">](pipeline_elements_20260421_041932_figures/top20_expr_genes.png)


---
## 4. Dimensional reduction


### Principal component analysis (PCA)


Scaled data with max variance cutoff 10

Calculated PCA with 30 components

[<img src="pipeline_elements_20260421_041932_figures/pca_variance_ratio_n=30.png" width="48%">](pipeline_elements_20260421_041932_figures/pca_variance_ratio_n=30.png)
[<img src="pipeline_elements_20260421_041932_figures/pca_scatter_n=30.png" width="48%">](pipeline_elements_20260421_041932_figures/pca_scatter_n=30.png)


---
## 5. UMAP


Calculated nearest 15 neighbours using 30 PCs

Calculated UMAP with min_dist 0.3 and spread 1.0

---
## 6. Clustering


Detected clusters with leiden at resolution 0.5

[<img src="pipeline_elements_20260421_041932_figures/umap_plot_n=15_mindist=0.3_res=0.5.png" width="50%">](pipeline_elements_20260421_041932_figures/umap_plot_n=15_mindist=0.3_res=0.5.png)


### Labeling states based on markers:


_Note: Showing at most 10 rows and 10 columns from data with 31 rows and 2 columns._

| index | gene | marker |
| --- | --- | --- |
| 0.00 | IFNG | cytotoxic |
| 1.00 | PRF1 | cytotoxic |
| 2.00 | GZMA | cytotoxic |
| 3.00 | GZMB | cytotoxic |
| 4.00 | GNLY | cytotoxic |
| 5.00 | RORC | effector |
| 6.00 | TBX21 | effector |
| 7.00 | CX3CR1 | effector |
| 8.00 | RORA | effector |
| 9.00 | ZEB2 | effector |


State marker genes in dataset: 27 present, 4 missing

Assigned state labels:

| index | state | n_cells |
| --- | --- | --- |
| 0.00 | cytotoxic | 33294.00 |
| 1.00 | naive/memory | 29498.00 |
| 2.00 | effector | 16136.00 |
| 3.00 | exhausted | 12170.00 |


[<img src="pipeline_elements_20260421_041932_figures/umap_plot_states.png" width="50%">](pipeline_elements_20260421_041932_figures/umap_plot_states.png)


---
## 6. Dotplot All Markers


### Finding significant genes with FindAllMarkers (via Seurat in R)


FindAllMarkers will use up to 32 R worker(s)

Ran FindAllMarkers in R, keeping up to 3 markers per cluster and selected 67 genes

Plotting 67 genes from adata.var[`findallmarkers_selected`]

[<img src="pipeline_elements_20260421_041932_figures/dotplot_leiden_res_0.5_allmarkers.png" width="50%">](pipeline_elements_20260421_041932_figures/dotplot_leiden_res_0.5_allmarkers.png)

