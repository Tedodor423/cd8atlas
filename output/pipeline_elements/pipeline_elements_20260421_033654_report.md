# pipeline_elements

- Timestamp: `2026-04-21 03:36:54`
- Source file: `/ceph/project/sharmalab/dnimrich/cd8atlas/code/pipeline_elements.py`


*Loading from [../../data/qc+subsampled_100000.h5ad](../../data/qc+subsampled_100000.h5ad)*

Loaded adata with with shape (91098, 14025)

Preserved existing `counts` layer from loaded adata

Labeled 13256 genes from protein-coding_gene.txt

---
## 3. Feature selection


### 3.3 Highly Variable Gene selection


HVGs selected: 4000

[<img src="pipeline_elements_20260421_033654_figures/hvg_before_seurat_4000.png" width="48%">](pipeline_elements_20260421_033654_figures/hvg_before_seurat_4000.png)
[<img src="pipeline_elements_20260421_033654_figures/hvg_after_seurat_4000.png" width="48%">](pipeline_elements_20260421_033654_figures/hvg_after_seurat_4000.png)


### Top 20 expressed genes after selection:


[<img src="pipeline_elements_20260421_033654_figures/top20_expr_genes.png" width="50%">](pipeline_elements_20260421_033654_figures/top20_expr_genes.png)


---
## 4. Dimensional reduction


### Principal component analysis (PCA)


Scaled data with max variance cutoff 10

Calculated PCA with 30 components

[<img src="pipeline_elements_20260421_033654_figures/pca_variance_ratio_n=30.png" width="48%">](pipeline_elements_20260421_033654_figures/pca_variance_ratio_n=30.png)
[<img src="pipeline_elements_20260421_033654_figures/pca_scatter_n=30.png" width="48%">](pipeline_elements_20260421_033654_figures/pca_scatter_n=30.png)


---
## 5. UMAP


Calculated nearest 15 neighbours using 30 PCs

Calculated UMAP with min_dist 0.3 and spread 1.0

---
## 6. Clustering


Detected clusters with leiden at resolution 0.5

[<img src="pipeline_elements_20260421_033654_figures/umap_plot_n=15_mindist=0.3_res=0.5.png" width="50%">](pipeline_elements_20260421_033654_figures/umap_plot_n=15_mindist=0.3_res=0.5.png)


*Saving into [pipeline_elements_20260421_033654_data/adata_umap_clustering_n=15_mindist=0.3_res=0.5.h5ad](pipeline_elements_20260421_033654_data/adata_umap_clustering_n=15_mindist=0.3_res=0.5.h5ad)*


> Saved adata with shape (91098, 4000)

### State labeling from state_markers.csv


| state | n_markers_listed | n_markers_matched | n_markers_missing | matched_markers | missing_markers |
| --- | --- | --- | --- | --- | --- |
| cytotoxic | 5.00 | 4.00 | 1.00 | GNLY, GZMA, GZMB, PRF1 | IFNG |
| effector | 9.00 | 8.00 | 1.00 | BHLHE40, CX3CR1, KLRG1, PRDM1, RO... | 1D2 |
| exhausted | 9.00 | 8.00 | 1.00 | BATF, CTLA4, ENTPD1, LAG3, PDCD1,... | TIGID |
| naive/memory | 8.00 | 7.00 | 1.00 | CCR7, CD28, CXCR3, LEF1, SELL, TC... | 1D3 |


State marker coverage in current dataset

| n_cells | n_genes_in_dataset | n_unique_state_markers_listed | n_unique_state_markers_matched | n_unique_state_markers_missing | pct_state_markers_matched |
| --- | --- | --- | --- | --- | --- |
| 91098.00 | 4000.00 | 31.00 | 27.00 | 4.00 | 87.10 |


Missing state markers: 1D2, 1D3, IFNG, TIGID

Assigned per-cell state labels

| index | state | n_cells |
| --- | --- | --- |
| 0.00 | cytotoxic | 33294.00 |
| 1.00 | naive/memory | 29498.00 |
| 2.00 | effector | 16136.00 |
| 3.00 | exhausted | 12170.00 |


[<img src="pipeline_elements_20260421_033654_figures/umap_plot_states.png" width="50%">](pipeline_elements_20260421_033654_figures/umap_plot_states.png)


---
## 6. Dotplot All Markers


### Finding significant genes with FindAllMarkers (via Seurat in R)


FindAllMarkers will use up to 32 R worker(s)

Ran FindAllMarkers in R, keeping up to 3 markers per cluster and selected 67 genes

Plotting 67 genes from adata.var[`findallmarkers_selected`]

[<img src="pipeline_elements_20260421_033654_figures/dotplot_leiden_res_0.5_allmarkers.png" width="50%">](pipeline_elements_20260421_033654_figures/dotplot_leiden_res_0.5_allmarkers.png)

