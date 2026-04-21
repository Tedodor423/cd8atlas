# pipeline_elements

- Timestamp: `2026-04-21 02:19:00`
- Source file: `/ceph/project/sharmalab/dnimrich/cd8atlas/code/pipeline_elements.py`


*Loading from [../../data/qc+subsampled_100000.h5ad](../../data/qc+subsampled_100000.h5ad)*

Loaded adata with with shape (91098, 14025)

Preserved existing `counts` layer from loaded adata

Labeled 13256 genes from protein-coding_gene.txt

---
## 3. Feature selection


### 3.3 Highly Variable Gene selection


HVGs selected: 4000

[<img src="pipeline_elements_20260421_021900_figures/hvg_before_seurat_4000.png" width="48%">](pipeline_elements_20260421_021900_figures/hvg_before_seurat_4000.png)
[<img src="pipeline_elements_20260421_021900_figures/hvg_after_seurat_4000.png" width="48%">](pipeline_elements_20260421_021900_figures/hvg_after_seurat_4000.png)


### Top 20 expressed genes after selection:


[<img src="pipeline_elements_20260421_021900_figures/top20_expr_genes.png" width="50%">](pipeline_elements_20260421_021900_figures/top20_expr_genes.png)


---
## 4. Dimensional reduction


### Principal component analysis (PCA)


Scaled data with max variance cutoff 10

Calculated PCA with 30 components

[<img src="pipeline_elements_20260421_021900_figures/pca_variance_ratio_n=30.png" width="48%">](pipeline_elements_20260421_021900_figures/pca_variance_ratio_n=30.png)
[<img src="pipeline_elements_20260421_021900_figures/pca_scatter_n=30.png" width="48%">](pipeline_elements_20260421_021900_figures/pca_scatter_n=30.png)


---
## 5. Clustering


Calculated nearest 15 neighbours using 30 PCs

Calculated UMAP with min_dist 0.3 and spread 1.0

Detected clusters with leiden at resolution 0.5

[<img src="pipeline_elements_20260421_021900_figures/umap_clustering_overview_n=15_mindist=0.3_res=0.5.png" width="50%">](pipeline_elements_20260421_021900_figures/umap_clustering_overview_n=15_mindist=0.3_res=0.5.png)


*Saving into [pipeline_elements_20260421_021900_data/adata_umap_clustering_n=15_mindist=0.3_res=0.5.h5ad](pipeline_elements_20260421_021900_data/adata_umap_clustering_n=15_mindist=0.3_res=0.5.h5ad)*


> Saved adata with shape (91098, 4000)

---
## 6. Dotplot All Markers


### Finding significant genes with FindAllMarkers (via Seurat in R)


FindAllMarkers will use up to 32 R worker(s)

Ran FindAllMarkers in R, keeping up to 3 markers per cluster and selected 67 genes

Plotting 67 genes from adata.var[`findallmarkers_selected`]

[<img src="pipeline_elements_20260421_021900_figures/dotplot_leiden_res_0.5_allmarkers.png" width="50%">](pipeline_elements_20260421_021900_figures/dotplot_leiden_res_0.5_allmarkers.png)

