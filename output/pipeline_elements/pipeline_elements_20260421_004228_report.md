# pipeline_elements

- Timestamp: `2026-04-21 00:42:28`
- Source file: `/ceph/project/sharmalab/dnimrich/cd8atlas/code/pipeline_elements.py`


*Loading from [../../data/adata_umap_clustering_n=15_mindist=0.3_res=0.5.h5ad](../../data/adata_umap_clustering_n=15_mindist=0.3_res=0.5.h5ad)*

Loaded adata with with shape (91098, 4000)

Preserved existing `counts` layer from loaded adata

---
## 6. Dotplot All Markers


### Finding significant genes with FindAllMarkers (via Seurat in R)


FindAllMarkers will use up to 32 R worker(s)

Ran FindAllMarkers in R, selected 1212 genes

Plotting 1212 genes from adata.var[`findallmarkers_selected`]

[<img src="pipeline_elements_20260421_004228_figures/dotplot_leiden_res_0.5_allmarkers.png" width="50%">](pipeline_elements_20260421_004228_figures/dotplot_leiden_res_0.5_allmarkers.png)

