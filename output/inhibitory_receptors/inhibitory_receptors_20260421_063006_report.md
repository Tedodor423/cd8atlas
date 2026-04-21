# inhibitory_receptors

- Timestamp: `2026-04-21 06:30:06`
- Source file: `/ceph/project/sharmalab/dnimrich/cd8atlas/code/pipeline_elements.py`


*Loading from [../../data/qc+subsampled_100000.h5ad](../../data/qc+subsampled_100000.h5ad)*

Loaded adata with with shape (91098, 14025)

Preserved existing `counts` layer from loaded adata

Labeled 13256 genes from protein-coding_gene.txt

---
## 3. Feature selection


### 3.3 Highly Variable Gene selection


HVGs selected: 4000 (including 27 whitelisted genes)

[<img src="inhibitory_receptors_20260421_063006_figures/hvg_before_seurat_4000.png" width="48%">](inhibitory_receptors_20260421_063006_figures/hvg_before_seurat_4000.png)
[<img src="inhibitory_receptors_20260421_063006_figures/hvg_after_seurat_4000.png" width="48%">](inhibitory_receptors_20260421_063006_figures/hvg_after_seurat_4000.png)


### Top 20 expressed genes after selection:


[<img src="inhibitory_receptors_20260421_063006_figures/top20_expr_genes.png" width="50%">](inhibitory_receptors_20260421_063006_figures/top20_expr_genes.png)


---
## 4. Dimensional reduction


### Principal component analysis (PCA)


Scaled data with max variance cutoff 10

Calculated PCA with 25 components

[<img src="inhibitory_receptors_20260421_063006_figures/pca_variance_ratio_n=25.png" width="48%">](inhibitory_receptors_20260421_063006_figures/pca_variance_ratio_n=25.png)
[<img src="inhibitory_receptors_20260421_063006_figures/pca_scatter_n=25.png" width="48%">](inhibitory_receptors_20260421_063006_figures/pca_scatter_n=25.png)


### UMAP


Calculated nearest 12 neighbours using 25 PCs

Calculated UMAP with min_dist 0.1 and spread 1.0

---
## 5. Clustering

