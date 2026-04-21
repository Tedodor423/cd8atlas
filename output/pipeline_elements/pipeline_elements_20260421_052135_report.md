# pipeline_elements

- Timestamp: `2026-04-21 05:21:35`
- Source file: `/ceph/project/sharmalab/dnimrich/cd8atlas/code/pipeline_elements.py`


*Loading from [../../data/qc+subsampled_100000.h5ad](../../data/qc+subsampled_100000.h5ad)*

Loaded adata with with shape (91098, 14025)

Preserved existing `counts` layer from loaded adata

Labeled 13256 genes from protein-coding_gene.txt

---
## 3. Feature selection


### 3.3 Highly Variable Gene selection


HVGs selected: 4000 (including 27 whitelisted genes)

[<img src="pipeline_elements_20260421_052135_figures/hvg_before_seurat_4000.png" width="48%">](pipeline_elements_20260421_052135_figures/hvg_before_seurat_4000.png)
[<img src="pipeline_elements_20260421_052135_figures/hvg_after_seurat_4000.png" width="48%">](pipeline_elements_20260421_052135_figures/hvg_after_seurat_4000.png)


### Top 20 expressed genes after selection:


[<img src="pipeline_elements_20260421_052135_figures/top20_expr_genes.png" width="50%">](pipeline_elements_20260421_052135_figures/top20_expr_genes.png)


---
## 4. Dimensional reduction


### Principal component analysis (PCA)


Scaled data with max variance cutoff 10
