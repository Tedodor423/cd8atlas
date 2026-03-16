import os, time, matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
import anndata as ad
from anndata import AnnData
import scipy.io
import scipy.sparse as sp
import scrublet as scr
import decoupler as dc
import celltypist
from celltypist import models
import PyComplexHeatmap as pch

from pathlib import Path

plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.figsize"] = [5, 4]

root_dir = Path("/ceph/project/sharmalab/dnimrich/cd8atlas")

def run_pipeline(data_dir = root_dir / "data/",
                 input_h5ad = root_dir / "data/" / "huARdb_v2_GEX.CD8.hvg4k.h5ad",
                 output_dir = root_dir / "output/test_pipeline",
                 max_cells: int = 100_000,
                 random_seed: int = 42,
                 min_genes: int = 200,
                 max_mito_pct: float = 20.0,
            ) -> Path:

    # load data
    adata = sc.read_h5ad(str(input_h5ad))
    print(f"Loaded: {input_h5ad} with shape {adata.shape}")

    if max_cells and adata.n_obs > max_cells:
        idx = np.random.default_rng(random_seed).choice(adata.n_obs, size=max_cells, replace=False)
        adata = adata[idx, :].copy()
        print(f"After subsampling: {adata.shape}")

    # QC on raw counts: mitochondrial content + low-complexity cells
    var_upper = adata.var_names.str.upper()
    adata.var["mt"] = var_upper.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    before_qc = adata.n_obs
    keep_qc = (adata.obs["n_genes_by_counts"] >= min_genes) & (adata.obs["pct_counts_mt"] <= max_mito_pct)
    adata = adata[keep_qc, :].copy()
    print(f"After mito/gene QC: {adata.shape} (kept {adata.n_obs}/{before_qc})")

    # Doublet detection/removal (Scrublet)
    before_doublet = adata.n_obs
    try:
        counts = adata.X
        if sp.issparse(counts):
            counts = counts.tocsc()
        else:
            counts = np.asarray(counts)

        scrub = scr.Scrublet(counts)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        adata.obs["scrublet_score"] = doublet_scores
        adata.obs["predicted_doublet"] = predicted_doublets
        adata = adata[~adata.obs["predicted_doublet"].astype(bool), :].copy()
        print(f"After doublet removal: {adata.shape} (kept {adata.n_obs}/{before_doublet})")
    except Exception as e:
        print(f"Warning: Scrublet failed ({e}); continuing without doublet filtering.")

    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
    sc.pp.log1p(adata)
    adata.X = adata.X.astype(np.float32)

    if "highly_variable" in adata.var:
        n_hvg = int(adata.var["highly_variable"].sum())
        if n_hvg > 0:
            adata = adata[:, adata.var["highly_variable"]].copy()
            print(f"Using pre-annotated HVGs: {n_hvg}")
    else:
        print("No HVG mask found; using current features as-is.")

    # Embedding and clustering
    n_pcs = int(min(20, adata.n_vars - 1)) if adata.n_vars > 1 else 1
    sc.tl.pca(adata, n_comps=n_pcs, zero_center=False, svd_solver="randomized")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_pcs, random_state=random_seed)
    sc.tl.umap(adata, min_dist=0.5, spread=1.0, maxiter=200, random_state=random_seed)
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_res_05")

    # save figures
    output_dir.mkdir(parents=True, exist_ok=True)

    sc.pl.umap(adata, color=["leiden_res_05"], show=False)
    plt.savefig(output_dir / "umap_leiden_res_05.png", bbox_inches="tight")
    plt.close()

    qc_cols = [c for c in ["n_genes_by_counts", "total_counts", "pct_counts_mt"] if c in adata.obs]
    if qc_cols:
        sc.pl.violin(adata, qc_cols, jitter=0.2, multi_panel=True, show=False)
        plt.savefig(output_dir / "qc_violin.png", bbox_inches="tight")
        plt.close()

    if "scrublet_score" in adata.obs:
        fig, ax = plt.subplots(figsize=(5, 4))
        sns.histplot(adata.obs["scrublet_score"], bins=40, ax=ax)
        ax.set_title("Scrublet score distribution")
        ax.set_xlabel("scrublet_score")
        plt.tight_layout()
        plt.savefig(output_dir / "scrublet_scores_hist.png", bbox_inches="tight")
        plt.close(fig)

    # save processed data
    out_h5ad = data_dir / "adata_processed.h5ad"
    adata.write_h5ad(out_h5ad)
    
    return out_h5ad


if __name__ == "__main__":
    print(f"Saved: {run_pipeline()}")
