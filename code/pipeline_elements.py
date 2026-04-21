import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scrublet as scr
# import celltypist
# from celltypist import models
import scipy.stats as stats
from pathlib import Path
import bioreport
import subprocess
import re

random_seed = 423
root_dir = Path("/ceph/project/sharmalab/dnimrich/cd8atlas")
LEIDEN_COL = "leiden"

def init(output_name: str):
    OUTPUT_NAME = output_name
    bioreport.start(OUTPUT_NAME, path=root_dir / "output", automatic_titles=True, git_image_links='Tedodor423/cd8atlas/refs/heads/main')

def load_data(h5ad_file: str = "huARdb_v2_GEX.CD8.all_genes.h5ad", root=True) -> sc.AnnData:
    h5ad_file = Path(h5ad_file)
    if root:
        adata = sc.read_h5ad(bioreport.file_load(h5ad_file, data_dir = root_dir / "data/"))
    else:
        adata = sc.read_h5ad(bioreport.file_load(h5ad_file))
    bioreport.log(f"Loaded adata with with shape {adata.shape}")

    # Preserve raw counts stored in the file; only synthesize them from X for
    # raw inputs that do not already carry a counts layer.
    if "counts" in adata.layers:
        bioreport.log("Preserved existing `counts` layer from loaded adata")
    else:
        adata.layers["counts"] = adata.X.copy()
        bioreport.log("Initialized `counts` layer from adata.X")
    return adata

def save_data(adata: sc.AnnData,
              h5ad_file: str = "adata_unnamed.h5ad",
              root=False) -> None:
    if root:
        adata.write_h5ad(bioreport.file_save(h5ad_file, data_dir = root_dir / "output"))
    else:
        adata.write_h5ad(bioreport.file_save(h5ad_file))
    bioreport.log("\n> Saved adata with shape " + str(adata.shape))


# subsample into adata_subsampled
def subsample(adata: sc.AnnData, subsample_amount: int = 100_000) -> sc.AnnData:
    adata = adata
    if subsample_amount and adata.n_obs > subsample_amount:
        idx = np.random.default_rng(random_seed).choice(adata.n_obs, size=subsample_amount, replace=False)
        adata = adata[idx, :].copy()
        bioreport.log(f"Subsampled into: {adata.shape}")
    # save subsampled data
    save_data(adata, h5ad_file=f"subsampled_{subsample_amount}.h5ad")


 
def qc_ambient_rna(adata: sc.AnnData,
                   low_count_treshold: int = 1_000,
                   frac_treshold: float = 0.05):
    bioreport.log("1.1 Ambient RNA QC", heading=3)
    bioreport.log(f"Checking if a cleanup is necessary, flagging cells with less than {low_count_treshold} counts as potential ambient RNA contamination")

    low_count_cells = np.asarray(adata.X.sum(axis=1)) < low_count_treshold
    if low_count_cells.mean() > frac_treshold:
        bioreport.log("Potential RNA contamination detected, cleaning up the data")
        bioreport.log("no cleanup code exists yet, exiting for now")
        exit()
    else:
        bioreport.log("No significant RNA contamination detected, skipping cleanup")

    return adata

def qc_counts(adata: sc.AnnData,
              mitochondria_cutoff_percentage: float = 8,
              mad_outlier_var: int = 5,
              min_cells_per_gene: int = 20) -> sc.AnnData:
    bioreport.log("1.2 Count QC", heading=3)
    ### QC metrics
 
    #### mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    #### ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    #### hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")


    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt","ribo", "hb"],  # what percentage of counts is mito, ribosomal, hemoglobin
        percent_top=[20, 50, 100],    # what percenteage of counts is top 20, 50, 100 genes
        inplace=True
    )
    bioreport.table(adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb", "pct_counts_in_top_20_genes"]].describe().T)


    ### plot before QC
    sns.displot(adata.obs["total_counts"], bins=100, kde=True)
    bioreport.figure("qc_before_total_counts_hist")

    sns.displot(adata.obs["pct_counts_mt"], bins=100, kde=True)
    bioreport.figure("qc_before_pct_counts_mt_hist", sameline=True)

    sc.pl.violin(adata, "pct_counts_mt", show=False)
    bioreport.figure("qc_before_pct_counts_mt_violin", sameline=True)

    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
    bioreport.figure("qc_before_scatter_genes_total_mt", sameline=True)

    ### Get outliers with MAD

    def outlier_MAD(adata, metric: str, nmads: int) -> pd.Series:
        M = adata.obs[metric]
        outlier = (M < np.median(M) - nmads * stats.median_abs_deviation(M)) | (
            np.median(M) + nmads * stats.median_abs_deviation(M) < M
        )
        return outlier


    adata.obs["outlier"] = (
        outlier_MAD(adata, "log1p_total_counts", mad_outlier_var)
        | outlier_MAD(adata, "log1p_n_genes_by_counts", mad_outlier_var)
        | outlier_MAD(adata, "pct_counts_in_top_20_genes", mad_outlier_var)
    )
    bioreport.log(f"Number of count outliers with `MAD > {mad_outlier_var}`: " + str(adata.obs["outlier"].value_counts()[True]))

    adata.obs["mt_outlier"] = outlier_MAD(adata, "pct_counts_mt", mad_outlier_var) | (
        adata.obs["pct_counts_mt"] > mitochondria_cutoff_percentage
    )
    bioreport.log(f"Number of mitochondrial outliers with `MAD > {mad_outlier_var}` or `pct_counts_mt > {mitochondria_cutoff_percentage}`: " + str(adata.obs["mt_outlier"].value_counts()[True]))
   
    ### Filter based on counts
    before_filter_counts = adata.shape[0]  # statistics
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)]
    bioreport.log(f"Filtered out outliers: {before_filter_counts - adata.shape[0]} ({1 - adata.shape[0]/before_filter_counts:.2%})")

    ### Filter genes with too few detected cells before HVG selection
    genes_before_filter = adata.n_vars # statistics
    bioreport.log(f"Filtering out genes expressed in less than {min_cells_per_gene} cells")
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    bioreport.log(f"Removed {adata.n_vars} genes ({1 - adata.n_vars / genes_before_filter:.2%})")

    ### plot after count QC
    bioreport.log("After count QC:", heading=4)
    sns.displot(adata.obs["total_counts"], bins=100, kde=True)
    bioreport.figure(f"qc_after_total_counts_MAD={mad_outlier_var}", sameline=True)

    sns.displot(adata.obs["pct_counts_mt"], bins=100, kde=True)
    bioreport.figure(f"qc_after_pct_counts_mt_MAD={mad_outlier_var}", sameline=True)

    sc.pl.violin(adata, "pct_counts_mt", show=False)
    bioreport.figure(f"qc_after_pct_counts_mt_violin_MAD={mad_outlier_var}", sameline=True)

    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
    bioreport.figure(f"qc_after_scatter_genes_total_mt_MAD={mad_outlier_var}", sameline=True)

    return adata

def qc_doublets(adata: sc.AnnData,
                expected_doublet_rate: float = 0.05) -> sc.AnnData:
    bioreport.log("1.3 Doublet QC", heading=3)
    n_obs_before_doublet_qc = adata.n_obs  # statistics
    ### prepare input
    counts_csr = adata.layers["counts"].tocsr() if not sp.issparse(adata.layers["counts"]) else adata.layers["counts"].tocsr()
    scrublet = scr.Scrublet(counts_matrix=counts_csr, expected_doublet_rate=expected_doublet_rate, random_state=random_seed)

    ### run scrublet
    bioreport.log(f"Running scrublet with expected doublet rate {expected_doublet_rate}")
    doublet_scores, predicted_doublets = scrublet.scrub_doublets()

    ### save scores, calls, histogram
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["doublet_call"] = pd.Categorical(predicted_doublets)
    adata.uns["scrublet_sim_scores"] = scrublet.doublet_scores_sim_

    scrublet.plot_histogram()
    bioreport.figure("qc_scrublet_histogram")

    ### Filter if necessary
    if adata.obs["doublet_call"].nunique() > 1:
        bioreport.log(f"There are {adata.obs['doublet_call'].nunique()} predicted doublets, applying filter")
        adata = adata[~adata.obs['doublet_call'].astype(bool).to_numpy(), :]
    else:
        bioreport.log("No doublets predicted")
    bioreport.log(f"After doublet step: {adata.shape} (removed {1 - adata.n_obs / n_obs_before_doublet_qc :.2%} cells)")

    return adata

def qc_gene_exclusion(adata: sc.AnnData,
    excluded_gene_patterns = {
        "mt_genes": [r"MT-"],
        "ribo_genes": [r"RPS", r"RPL"],
        "hb_genes": [r"^HBA|^HBB"],
        "stress_genes": ["G0S2", r"^JUN", r"^FOS", "CDKN1A"],
        "IFN_genes": [r"^IRF", r"^IFN"],
        "ccl_genes": [r"CCL"],
        "MHC_genes": [r"^HLA-[A-Z]"],
        "hist_genes": [r"HIST"],
        "comp_genes": [r"^C1Q"],
        "ig_genes": [r"^IGH", r"^IGK", r"^IGL"],
        "miscelaneous_bad_features": ["ACTB"]
        #"cc_genes": [r"^MCM", r"^CDC", r"^CDK"],
        # "cc_genes": [
        #     # Source: https://www.nature.com/articles/s41592-025-02793-1#Sec30 
        #     # https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-025-02793-1/MediaObjects/41592_2025_2793_MOESM5_ESM.xlsx
        #     # - top 50 genes from the gene expression programmes CellCycle-G2M, CellCycle-S, CellCycle-Late-S
        #     "CDC20","CCNB2","DLGAP5","PLK1","UBE2C","CCNB1","HMMR","BIRC5","CENPF","ASPM","TPX2","TOP2A","TROAP","KIF14","CENPA","CDKN3","MKI67","KIF20A","CEP55","CENPE","HMGB2","GTSE1","PTTG1","NEK2","PIMREG","KIF2C","AURKB","FAM64A","CDCA8","CDCA3","KIF4A","KIF23","DEPDC1","NUF2","AURKA","CCNA2","BUB1","DEPDC1B","CDCA2","STMN1","CKAP2L", #"CDK1","NUSAP1","TUBA1B","HMGB3","TTK","KNL1","KIF15","HJURP","KIFC1","KPNA2","PRR11","NCAPG","CKS2","CKS1B","CENPW","TUBB","CDC25C","KIF11","HMGN2","TUBA1C","SGO2","SMC4","FOXM1","CIT","ARL6IP1","SHCBP1","SGO1","TUBB4B","PBK","PSRC1","SAPCD2","H2AFZ","KNSTRN","HN1","UBE2S","H4C3","PRC1","JPT1","KIF18B","MXD3","SPC25","SPC24","ANLN","RACGAP1","ECT2","ARHGAP11A","AB_CD71","BUB1B","MAD2L1","PIF1","SPAG5","KIAA1524","NCAPH","H2AZ1","ANP32E","H2AFV","CCNF","SMC2","PCLAF","KIF20B","CENPN","KIF18A","SKA3","NDC80","PARPBP","POC1A","FAM72B","RAD21","HMGB1","FAM72C","ESPL1","DIAPH3","PHF19","OIP5","APOBEC3B","SKA1","SPDL1","INCENP","ZWINT","NUCKS1","MTFR2","FAM72D","NCAPD2","CIP2A","CKAP5","NEIL3","FBXO43","CKAP2","KIF22","MELK","CENPM","TRIP13","TMPO","GPSM2","H2AFX","SKA2","ARHGEF39","CALM3","DTYMK","LMNB1","UBE2T","GAS2L3","BCL2L12","STIL","TYMS","TK1","H2AX","PLK4","RPL39L","SMTN","FAM83D","CDC25B","DCTN3","POLQ","LSM5","GGH","HNRNPA2B1","CALM2","HDGF","TICRR","AB_CD61","TRAIP","AB_CD278","TACC3","RAN","MND1","ANP32B","AB_CD9","H2AC13","SPTBN1","BORA","MZT1","RDM1","C14orf80","H3C8","CDCA5","NCAPG2","ASF1B","TBXT","GSG2","HP1BP3","LBR","GINS1","MYBL2","BRCA2","EMC9","ARHGAP11B","DEK","GEN1","DDX39A","TPM4","H2AC14","FANCI","UBALD2","BUB3","NUDT1","REEP4","LGALS1","ERCC6L",
        #     "GINS2","TK1","TYMS","DTL","UHRF1","PCLAF","CLSPN","CDC45","MCM10","FAM111B","CDT1","CDC6","E2F1","PCNA","RRM2","MCM4","ASF1B","MCM2","CENPU","STMN1","MCM7","DUT","PKMYT1","MYBL2","MCM5","ZWINT","CHEK1","DHFR","HELLS","FEN1","MCM3","BRCA1","CDCA7","EXO1","MCM6","ORC6","PAQR4","TUBA1B","C16orf59","RAD51AP1","DSCC1", #"ORC1","CDCA5","RAD51","CENPM","ESCO2","MND1","CHAF1A","TCF19","RAD54L","AB_CD71","CCNE2","RMI2","E2F8","TEDC2","WDR34","NASP","CDC25A","GMNN","ATAD2","E2F2","SMC2","DMC1","CENPW","RBBP8","DEK","ATAD5","TMEM106C","TIMELESS","H2AFZ","CTNNAL1","POLE2","WDR76","H2AZ1","FANCI","SLBP","BRCA2","HMGB2","GGH","RFC4","RRM1","DIAPH3","WDHD1","LIG1","CENPH","DNAJC9","TUBB","DNMT1","HMGN2","GINS1","E2F7","UBE2T","PTMA","UNG","CHAF1B","SGO1","BRIP1","AUNIP","TYMSOS","RPA3","MELK","MKI67","FBXO5","CARHSP1","NUDT1","CCNE1","SPC25","SHCBP1","POLQ","LMNB1","XRCC2","RFC2","MSH6","SKA3","CENPN","HMGB1","RIBC2","PSMC3IP","EZH2","CDCA7L","SNRNP25","C19orf48","POC1A","CDK1","POLD3","RANBP1","CENPK","MAD2L1","RECQL4","SIVA1","USP1","ACOT7","RFC5","ZNF367","DNA2","MTHFD1","GINS3","IFI27L1","H2AFV","GAPDH","IDH2","H4C3","HNRNPAB","SMC3","DSN1","NCAPH","NCAPG2","RFC3","CENPP","CDK2","CENPX","RNASEH2A","BIRC5","SPC24","NEIL3","H2AFY","ZGRF1","KIAA1524","OIP5","ANP32B","CEP55","CDK4","FANCG","SLC29A1","GINS4","KNTC1","SLC43A3","SLC25A5","PXMP2","GGCT","TFDP1","POLA1","HNRNPF","FABP5","CBX5","POLD1","MTFR2","KIFC1","TP73","SMC1A","FANCA","PAICS","PPIA","NSD2","C21orf58","PHGDH","STIL","RAN","SNRPB","TOPBP1","PLK4","DTYMK","STRIP2","FOXM1","BARD1","HMGA1","EME1","MSH2","SYCE2","MMS22L",
        #     "HIST1H1B","HIST1H4C","HIST1H3G","HIST1H3B","RRM2","SPC25","CDCA5","NUSAP1","HIST1H2AJ","HJURP","TOP2A","HIST1H3F","TYMS","KIF15","PKMYT1","ESCO2","CDK1","E2F8","KIFC1","UBE2C","DIAPH3","CKAP2L","MKI67","KIF23","TUBB","NCAPG","HIST1H2AH","PBK","GTSE1","AURKB","H4C3","ANLN","MELK","KNL1","CDCA2","C2orf48","STMN1","KIF11","HIST1H2AL","CDCA3","SPC24", #"TUBA1B","MND1","CCNA2","RAD51AP1","KIF2C","PCLAF","SKA3","FBXO5","TK1","ZWINT","HIST1H3D","CDCA8","SKA1","BUB1B","POLQ","KIF4A","H2AFX","SHCBP1","KIF18B","ASF1B","ASPM","HIST1H2AI","MYBL2","TPX2","HIST1H1C","ESPL1","HIST1H1D","MXD3","HIST1H3C","E2F7","CDC25C","HIST1H2AM","NCAPH","FOXM1","FANCI","FAM111B","BIRC5","HIST2H2AC","RAD54L","UBE2T","FBXO43","HIST1H2BH","STIL","NEIL3","CENPM","DLGAP5","EZH2","PCNA","HIST1H1E","HMMR","CLSPN","CENPF","NDC80","CKS1B","SAC3D1","CDT1","CENPU","HMGB1","CKS2","WDR34","TTK","ATAD2","TICRR","SMC4","CENPA","HMGB2","CENPE","PRC1","KIF14","C21orf58","BUB1","AB_CD38-2","DEPDC1","HMGN2","HASPIN","SGO1","RACGAP1","DSCC1","TMPO","GGH","HIST1H4H","TRIP13","ARHGAP11A","HIST1H2BM","HIST1H2AG","HIST2H4B","HIST1H2BO","CDC45","NCAPG2","ECT2","KIF22","H2AFZ","DUT","UHRF1","KIAA1524","EXO1","CENPW","SMC2","HIST1H2BL","SPAG5","POC1A","CIT","CIP2A","TUBB4B","RAD51","E2F2","CCNF","SGO2","RRM1","HIST1H4D","KIF24","BRCA1","TMEM106C","HIST1H2BF","LMNB1","KIF20A","DEPDC1B","CDKN2C","H2AFV","CENPN","MTFR2","APOBEC3B","ERCC6L","BRCA2","CEP55","TACC3","AURKA","HIST1H2AB","CCDC150","FEN1","EME1","GSG2","HIST2H2AB","PLK4","WDR62","NUF2","TCF19","INCENP","CCDC34","HIST1H1A","C1QL1","HIST1H4I","C16orf59","ORC6","H2AZ1","MCM10","DHFR","CDCA4","CARHSP1","TYMSOS","H1-1","HIST1H2BI","H1-5","CDC25A","RECQL4","C17orf53","HIST1H2BB","ARHGEF39","DNAJC9",
        # ],
                                }) -> sc.AnnData:

    bioreport.log("1.4 Gene exclusion", heading=3)


    bioreport.log("Excluding genes matching the following patterns:")
    bioreport.table(excluded_gene_patterns)

    ### join all the lists together
    all_excluded_gene_patterns = [pattern for patterns in excluded_gene_patterns.values() for pattern in patterns]

    ### create the mask with regex
    excluded_gene_mask = pd.Series(False, index=adata.var_names, dtype=bool)
    for pattern in all_excluded_gene_patterns:
        excluded_gene_mask |= adata.var_names.str.upper().str.contains(pattern.upper(), regex=True)

    ### Filter
    bioreport.log(f"Removing {int(excluded_gene_mask.sum())} genes ({excluded_gene_mask.sum() / len(adata.var_names) * 100:.2f}%)")
    adata = adata[:, ~excluded_gene_mask.to_numpy()].copy()
    bioreport.log(f"Shape after exclusion panel filtering: {adata.shape}")

    return adata

def qc_full(adata: sc.AnnData,
            low_count_treshold: int = 1_000,
            frac_treshold: float = 0.05,
            mitochondria_cutoff_percentage: float = 8,
            mad_outlier_var: int = 5,
            min_cells_per_gene: int = 20,
            expected_doublet_rate: float = 0.05,
            excluded_gene_patterns: dict[str, list[str]] | None = None):
    bioreport.log("1. Quality Control", heading=2)
    adata = qc_ambient_rna(adata, low_count_treshold, frac_treshold)
    adata = qc_counts(adata, mitochondria_cutoff_percentage, mad_outlier_var, min_cells_per_gene)
    adata = qc_doublets(adata, expected_doublet_rate)
    if excluded_gene_patterns is None:
        adata = qc_gene_exclusion(adata)
    else:
        adata = qc_gene_exclusion(adata, excluded_gene_patterns)

    return adata

def normalise(adata: sc.AnnData) -> sc.AnnData:
    bioreport.log("2. Normalisation", heading=2)

    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
    adata.layers["normalised"] = adata.X.copy()
    bioreport.log("Normalised in reference to total counts")

    sc.pp.log1p(adata)
    adata.X = adata.X.astype(np.float32)  # for better storage
    adata.layers["log1p"] = adata.X.copy()
    bioreport.log("Transformed with log1p")

    log1p_library_size = np.asarray(adata.layers["log1p"].sum(axis=1)).ravel()

    sns.histplot(adata.obs["total_counts"], bins=100, kde=False)
    bioreport.figure("normalisation_total_counts_hist")
    sns.histplot(log1p_library_size, bins=100, kde=False)
    bioreport.figure("normalisation_log1p_hist", sameline=True)

    return adata

# labeling, CC analysis

# 3. Feature selection

def selection_cells(adata: sc.AnnData) -> sc.AnnData:
    return adata

def label_genes(adata: sc.AnnData) -> sc.AnnData:
    label_file = root_dir / "data" / "protein-coding_gene.txt"
    label_df = pd.read_csv(label_file, sep="\t", dtype=str).fillna("")

    if "ensembl_gene_id" not in label_df.columns or "symbol" not in label_df.columns:
        raise KeyError("`protein-coding_gene.txt` must contain `ensembl_gene_id` and `symbol` columns")

    original_var_names = pd.Index(adata.var_names.astype(str))
    stripped_var_names = original_var_names.str.replace(r"\.\d+$", "", regex=True)
    stripped_var_series = pd.Series(stripped_var_names, index=adata.var_names)

    label_df = label_df.drop_duplicates(subset="ensembl_gene_id", keep="first")
    symbol_map_by_ensembl = label_df.set_index("ensembl_gene_id")["symbol"]
    name_map_by_ensembl = (
        label_df.set_index("ensembl_gene_id")["name"] if "name" in label_df.columns else pd.Series(dtype=str)
    )
    name_map_by_symbol = (
        label_df.drop_duplicates(subset="symbol", keep="first").set_index("symbol")["name"]
        if "name" in label_df.columns
        else pd.Series(dtype=str)
    )

    adata.var["gene_id_original"] = original_var_names.to_numpy()
    adata.var["gene_id_stripped"] = stripped_var_names.to_numpy()

    gene_symbols = stripped_var_series.map(symbol_map_by_ensembl)
    symbol_mask = stripped_var_series.isin(name_map_by_symbol.index)
    gene_symbols = gene_symbols.where(gene_symbols.notna(), stripped_var_series.where(symbol_mask, ""))

    gene_names = stripped_var_series.map(name_map_by_ensembl)
    gene_names = gene_names.where(gene_names.notna(), gene_symbols.map(name_map_by_symbol))

    adata.var["gene_symbol"] = gene_symbols.fillna("").to_numpy()
    adata.var["gene_name_human_readable"] = gene_names.fillna("").to_numpy()
    adata.var["gene_label_plot"] = [
        f"{symbol} ({human_name[:25]})" if symbol and human_name
        else symbol if symbol
        else original_id
        for original_id, symbol, human_name in zip(
            original_var_names,
            adata.var["gene_symbol"],
            adata.var["gene_name_human_readable"],
            strict=False,
        )
    ]

    n_relabeled = int((adata.var["gene_symbol"] != "").sum())
    bioreport.log(f"Labeled {n_relabeled} genes from {label_file.name}")
    return adata


def plot_var_names(adata: sc.AnnData) -> pd.Index:
    if "gene_label_plot" not in adata.var.columns:
        return pd.Index(adata.var_names.astype(str))

    plot_names = [
        plot_name if plot_name else original_id
        for original_id, plot_name in zip(
            adata.var_names.astype(str),
            adata.var["gene_label_plot"].astype(str),
            strict=False,
        )
    ]
    unique_plot_names = pd.Index(
        [
            plot_name if plot_names[:i].count(plot_name) == 0 else f"{plot_name}-{plot_names[:i + 1].count(plot_name)}"
            for i, plot_name in enumerate(plot_names)
        ]
    )
    return unique_plot_names


def label_states(adata: sc.AnnData,
                 marker_csv_name: str = 'state_markers.csv'  # src https://www.biorxiv.org/content/10.1101/2024.07.17.603780v1.full
                 ) -> sc.AnnData:
    marker_file = root_dir / "data" / marker_csv_name
    marker_df = pd.read_csv(marker_file, sep=",", quotechar="'", skipinitialspace=True, dtype=str).fillna("")
    dataset_genes = set(adata.var_names)

    bioreport.log("Labeling states based on markers:", heading=3)
    marker_table = pd.DataFrame(
        {
            state_name: pd.Series(state_rows["gene"].tolist())
            for state_name, state_rows in marker_df.groupby("marker", sort=False)
        }
    ).fillna("")
    bioreport.table(marker_table.to_dict(orient="records"))

    score_columns: dict[str, str] = {}
    score_layer = "log1p"
    all_state_genes: set[str] = set()
    all_matched_genes: set[str] = set()

    for state_name, state_rows in marker_df.groupby("marker", sort=False):
        state_genes = state_rows["gene"].tolist()
        matched_features = [gene for gene in state_genes if gene in dataset_genes]
        all_state_genes.update(state_genes)
        all_matched_genes.update(matched_features)

        col_suffix = re.sub(r"[^A-Za-z0-9]+", "_", state_name).strip("_").lower()
        score_col = f"state_score_{col_suffix}"
        score_columns[state_name] = score_col

        sc.tl.score_genes(
            adata,
            gene_list=matched_features,
            score_name=score_col,
            ctrl_size=max(1, min(50, len(matched_features))),
            random_state=random_seed,
            use_raw=False,
            layer=score_layer,
        )

    overall_missing_genes = sorted(all_state_genes.difference(all_matched_genes))
    n_matched_markers = int(len(all_matched_genes))
    n_missing_markers = int(len(overall_missing_genes))
    bioreport.log(
        f"State marker genes in dataset: {n_matched_markers} present, {n_missing_markers} missing")
    if overall_missing_genes:
        bioreport.log(f"Omitted state markers: {', '.join(overall_missing_genes)}")

    adata.uns["state_marker_coverage"] = {
        "n_present": n_matched_markers,
        "n_missing": n_missing_markers,
    }

    score_df = pd.DataFrame({state: adata.obs[col] for state, col in score_columns.items()}, index=adata.obs_names)
    adata.obs["state_label"] = score_df.idxmax(axis=1)
    adata.obs["state_score_max"] = score_df.max(axis=1)

    sorted_scores = np.sort(score_df.to_numpy(copy=True), axis=1)
    if score_df.shape[1] >= 2:
        adata.obs["state_score_margin"] = sorted_scores[:, -1] - sorted_scores[:, -2]
    else:
        adata.obs["state_score_margin"] = np.nan

    bioreport.log("Assigned state labels:", heading=3)
    state_counts = adata.obs["state_label"].value_counts().rename_axis("state").reset_index(name="n_cells")
    state_counts["pct_cells"] = (100 * state_counts["n_cells"] / adata.n_obs).round(1)
    bioreport.table(state_counts.to_dict(orient="records"), maxrows=20)
    return adata

def selection_hvg(adata: sc.AnnData,
                  n_top_genes: int = 4000,
                  old_seurat: bool = True,
                  whitelist_genes: list[str] | None = None) -> sc.AnnData:
    bioreport.log("3.3 Highly Variable Gene selection", heading=3)
    ### using seurat as seurat_v3 with original counts raises an environment error
    if old_seurat:
        sc.pp.highly_variable_genes(adata, layer='log1p', flavor='seurat', n_top_genes=4000, inplace=True)
        y_axis = "dispersions_norm"
    else:
        sc.pp.highly_variable_genes(adata, layer='counts', flavor='seurat_v3', n_top_genes=n_top_genes, inplace=True)
        y_axis = "variances_norm"

    if whitelist_genes:
        whitelist_set = {gene for gene in whitelist_genes if gene}
        whitelist_mask = adata.var_names.isin(whitelist_set)
        if "gene_symbol" in adata.var.columns:
            whitelist_mask = whitelist_mask | adata.var["gene_symbol"].astype(str).isin(whitelist_set).to_numpy()
        n_whitelist_matched = int(whitelist_mask.sum())
        adata.var["highly_variable"] = adata.var["highly_variable"] | whitelist_mask
        bioreport.log(
            f"HVGs selected: {int(adata.var['highly_variable'].sum())} "
            f"(including {n_whitelist_matched} whitelisted genes)"
        )
    else:
        bioreport.log(f"HVGs selected: {int(adata.var['highly_variable'].sum())}")

    #sc.pl.highly_variable_genes(adata, show=False)
    sns.scatterplot(data=adata.var, x="means", y=y_axis, hue="highly_variable", s=5)
    bioreport.figure(f"hvg_before_{'seurat' if old_seurat else 'seurat_v3'}_{n_top_genes}", sameline=True)

    ### Filter
    adata = adata[:, adata.var.highly_variable]
    
    #sc.pl.highly_variable_genes(adata, show=False)
    sns.scatterplot(data=adata.var, x="means", y=y_axis, color="orange", s=5)
    bioreport.figure(f"hvg_after_{'seurat' if old_seurat else 'seurat_v3'}_{n_top_genes}", sameline=True)
    return adata


def display_top_genes(adata: sc.AnnData, top_n: int = 20):
# plot top 20 genes for sanity check
    bioreport.log(f"Top {top_n} expressed genes after selection:", heading=3)
    adata_plot = adata.copy()
    adata_plot.var_names = plot_var_names(adata_plot)
    sc.pl.highest_expr_genes(adata_plot, n_top=top_n, show=False)
    bioreport.figure(f"top{top_n}_expr_genes")


def marker_genes_from_csv(marker_csv_name: str) -> list[str]:
    marker_file = root_dir / "data" / marker_csv_name
    marker_df = pd.read_csv(marker_file, sep=",", quotechar="'", skipinitialspace=True, dtype=str).fillna("")
    return marker_df["gene"].tolist()


def marker_gene_groups_from_csv(marker_csv_name: str) -> dict[str, list[str]]:
    marker_file = root_dir / "data" / marker_csv_name
    marker_df = pd.read_csv(marker_file, sep=",", quotechar="'", skipinitialspace=True, dtype=str).fillna("")
    return {
        state_name: state_rows["gene"].tolist()
        for state_name, state_rows in marker_df.groupby("marker", sort=False)
    }


def feature_selection_full(adata: sc.AnnData,
                           n_top_genes: int = 4000,
                           old_seurat: bool = True,
                           include_markers: str | None = None,
                           whitelist_genes: list[str] | None = None,
                           display_top_n_genes: int = 20) -> sc.AnnData:
    bioreport.log("3. Feature selection", heading=2)
    if include_markers is not None:
        marker_genes = marker_genes_from_csv(include_markers)
        whitelist_genes = marker_genes if whitelist_genes is None else [*whitelist_genes, *marker_genes]

    adata = selection_cells(adata)
    adata = selection_hvg(
        adata,
        n_top_genes=n_top_genes,
        old_seurat=old_seurat,
        whitelist_genes=whitelist_genes,
    )
    display_top_genes(adata, top_n=display_top_n_genes)

    return adata





# region 4. Dimensional reduction (=Embedding) 
def pca(adata: sc.AnnData,
        n_pcs: int = 30,
        var_scale_cutoff: int = 10) -> sc.AnnData:
    bioreport.log("4. Dimensional reduction", heading=2)

    ## PCA (linear reduction)
    bioreport.log("Principal component analysis (PCA)", heading=3)
    sc.pp.scale(adata, max_value=var_scale_cutoff)  # scale to variance
    bioreport.log(f"Scaled data with max variance cutoff {var_scale_cutoff}")

    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", random_state=random_seed)
    bioreport.log(f"Calculated PCA with {n_pcs} components")

    ## examine PCA - plot variance ratio
    sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs, log=True, show=False)
    bioreport.figure(f"pca_variance_ratio_n={n_pcs}")

    ## examine PCA - plot scatter with different colours
    pca_color_cols = ["cell_subtype_3", "cell_subset", "disease_type", "meta_tissue_type", "study_name", "total_counts"]
    sc.pl.pca_scatter(adata, color=pca_color_cols, ncols=2, show=False)
    bioreport.figure(f"pca_scatter_n={n_pcs}", dpi=100, sameline=True)

    return adata

def umap(adata: sc.AnnData,
         n_pcs: int = 30,
         n_neighbors: int = 15,
         min_dist: float = 0.3,
         spread: float = 1.0) -> sc.AnnData:
    bioreport.log("UMAP", heading=3)

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_seed)
    bioreport.log(f"Calculated nearest {n_neighbors} neighbours using {n_pcs} PCs")
    sc.tl.umap(adata, min_dist=min_dist, spread=spread, random_state=random_seed)
    bioreport.log(f"Calculated UMAP with min_dist {min_dist} and spread {spread}")
    return adata


def label_clusters(adata: sc.AnnData,
                   leiden_resolution: float = 0.5) -> sc.AnnData:
    bioreport.log("5. Clustering", heading=2)

    leiden_key = LEIDEN_COL
    sc.tl.leiden(adata, resolution=leiden_resolution, key_added=leiden_key)
    n_clusters = adata.obs[leiden_key].nunique()
    bioreport.log(f"Detected {n_clusters} clusters with leiden at resolution {leiden_resolution}")
    return adata


def umap_plot(adata: sc.AnnData,
              umap_color_cols: list[str] | None = None,
              name: str = "default") -> None:
    default_leiden = LEIDEN_COL if LEIDEN_COL in adata.obs.columns else None
    if default_leiden is None:
        leiden_cols = [col for col in adata.obs.columns if col.startswith("leiden_res_")]
        default_leiden = leiden_cols[0] if leiden_cols else None
    if umap_color_cols is None:
        umap_color_cols = [
            col
            for col in [
                default_leiden,
                "cell_subtype_3",
                "disease_type",
                "meta_tissue_type",
                "state_score_max",
                "state_score_margin",
            ]
            if col is not None and col in adata.obs.columns
        ]

    sc.pl.umap(adata, color=umap_color_cols, ncols=2, wspace=0.35, show=False)
    bioreport.figure(f"umap_plot_{name}")


def pca_umap(adata: sc.AnnData,
             n_pcs: int = 30,
             n_neighbors: int = 15,
             min_dist: float = 0.3,
             spread: float = 1.0) -> sc.AnnData:
    adata = pca(adata, n_pcs=n_pcs)
    adata = umap(
        adata,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        spread=spread,
    )
    return adata


def umap_basic(adata: sc.AnnData,
               n_pcs: int = 30,
               n_neighbors: int = 15,
               min_dist: float = 0.3,
               leiden_resolution: float = 0.5,
               spread: float = 1.0,
               save: bool = True) -> sc.AnnData:
    adata = pca_umap(
        adata,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        spread=spread,
    )
    adata = label_clusters(adata, leiden_resolution=leiden_resolution)
    umap_plot(
        adata,
        name=f"n={n_neighbors}_mindist={min_dist}_res={leiden_resolution}",
    )
    if save:
        save_data(adata, h5ad_file=f"adata_umap_clustering_n={n_neighbors}_mindist={min_dist}_res={leiden_resolution}.h5ad")
    return adata


# test

def findallmarkers(adata: sc.AnnData,
                   cluster_col: str | None = LEIDEN_COL,
                   leiden_res: float | None = None,
                   n_workers: int | None = None,
                   n_markers_per_cluster: int = 5) -> sc.AnnData:
    import tempfile
    import os
    from scipy.io import mmwrite
    output_basename = "findallmarkers_temp"

    if leiden_res is not None:
        cluster_col = LEIDEN_COL
    elif cluster_col is None:
        raise ValueError("Either `cluster_col` or `leiden_res` must be provided")

    if n_workers is None:
        n_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    n_workers = max(1, int(n_workers))
    if n_markers_per_cluster < 1:
        raise ValueError("`n_markers_per_cluster` must be at least 1")

    bioreport.log("Finding significant genes with FindAllMarkers (via Seurat in R)", heading=3)
    bioreport.log(f"FindAllMarkers will use up to {n_workers} R worker(s)")
    output_csv = Path(root_dir) / "output" / f"{output_basename}_markers.csv"

    counts_matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    counts_matrix = counts_matrix if sp.issparse(counts_matrix) else sp.csr_matrix(counts_matrix)
    counts_matrix = counts_matrix.transpose().tocoo()

    r_script_content = '''
    args <- commandArgs(trailingOnly = TRUE)
    matrix_file <- args[1]
    genes_file <- args[2]
    cells_file <- args[3]
    metadata_file <- args[4]
    output_file <- args[5]
    cluster_col <- args[6]
    workers <- as.integer(args[7])
    suppressPackageStartupMessages({
      library(Matrix)
      library(Seurat)
    })
    counts <- readMM(matrix_file)
    genes <- read.csv(genes_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
    cells <- read.csv(cells_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
    metadata <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)
    rownames(counts) <- genes
    colnames(counts) <- cells
    seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    Idents(seurat_obj) <- cluster_col
    if (requireNamespace("future", quietly = TRUE) && workers > 1) {
      options(future.globals.maxSize = 32 * 1024^3)
      future::plan("multicore", workers = workers)
      on.exit(future::plan("sequential"), add = TRUE)
      message("Using future::multicore with ", workers, " workers for FindAllMarkers")
    } else {
      message("Running FindAllMarkers serially")
    }
    markers <- FindAllMarkers(
      seurat_obj,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25
    )
    if (nrow(markers) == 0) {
      markers <- data.frame(
        p_val = numeric(),
        avg_log2FC = numeric(),
        pct.1 = numeric(),
        pct.2 = numeric(),
        p_val_adj = numeric(),
        cluster = character(),
        gene = character(),
        stringsAsFactors = FALSE
      )
    } else if (!"gene" %in% colnames(markers)) {
      markers$gene <- rownames(markers)
    }
    write.csv(markers, file = output_file, row.names = FALSE)
    '''

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        matrix_file = temp_dir / "counts.mtx"
        genes_file = temp_dir / "genes.csv"
        cells_file = temp_dir / "cells.csv"
        metadata_file = temp_dir / "metadata.csv"
        r_script_path = temp_dir / "run_findallmarkers.R"

        mmwrite(str(matrix_file), counts_matrix)
        pd.Series(adata.var_names).to_csv(genes_file, index=False, header=False)
        pd.Series(adata.obs_names).to_csv(cells_file, index=False, header=False)
        adata.obs[[cluster_col]].copy().to_csv(metadata_file)
        r_script_path.write_text(r_script_content)

        subprocess.run([
            "Rscript",
            str(r_script_path),
            str(matrix_file),
            str(genes_file),
            str(cells_file),
            str(metadata_file),
            str(output_csv),
            cluster_col,
            str(n_workers),
        ], check=True)

    try:
        marker_df = pd.read_csv(output_csv)
    except pd.errors.EmptyDataError:
        marker_df = pd.DataFrame(columns=["gene", "cluster"])

    if "gene" not in marker_df.columns:
        marker_df["gene"] = pd.Series(dtype=str)
    if "cluster" not in marker_df.columns:
        marker_df["cluster"] = pd.Series(dtype=str)

    ranking_cols = []
    ascending = []
    for col, asc in [("cluster", True), ("avg_log2FC", False), ("pct.1", False), ("p_val_adj", True)]:
        if col in marker_df.columns:
            ranking_cols.append(col)
            ascending.append(asc)

    if not marker_df.empty and "cluster" in marker_df.columns:
        if ranking_cols:
            marker_df = marker_df.sort_values(ranking_cols, ascending=ascending)
        marker_df = marker_df.groupby("cluster", sort=False, group_keys=False).head(n_markers_per_cluster).copy()

    marker_genes = set(marker_df["gene"].dropna())
    adata.uns["findallmarkers"] = marker_df
    adata.uns["findallmarkers_cluster_col"] = cluster_col
    adata.var["findallmarkers_selected"] = adata.var_names.isin(marker_genes)
    bioreport.log(
        f"Ran FindAllMarkers in R, keeping up to {n_markers_per_cluster} markers per cluster "
        f"and selected {int(adata.var['findallmarkers_selected'].sum())} genes"
    )

    return adata

def dotplot(adata: sc.AnnData,
            cluster_col: str | None = None,
            genes_bool_col: str = "findallmarkers_selected",
            standard_scale: str | None = "var",
            gene_groups: dict[str, list[str]] | None = None,
            swap_axes: bool = False,
            figure_name: str | None = None) -> None:

    if cluster_col is not None and cluster_col not in adata.obs:
        raise KeyError(f"Column `{cluster_col}` not found in adata.obs")

    if genes_bool_col not in adata.var.columns:
        raise KeyError(f"Column `{genes_bool_col}` not found in adata.var")

    gene_mask = adata.var[genes_bool_col].fillna(False).astype(bool)
    dotplot_gene_ids = adata.var_names[gene_mask].tolist()

    if not dotplot_gene_ids:
        bioreport.log(f"No genes selected in adata.var[`{genes_bool_col}`], skipping dotplot")
        return

    if figure_name is None:
        figure_name = f"dotplot_{cluster_col}_{genes_bool_col}" if cluster_col is not None else f"dotplot_{genes_bool_col}"

    bioreport.log(f"Plotting {len(dotplot_gene_ids)} genes from adata.var[`{genes_bool_col}`]")
    dotplot_layer = "log1p" if "log1p" in adata.layers else None
    adata_plot = adata[:, dotplot_gene_ids].copy()
    adata_plot.var_names = plot_var_names(adata_plot)
    plot_name_by_gene_id = dict(zip(dotplot_gene_ids, adata_plot.var_names.tolist(), strict=False))
    if cluster_col is None:
        adata_plot.obs["_all_cells"] = "all"
        groupby_col = "_all_cells"
    else:
        groupby_col = cluster_col

    var_names = adata_plot.var_names.tolist()
    if gene_groups is not None:
        grouped_var_names = {
            group_name: [plot_name_by_gene_id[gene_id] for gene_id in gene_ids if gene_id in plot_name_by_gene_id]
            for group_name, gene_ids in gene_groups.items()
        }
        var_names = {group_name: genes for group_name, genes in grouped_var_names.items() if genes}

    plot = sc.pl.dotplot(
        adata_plot,
        var_names=var_names,
        groupby=groupby_col,
        use_raw=False,
        layer=dotplot_layer,
        standard_scale=standard_scale,
        swap_axes=swap_axes,
        show=False,
        return_fig=True,
    )
    plot.make_figure()
    bioreport.figure(figure_name)
    return

def dotplot_allmarkers(adata: sc.AnnData,
                       cluster_col: str | None = None,
                       standard_scale: str | None = "var",
                       n_markers_per_cluster: int = 5) -> sc.AnnData:
    bioreport.log("6. FindAllMarkers", heading=2)

    if cluster_col is None:
        cluster_col = adata.uns.get("findallmarkers_cluster_col")
    if cluster_col is None:
        if LEIDEN_COL in adata.obs.columns:
            cluster_col = LEIDEN_COL
        else:
            leiden_cols = [col for col in adata.obs.columns if col.startswith("leiden_res_")]
            cluster_col = leiden_cols[0] if leiden_cols else None

    if cluster_col is not None:
        adata = findallmarkers(
            adata,
            cluster_col=cluster_col,
            n_markers_per_cluster=n_markers_per_cluster,
        )
    dotplot(
        adata,
        cluster_col=cluster_col,
        genes_bool_col="findallmarkers_selected",
        standard_scale=standard_scale,
        figure_name=f"dotplot_{cluster_col}_allmarkers" if cluster_col is not None else "dotplot_allmarkers",
    )
    return adata


def export():
    bioreport.exportpdf()

def test_pipeline_full():
    init("pipeline_elements")
    # adata = subsample(adata, save=True)
    # adata = qc_full(adata)
    # adata = normalise(adata)

    # #save_data(adata, h5ad_file="qc+subsampled_100000.h5ad", root=True)
    
    adata = load_data('qc+subsampled_100000.h5ad')

    adata = label_genes(adata)
    adata = feature_selection_full(adata, include_markers='state_markers.csv')

    adata = umap_basic(
        adata,
        n_pcs=25,
        n_neighbors=10,
        min_dist=0.2,
        leiden_resolution=0.6,
        save=True
    )
    
    #adata = load_data('qc+subsampled+hvg+pca+umap.h5ad')

    adata = label_states(adata, "state_markers.csv") # from https://www.biorxiv.org/content/10.1101/2024.07.17.603780v1.full
    umap_plot(adata, umap_color_cols=["state_label", "state_score_max", "state_score_margin"], name="states")
    
    adata.var["state_markers_selected"] = adata.var_names.isin(marker_genes_from_csv("state_markers.csv"))
    dotplot(adata, cluster_col="state_label", genes_bool_col="state_markers_selected", figure_name="dotplot_state_label_X_state_markers")
    
    dotplot(adata, cluster_col="leiden", genes_bool_col="state_markers_selected", figure_name="dotplot_clustering_X_state_markers")

    # adata = load_data('adata_umap_clustering_n=15_mindist=0.3_res=0.5.h5ad')
    adata = dotplot_allmarkers(adata, n_markers_per_cluster=3, cluster_col="leiden")

    export()

if __name__ == "__main__":
    test_pipeline_full()
