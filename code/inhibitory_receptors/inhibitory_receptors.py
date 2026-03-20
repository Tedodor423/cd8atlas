import os

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import csv
import re
import time, matplotlib
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
import scipy.stats as stats
import statsmodels.formula.api as smf
from collections import defaultdict
from statsmodels.stats.multitest import multipletests
# import PyComplexHeatmap as pch
from pathlib import Path


# SETTINGS ####
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.figsize"] = [5, 4]
plt.rcParams["savefig.bbox"] = "tight"

OUTPUT_NAME = "inhibitory_receptors"

root_dir = Path("/ceph/project/sharmalab/dnimrich/cd8atlas")

input_h5ad = root_dir / "data/" / "huARdb_v2_GEX.CD8.all_genes.h5ad"  # "huARdb_v2_GEX.CD8.hvg4k.h5ad"
markers_csv = root_dir / "code" / "inhibitory_receptors" / "exhaustion_markers.csv"
receptors_csv = root_dir / "code" / "inhibitory_receptors" / "inhibitory_receptor_list.csv"

figures_dir = root_dir / "output" / OUTPUT_NAME / "figures"
data_dir = root_dir / "output" / OUTPUT_NAME / "data"

figures_dir.mkdir(parents=True, exist_ok=True)
data_dir.mkdir(parents=True, exist_ok=True)

random_seed = 423

reference_state_order = [
    "Tn",
    "Tcm",
    "Early Tcm/Tem",
    "GZMK+ Tem",
    "ZNF683+ Teff",
    "GZMK+ Tex",
    "ITGAE+ Tex",
    "XBP1+ Tex",
    "Tpex",
]

exhausted_reference_states = ["Tpex", "GZMK+ Tex", "ITGAE+ Tex", "XBP1+ Tex"]
major_disease_types_for_covariates = ["Healthy", "Infection", "Solid tumor"]
##############


def unique_preserve_order(values: list[str]) -> list[str]:
    seen = set()
    ordered_values = []

    for value in values:
        if value not in seen:
            seen.add(value)
            ordered_values.append(value)

    return ordered_values


def save_placeholder_figure(output_path: Path, message: str) -> None:
    plt.figure(figsize=(6, 4))
    plt.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    plt.axis("off")
    plt.savefig(output_path)
    plt.close()


def split_into_chunks(values: list[str], chunk_size: int) -> list[list[str]]:
    return [values[i : i + chunk_size] for i in range(0, len(values), chunk_size)]


def first_non_null(values: pd.Series):
    non_null_values = values.dropna()
    if non_null_values.empty:
        return np.nan
    return non_null_values.iloc[0]


def clean_sex_column(values: pd.Series) -> pd.Series:
    cleaned_values = values.astype(str).str.strip()
    output = pd.Series(pd.NA, index=values.index, dtype="object")
    output.loc[cleaned_values.str.fullmatch("Male", case=False, na=False)] = "Male"
    output.loc[cleaned_values.str.fullmatch("Female", case=False, na=False)] = "Female"
    return output


def clean_age_column(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values.astype(str).str.strip(), errors="coerce")


def copy_embedding_results(
    adata_source: AnnData,
    adata_target: AnnData,
    obs_cols: list[str],
    uns_keys: list[str],
) -> None:
    for obs_col in obs_cols:
        adata_target.obs[obs_col] = adata_source.obs[obs_col].copy()

    for obsm_key in ["X_pca", "X_umap"]:
        if obsm_key in adata_source.obsm:
            adata_target.obsm[obsm_key] = adata_source.obsm[obsm_key].copy()

    for obsp_key in ["connectivities", "distances"]:
        if obsp_key in adata_source.obsp:
            adata_target.obsp[obsp_key] = adata_source.obsp[obsp_key].copy()

    for uns_key in uns_keys:
        if uns_key in adata_source.uns:
            uns_value = adata_source.uns[uns_key]
            adata_target.uns[uns_key] = uns_value.copy() if hasattr(uns_value, "copy") else uns_value


def save_dotplot(
    adata_plot: AnnData,
    genes: list[str],
    groupby: str,
    output_path: Path,
    categories_order: list[str] | None = None,
    title: str | None = None,
) -> None:
    if not genes:
        save_placeholder_figure(output_path, "No genes available for dotplot.")
        return

    plot = sc.pl.dotplot(
        adata_plot,
        var_names=genes,
        groupby=groupby,
        categories_order=categories_order,
        use_raw=False,
        show=False,
        return_fig=True,
        title=title,
    )
    plot.savefig(output_path)
    plt.close("all")


def save_matrixplot(
    adata_plot: AnnData,
    genes: list[str],
    groupby: str,
    output_path: Path,
    categories_order: list[str] | None = None,
    title: str | None = None,
    standard_scale: str | None = "var",
) -> None:
    if not genes:
        save_placeholder_figure(output_path, "No genes available for matrixplot.")
        return

    plot = sc.pl.matrixplot(
        adata_plot,
        var_names=genes,
        groupby=groupby,
        categories_order=categories_order,
        use_raw=False,
        show=False,
        return_fig=True,
        title=title,
        standard_scale=standard_scale,
        cmap="viridis",
    )
    plot.savefig(output_path)
    plt.close("all")


def save_heatmap(
    adata_plot: AnnData,
    genes: list[str],
    groupby: str,
    output_path: Path,
    title: str | None = None,
) -> None:
    if not genes:
        save_placeholder_figure(output_path, "No genes available for heatmap.")
        return

    sc.pl.heatmap(
        adata_plot,
        var_names=genes,
        groupby=groupby,
        use_raw=False,
        show=False,
        swap_axes=True,
        standard_scale="var",
        cmap="viridis",
        dendrogram=False,
    )
    if title is not None:
        plt.suptitle(title)
    plt.savefig(output_path)
    plt.close()


def normalise_alias(value: str) -> str:
    alias = value.replace("–", "-").replace("—", "-").replace("−", "-").strip()
    alias = re.sub(r"\s+", " ", alias)
    return alias.upper()


def parse_candidate_gene_names(value: str) -> list[str]:
    if value is None or str(value).strip() == "":
        return []
    return [token.strip() for token in re.split(r"[,/; ]+", str(value)) if token.strip()]


def parse_parenthetical_aliases(value: str) -> list[str]:
    aliases = []
    for alias_group in re.findall(r"\(([^()]*)\)", str(value)):
        aliases.extend([token.strip() for token in re.split(r"[,/]", alias_group) if token.strip()])
    return aliases


def parse_name_tokens(value: str) -> list[str]:
    tokens = []
    name = str(value).strip()
    if name and len(name) <= 40:
        tokens.append(name)
    tokens.extend(re.findall(r"\b[A-Z][A-Za-z0-9\-\?]{1,20}\b", name))
    return unique_preserve_order(tokens)


def direct_gene_match(candidates: list[str], available_genes: set[str]) -> tuple[str | None, list[str]]:
    cleaned_candidates = []

    for candidate in candidates:
        cleaned_candidate = normalise_alias(candidate)
        if cleaned_candidate == "":
            continue
        cleaned_candidates.append(cleaned_candidate)
        if cleaned_candidate in available_genes:
            return cleaned_candidate, unique_preserve_order(cleaned_candidates)

    return None, unique_preserve_order(cleaned_candidates)


def load_exhaustion_markers(markers_csv: Path, available_genes: set[str]) -> list[str]:
    marker_aliases = {"1D2": "ID2", "1D3": "ID3", "TIGID": "TIGIT"}
    corrected_markers = []

    with open(markers_csv, newline="", encoding="utf-8-sig") as handle:
        for raw_line in handle:
            marker = raw_line.strip().strip(";")
            if marker == "":
                continue
            marker = marker_aliases.get(marker.upper(), marker.upper())
            if marker in available_genes:
                corrected_markers.append(marker)

    corrected_markers = unique_preserve_order(corrected_markers)
    print(f"Loaded {len(corrected_markers)} corrected exhaustion markers present in the atlas")
    return corrected_markers


def load_inhibitory_receptors(receptors_csv: Path, available_genes: set[str]) -> pd.DataFrame:
    override_map = {
        "PD-1": "PDCD1",
        "CTLA-4": "CTLA4",
        "TIM-3": "HAVCR2",
        "TIM-1": "HAVCR1",
        "41BB": "TNFRSF9",
        "GITR": "TNFRSF18",
        "CD22": "SIGLEC2",
        "CD33": "SIGLEC3",
        "SIGLEC-5": "SIGLEC5",
        "SIGLEC-6": "SIGLEC6",
        "SIGLEC-7": "SIGLEC7",
        "SIGLEC-8": "SIGLEC8",
        "SIGLEC-9": "SIGLEC9",
        "SIGLEC-10": "SIGLEC10",
        "SIGLEC-11": "SIGLEC11",
        "SIGLEC-12": "SIGLEC12",
        "NKP44": "NCR2",
        "PILR A": "PILRA",
        "SIRP": "SIRPA",
        "SIRP?": "SIRPA",
        "CD21": "CR2",
        "OXIDIZED LOW-DENSITY LIPOPROTEIN RECEPTOR 1": "OLR1",
        "FC RECEPTOR-LIKE PROTEIN 1": "FCRL1",
        "FC RECEPTOR-LIKE PROTEIN 2": "FCRL2",
        "FC RECEPTOR-LIKE PROTEIN 3": "FCRL3",
        "FC RECEPTOR-LIKE PROTEIN 4": "FCRL4",
        "FC RECEPTOR-LIKE PROTEIN 5": "FCRL5",
        "FC RECEPTOR-LIKE PROTEIN 6": "FCRL6",
        "CD32A": "FCGR2A",
        "CD32B": "FCGR2B",
        "CD32C": "FCGR2C",
        "CD300A": "CD300A",
        "CD300F": "CD300LF",
        "CD28H": "TMIGD2",
        "MEGAKARYOCYTE AND PLATELET INHIBITORY RECEPTOR G6B": "MPIG6B",
    }

    mapping_rows = []

    with open(receptors_csv, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=";")

        for row in reader:
            receptor_name = str(row.get("Name", "")).strip()
            receptor_family = str(row.get("NTR family", "")).strip()
            receptor_uniprot = str(row.get("Uniprot", "")).strip()
            receptor_gene_name = str(row.get("Gene_name", "")).strip()

            is_control = receptor_uniprot.upper() == "CONTROL" or "TRUNCATED" in receptor_name.upper()

            gene_name_candidates = parse_candidate_gene_names(receptor_gene_name)
            parenthetical_candidates = parse_parenthetical_aliases(receptor_name)
            name_token_candidates = parse_name_tokens(receptor_name)

            mapped_gene = None
            mapping_source = "unresolved"

            mapped_gene, gene_name_candidates = direct_gene_match(gene_name_candidates, available_genes)
            if mapped_gene is not None:
                mapping_source = "gene_name"

            if mapped_gene is None:
                mapped_gene, parenthetical_candidates = direct_gene_match(parenthetical_candidates, available_genes)
                if mapped_gene is not None:
                    mapping_source = "parenthetical_alias"

            if mapped_gene is None:
                mapped_gene, name_token_candidates = direct_gene_match(name_token_candidates, available_genes)
                if mapped_gene is not None:
                    mapping_source = "name_token"

            override_candidates = []
            if mapped_gene is None:
                raw_override_candidates = [receptor_name] + gene_name_candidates + parenthetical_candidates + name_token_candidates
                raw_override_candidates.extend(re.findall(r"CD32[a-z]?", receptor_name, flags=re.IGNORECASE))

                for override_candidate in raw_override_candidates:
                    normalised_candidate = normalise_alias(str(override_candidate))
                    override_candidates.append(normalised_candidate)
                    if normalised_candidate in override_map and override_map[normalised_candidate] in available_genes:
                        mapped_gene = override_map[normalised_candidate]
                        mapping_source = "override"
                        break

            used_for_plots = mapped_gene is not None and not is_control

            mapping_rows.append(
                {
                    "receptor_family": receptor_family,
                    "receptor_name": receptor_name,
                    "uniprot": receptor_uniprot,
                    "gene_name_input": receptor_gene_name,
                    "mapped_gene": mapped_gene,
                    "mapping_source": mapping_source,
                    "gene_name_candidates": ";".join(gene_name_candidates),
                    "parenthetical_candidates": ";".join(parenthetical_candidates),
                    "name_token_candidates": ";".join(name_token_candidates),
                    "override_candidates": ";".join(unique_preserve_order(override_candidates)),
                    "is_control": is_control,
                    "used_for_plots": used_for_plots,
                }
            )

    mapping_df = pd.DataFrame(mapping_rows)
    mapping_df.to_csv(data_dir / "receptor_mapping_table.csv", index=False)
    print(
        "Loaded receptor table:",
        mapping_df.shape[0],
        "rows,",
        int(mapping_df["used_for_plots"].sum()),
        "receptors available for plots",
    )
    return mapping_df


def build_pseudobulk_table(
    expression_df: pd.DataFrame,
    obs_df: pd.DataFrame,
    group_cols: list[str],
) -> pd.DataFrame:
    work_df = pd.concat([obs_df[group_cols + ["sex_clean", "age_clean"]].copy(), expression_df.copy()], axis=1)

    grouped_expression = work_df.groupby(group_cols, observed=False)[expression_df.columns.tolist()].mean()
    grouped_metadata = work_df.groupby(group_cols, observed=False).agg(
        n_cells=(group_cols[0], "size"),
        sex_clean=("sex_clean", first_non_null),
        age_clean=("age_clean", first_non_null),
    )

    pseudobulk_df = grouped_metadata.join(grouped_expression).reset_index()
    return pseudobulk_df


def apply_fdr_by_group(results_df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    if results_df.empty:
        results_df["q_value"] = pd.Series(dtype=float)
        return results_df

    results_df = results_df.copy()
    results_df["q_value"] = np.nan

    for _, idx in results_df.groupby(group_cols, observed=False).groups.items():
        index_values = list(idx)
        valid_mask = results_df.loc[index_values, "p_value"].notna()
        valid_index = results_df.loc[index_values, :].index[valid_mask]

        if len(valid_index) == 0:
            continue

        _, q_values, _, _ = multipletests(results_df.loc[valid_index, "p_value"], method="fdr_bh")
        results_df.loc[valid_index, "q_value"] = q_values

    return results_df


def run_binary_context_models(
    expression_df: pd.DataFrame,
    obs_df: pd.DataFrame,
    subset_type: str,
    subset_label: str,
    context_variable: str,
) -> tuple[pd.DataFrame, list[dict]]:
    group_stats = (
        obs_df.groupby(context_variable, observed=False)
        .agg(
            n_cells=("individual_id", "size"),
            n_donors=("individual_id", "nunique"),
            n_studies=("study_name", "nunique"),
        )
        .reset_index()
    )

    pseudobulk_df = build_pseudobulk_table(
        expression_df=expression_df,
        obs_df=obs_df,
        group_cols=["individual_id", "study_name", "disease_type", "meta_tissue_type"],
    )

    model_results = []
    skipped_rows = []

    for _, group_row in group_stats.iterrows():
        group_name = group_row[context_variable]

        if pd.isna(group_name):
            continue

        if group_row["n_donors"] < 10 or group_row["n_cells"] < 200 or group_row["n_studies"] < 2:
            skipped_rows.append(
                {
                    "analysis_type": context_variable,
                    "subset_type": subset_type,
                    "subset_label": subset_label,
                    "group_name": group_name,
                    "reason": "underpowered_context_group",
                    "n_cells": int(group_row["n_cells"]),
                    "n_donors": int(group_row["n_donors"]),
                    "n_studies": int(group_row["n_studies"]),
                }
            )
            continue

        fit_df = pseudobulk_df.copy()
        fit_df["is_group"] = (fit_df[context_variable].astype(str) == str(group_name)).astype(int)

        if fit_df["is_group"].nunique() < 2:
            skipped_rows.append(
                {
                    "analysis_type": context_variable,
                    "subset_type": subset_type,
                    "subset_label": subset_label,
                    "group_name": group_name,
                    "reason": "no_binary_contrast_available",
                    "n_cells": int(group_row["n_cells"]),
                    "n_donors": int(group_row["n_donors"]),
                    "n_studies": int(group_row["n_studies"]),
                }
            )
            continue

        for gene in expression_df.columns:
            gene_fit_df = fit_df[[gene, "is_group", "study_name"]].rename(columns={gene: "expr"}).dropna()

            if gene_fit_df["is_group"].nunique() < 2:
                continue

            try:
                model = smf.ols("expr ~ is_group + C(study_name)", data=gene_fit_df).fit()
                model_results.append(
                    {
                        "analysis_type": context_variable,
                        "subset_type": subset_type,
                        "subset_label": subset_label,
                        "group_name": group_name,
                        "gene": gene,
                        "coefficient": model.params.get("is_group", np.nan),
                        "p_value": model.pvalues.get("is_group", np.nan),
                        "n_rows": int(gene_fit_df.shape[0]),
                        "n_group_rows": int(gene_fit_df["is_group"].sum()),
                        "n_group_donors": int(group_row["n_donors"]),
                        "n_group_cells": int(group_row["n_cells"]),
                        "n_group_studies": int(group_row["n_studies"]),
                    }
                )
            except Exception:
                model_results.append(
                    {
                        "analysis_type": context_variable,
                        "subset_type": subset_type,
                        "subset_label": subset_label,
                        "group_name": group_name,
                        "gene": gene,
                        "coefficient": np.nan,
                        "p_value": np.nan,
                        "n_rows": int(gene_fit_df.shape[0]),
                        "n_group_rows": int(gene_fit_df["is_group"].sum()),
                        "n_group_donors": int(group_row["n_donors"]),
                        "n_group_cells": int(group_row["n_cells"]),
                        "n_group_studies": int(group_row["n_studies"]),
                    }
                )

    model_results_df = pd.DataFrame(model_results)
    if not model_results_df.empty:
        model_results_df = apply_fdr_by_group(
            model_results_df,
            group_cols=["analysis_type", "subset_type", "subset_label", "group_name"],
        )

    return model_results_df, skipped_rows


def run_age_models(
    expression_df: pd.DataFrame,
    obs_df: pd.DataFrame,
    subset_type: str,
    subset_label: str,
) -> tuple[pd.DataFrame, list[dict]]:
    pseudobulk_df = build_pseudobulk_table(
        expression_df=expression_df,
        obs_df=obs_df,
        group_cols=["individual_id", "study_name", "disease_type"],
    )

    model_results = []
    skipped_rows = []

    for disease_type in major_disease_types_for_covariates:
        fit_df = pseudobulk_df.loc[pseudobulk_df["disease_type"] == disease_type].copy()
        fit_df = fit_df.loc[fit_df["age_clean"].notna()].copy()

        if fit_df.shape[0] < 20:
            skipped_rows.append(
                {
                    "analysis_type": "age",
                    "subset_type": subset_type,
                    "subset_label": subset_label,
                    "group_name": disease_type,
                    "reason": "fewer_than_20_age_annotated_donor_pseudobulks",
                    "n_cells": int(obs_df.loc[obs_df["disease_type"] == disease_type].shape[0]),
                    "n_donors": int(fit_df["individual_id"].nunique()),
                    "n_studies": int(fit_df["study_name"].nunique()),
                }
            )
            continue

        for gene in expression_df.columns:
            gene_fit_df = fit_df[[gene, "age_clean", "study_name"]].rename(columns={gene: "expr"}).dropna()
            if gene_fit_df.shape[0] < 20:
                continue

            try:
                model = smf.ols("expr ~ age_clean + C(study_name)", data=gene_fit_df).fit()
                model_results.append(
                    {
                        "analysis_type": "age",
                        "subset_type": subset_type,
                        "subset_label": subset_label,
                        "disease_type": disease_type,
                        "gene": gene,
                        "coefficient": model.params.get("age_clean", np.nan),
                        "p_value": model.pvalues.get("age_clean", np.nan),
                        "n_rows": int(gene_fit_df.shape[0]),
                    }
                )
            except Exception:
                model_results.append(
                    {
                        "analysis_type": "age",
                        "subset_type": subset_type,
                        "subset_label": subset_label,
                        "disease_type": disease_type,
                        "gene": gene,
                        "coefficient": np.nan,
                        "p_value": np.nan,
                        "n_rows": int(gene_fit_df.shape[0]),
                    }
                )

    model_results_df = pd.DataFrame(model_results)
    if not model_results_df.empty:
        model_results_df = apply_fdr_by_group(
            model_results_df,
            group_cols=["analysis_type", "subset_type", "subset_label", "disease_type"],
        )

    return model_results_df, skipped_rows


def run_sex_models(
    expression_df: pd.DataFrame,
    obs_df: pd.DataFrame,
    subset_type: str,
    subset_label: str,
) -> tuple[pd.DataFrame, list[dict]]:
    pseudobulk_df = build_pseudobulk_table(
        expression_df=expression_df,
        obs_df=obs_df,
        group_cols=["individual_id", "study_name", "disease_type"],
    )

    model_results = []
    skipped_rows = []

    for disease_type in major_disease_types_for_covariates:
        fit_df = pseudobulk_df.loc[pseudobulk_df["disease_type"] == disease_type].copy()
        fit_df = fit_df.loc[fit_df["sex_clean"].isin(["Male", "Female"])].copy()

        male_n = int((fit_df["sex_clean"] == "Male").sum())
        female_n = int((fit_df["sex_clean"] == "Female").sum())

        if male_n < 10 or female_n < 10:
            skipped_rows.append(
                {
                    "analysis_type": "sex",
                    "subset_type": subset_type,
                    "subset_label": subset_label,
                    "group_name": disease_type,
                    "reason": "fewer_than_10_male_or_female_donor_pseudobulks",
                    "n_cells": int(obs_df.loc[obs_df["disease_type"] == disease_type].shape[0]),
                    "n_donors": int(fit_df["individual_id"].nunique()),
                    "n_studies": int(fit_df["study_name"].nunique()),
                }
            )
            continue

        fit_df["sex_clean"] = pd.Categorical(fit_df["sex_clean"], categories=["Male", "Female"])

        for gene in expression_df.columns:
            gene_fit_df = fit_df[[gene, "sex_clean", "study_name"]].rename(columns={gene: "expr"}).dropna()
            if gene_fit_df["sex_clean"].nunique() < 2:
                continue

            try:
                model = smf.ols("expr ~ C(sex_clean) + C(study_name)", data=gene_fit_df).fit()
                model_results.append(
                    {
                        "analysis_type": "sex",
                        "subset_type": subset_type,
                        "subset_label": subset_label,
                        "disease_type": disease_type,
                        "gene": gene,
                        "coefficient": model.params.get("C(sex_clean)[T.Female]", np.nan),
                        "p_value": model.pvalues.get("C(sex_clean)[T.Female]", np.nan),
                        "n_rows": int(gene_fit_df.shape[0]),
                        "n_male": male_n,
                        "n_female": female_n,
                    }
                )
            except Exception:
                model_results.append(
                    {
                        "analysis_type": "sex",
                        "subset_type": subset_type,
                        "subset_label": subset_label,
                        "disease_type": disease_type,
                        "gene": gene,
                        "coefficient": np.nan,
                        "p_value": np.nan,
                        "n_rows": int(gene_fit_df.shape[0]),
                        "n_male": male_n,
                        "n_female": female_n,
                    }
                )

    model_results_df = pd.DataFrame(model_results)
    if not model_results_df.empty:
        model_results_df = apply_fdr_by_group(
            model_results_df,
            group_cols=["analysis_type", "subset_type", "subset_label", "disease_type"],
        )

    return model_results_df, skipped_rows


# 0. load data
adata = sc.read_h5ad(str(input_h5ad))
print(f"Loaded: {input_h5ad} with shape {adata.shape}")

## subsampling

max_cells = 100_000
if max_cells and adata.n_obs > max_cells:
    idx = np.random.default_rng(random_seed).choice(adata.n_obs, size=max_cells, replace=False)
    adata = adata[idx, :].copy()
    print(f"After subsampling: {adata.shape}")


## preserve original counts
adata.layers["counts"] = adata.X.copy()

# 1. QC

## 1.1 ambient RNA QC
### check if necessary

low_count_treshold = 1_000
frac_treshold = 0.05

low_count_cells = np.asarray(adata.X.sum(axis=1)).ravel() < low_count_treshold
if low_count_cells.mean() > frac_treshold:
    print("potential RNA contamination, clean up the data")
    exit()


## 1.2 count QC
### QC metrics

#### mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
#### ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
#### hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],  # what percentage of counts is mito, ribosomal, hemoglobin
    percent_top=[20, 50, 100],  # what percenteage of counts is top 20, 50, 100 genes
    inplace=True,
)

### Get outliers with MAD


def outlier_MAD(adata, metric: str, nmads: int) -> pd.Series:
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * stats.median_abs_deviation(M)) | (
        np.median(M) + nmads * stats.median_abs_deviation(M) < M
    )
    return outlier


adata.obs["outlier"] = (
    outlier_MAD(adata, "log1p_total_counts", 5)
    | outlier_MAD(adata, "log1p_n_genes_by_counts", 5)
    | outlier_MAD(adata, "pct_counts_in_top_20_genes", 5)
)
adata.obs.outlier.value_counts()
adata.obs.outlier.value_counts(normalize=True)

adata.obs["mt_outlier"] = outlier_MAD(adata, "pct_counts_mt", 3) | (adata.obs["pct_counts_mt"] > 8)
adata.obs.mt_outlier.value_counts()
adata.obs.mt_outlier.value_counts(normalize=True)

### plot before QC
sns.displot(adata.obs["total_counts"], bins=100, kde=True)
plt.savefig(figures_dir / "qc_before_total_counts_hist.png")
plt.close()
sns.displot(adata.obs["pct_counts_mt"], bins=100, kde=True)
plt.savefig(figures_dir / "qc_before_pct_counts_mt_hist.png")
plt.close()

sc.pl.violin(adata, "pct_counts_mt", show=False)
plt.savefig(figures_dir / "qc_before_pct_counts_mt_violin.png", bbox_inches="tight")
plt.close()
qc_scatter_genes = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
plt.savefig(figures_dir / "qc_before_scatter_genes_total_mt.png")
plt.close()


### Filter
print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

### plot after count QC
sns.displot(adata.obs["total_counts"], bins=100, kde=True).savefig(figures_dir / "qc_after_total_counts_hist.png")
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
plt.savefig(figures_dir / "qc_after_scatter_genes_total_mt.png")
plt.close()


## 1.3 Doublet QC
n_obs_before_doublet_qc = adata.n_obs  # for statistics
### prepare input
counts_csr = adata.layers["counts"].tocsr() if not sp.issparse(adata.layers["counts"]) else adata.layers["counts"].tocsr()
scrublet = scr.Scrublet(counts_matrix=counts_csr, expected_doublet_rate=0.05)

### run scrublet
doublet_scores, predicted_doublets = scrublet.scrub_doublets()

### save scores, calls, histogram
adata.obs["doublet_score"] = doublet_scores
adata.obs["doublet_call"] = pd.Categorical(predicted_doublets)
adata.uns["scrublet_sim_scores"] = scrublet.doublet_scores_sim_

scrublet.plot_histogram()
plt.savefig(figures_dir / "qc_scrublet_histogram.png", dpi=300)
plt.close()

### Filter if necessary
if adata.obs["doublet_call"].nunique() > 1:
    print("There are predicted doublets, applying filter")
    adata = adata[~adata.obs["doublet_call"].astype(bool).to_numpy(), :].copy()
else:
    print("No doublets predicted")
print(f"After doublet step: {adata.shape} (kept {adata.n_obs}/{n_obs_before_doublet_qc} cells)")


# 2. Normalisation

# ## retrieve count values - not necessary now
# adata.X = adata.layers["counts"].copy()

sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
adata.layers["normalised"] = adata.X.copy()

sc.pp.log1p(adata)
adata.X = adata.X.astype(np.float32)  # for better storage
adata.layers["log1p"] = adata.X.copy()

log1p_library_size = np.asarray(adata.layers["log1p"].sum(axis=1)).ravel()
log1p_library_size


sns.histplot(adata.obs["total_counts"], bins=100, kde=False)
plt.savefig(figures_dir / "normalisation_total_counts_hist.png")
plt.close()
sns.histplot(log1p_library_size, bins=100, kde=False)
plt.savefig(figures_dir / "normalisation_log1p_hist.png")
plt.close()


# ## save progress
# adata.write_h5ad(data_dir / "qc+normalised.h5ad")

# ## retrieve progress
# adata = sc.read_h5ad(data_dir / "qc+normalised.h5ad")
# print(f"Loaded data with shape {adata.shape}")


## preserve full-gene data
adata_full = adata.copy()
available_genes = set(adata_full.var_names.tolist())

exhaustion_markers = load_exhaustion_markers(markers_csv, available_genes)
receptor_mapping_df = load_inhibitory_receptors(receptors_csv, available_genes)
receptor_genes = unique_preserve_order(
    receptor_mapping_df.loc[receptor_mapping_df["used_for_plots"], "mapped_gene"].dropna().tolist()
)
receptors_not_in_hvg = []
print(f"Mapped receptor genes available for plots: {len(receptor_genes)}")


# 3. Feature selection
## 3.1 Cell selection - none here

## 3.2. HVG selection
### using seurat as seurat_v3 with original counts raises an environment error
adata_hvg = adata_full.copy()
sc.pp.highly_variable_genes(adata_hvg, layer="log1p", flavor="seurat", n_top_genes=4000, inplace=True)
print("HVG selected:", int(adata_hvg.var["highly_variable"].sum()))

# sc.pl.highly_variable_genes(adata_hvg, show=False)
sns.scatterplot(data=adata_hvg.var, x="means", y="dispersions", hue="highly_variable", s=5)
plt.savefig(figures_dir / "hvg_before.png")
plt.close()

### Filter
adata_hvg = adata_hvg[:, adata_hvg.var.highly_variable].copy()

# sc.pl.highly_variable_genes(adata_hvg, show=False)
sns.scatterplot(data=adata_hvg.var, x="means", y="dispersions", hue="highly_variable", s=5)
plt.savefig(figures_dir / "hvg_after.png")
plt.close()

# ## save progress
# adata_hvg.write_h5ad(data_dir / "qc+normalised+hvg.h5ad")
# ## retrieve progress
# adata_hvg = sc.read_h5ad(data_dir / "qc+normalised+hvg.h5ad")
# print(f"Loaded data with shape {adata_hvg.shape}")

# plot top 20 genes for sanity check
sc.pl.highest_expr_genes(adata_hvg, n_top=20, show=False)
plt.savefig(figures_dir / "top20_expr_genes.png")
plt.close()


# 4. Dimension reduction (=Embedding)

## PCA (linear reduction)
n_pcs = 30
sc.pp.scale(adata_hvg, max_value=10)  # scale to variance
sc.tl.pca(adata_hvg, n_comps=n_pcs, svd_solver="arpack", random_state=random_seed)

## examine PCA - plot variance ratio
sc.pl.pca_variance_ratio(adata_hvg, n_pcs=n_pcs, log=True, show=False)
plt.savefig(figures_dir / "pca_variance_ratio.png")
plt.close()

## examine PCA - plot scatter with different colours
pca_color_cols = ["cell_subtype_3", "disease_type", "meta_tissue_type", "study_name", "total_counts"]
sc.pl.pca_scatter(adata_hvg, color=pca_color_cols, ncols=2, show=False)
plt.savefig(figures_dir / "pca_scatter.png", dpi=100)
plt.close()


# 5. Clustering
## create k neighbour graph
n_neighbors = 15
sc.pp.neighbors(adata_hvg, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_seed)
sc.tl.umap(adata_hvg, min_dist=0.3, spread=1.0, random_state=random_seed)

## detect clusters at one resolution (test pipeline default)
leiden_resolution = 0.5
leiden_key = "leiden_res_0.5"
sc.tl.leiden(adata_hvg, resolution=leiden_resolution, key_added=leiden_key)

## save umap with different colourings
umap_color_cols = [leiden_key, "cell_subtype_3", "disease_type", "meta_tissue_type"]
sc.pl.umap(adata_hvg, color=umap_color_cols, ncols=2, wspace=0.35, show=False)
plt.savefig(figures_dir / "umap_clustering_overview.png")
plt.close()

## transfer embedding back to full-gene data
copy_embedding_results(
    adata_source=adata_hvg,
    adata_target=adata_full,
    obs_cols=[leiden_key],
    uns_keys=["neighbors", "pca", "umap", leiden_key, f"{leiden_key}_colors"],
)

adata_full.obs["is_exhausted_like"] = pd.Categorical(
    np.where(adata_full.obs["cell_subtype_3"].isin(exhausted_reference_states), "Exhausted-like", "Other")
)

receptors_not_in_hvg = [gene for gene in receptor_genes if gene not in adata_hvg.var_names]
pd.DataFrame({"mapped_receptor_gene": receptor_genes, "in_hvg_object": [gene in adata_hvg.var_names for gene in receptor_genes]}).to_csv(
    data_dir / "receptor_hvg_overlap.csv", index=False
)
print(f"Mapped receptor genes not retained in the HVG object: {len(receptors_not_in_hvg)}")


# 6. Preliminary plots
## reference-state overview
sc.pl.umap(
    adata_full,
    color=["cell_subtype_3", "disease_type", "meta_tissue_type", "is_exhausted_like"],
    ncols=2,
    wspace=0.35,
    show=False,
)
plt.savefig(figures_dir / "umap_reference_overview.png")
plt.close()

reference_states_present = [state for state in reference_state_order if state in adata_full.obs["cell_subtype_3"].astype(str).unique()]
adata_reference = adata_full[adata_full.obs["cell_subtype_3"].isin(reference_states_present), :].copy()
adata_reference.obs["cell_subtype_3"] = pd.Categorical(
    adata_reference.obs["cell_subtype_3"].astype(str),
    categories=reference_states_present,
    ordered=True,
)

save_dotplot(
    adata_plot=adata_reference,
    genes=exhaustion_markers,
    groupby="cell_subtype_3",
    output_path=figures_dir / "dotplot_exhaustion_markers_by_reference_state.png",
    categories_order=reference_states_present,
    title="Corrected exhaustion markers across reference states",
)

if len(receptor_genes) <= 70:
    save_matrixplot(
        adata_plot=adata_reference,
        genes=receptor_genes,
        groupby="cell_subtype_3",
        output_path=figures_dir / "matrixplot_inhibitory_receptors_by_reference_state.png",
        categories_order=reference_states_present,
        title="Mapped inhibitory receptors across reference states",
    )
else:
    receptor_gene_chunks = split_into_chunks(receptor_genes, 30)
    for chunk_idx, receptor_gene_chunk in enumerate(receptor_gene_chunks):
        file_name = (
            "matrixplot_inhibitory_receptors_by_reference_state.png"
            if chunk_idx == 0
            else f"matrixplot_inhibitory_receptors_by_reference_state_page_{chunk_idx + 1}.png"
        )
        save_matrixplot(
            adata_plot=adata_reference,
            genes=receptor_gene_chunk,
            groupby="cell_subtype_3",
            output_path=figures_dir / file_name,
            categories_order=reference_states_present,
            title=f"Mapped inhibitory receptors across reference states (page {chunk_idx + 1})",
        )


# 7. Exhausted-only clustering and modules
exhausted_mask = adata_full.obs["cell_subtype_3"].isin(exhausted_reference_states)
adata_exh = adata_full[exhausted_mask, :].copy()
print(f"Exhausted-like object shape: {adata_exh.shape}")

adata_exh.obs["sex_clean"] = clean_sex_column(adata_exh.obs["donor_sex"])
adata_exh.obs["age_clean"] = clean_age_column(adata_exh.obs["donor_age"])

adata_exh_hvg = adata_exh.copy()
sc.pp.highly_variable_genes(adata_exh_hvg, layer="log1p", flavor="seurat", n_top_genes=4000, inplace=True)
adata_exh_hvg = adata_exh_hvg[:, adata_exh_hvg.var["highly_variable"]].copy()
sc.pp.scale(adata_exh_hvg, max_value=10)
sc.tl.pca(adata_exh_hvg, n_comps=20, svd_solver="arpack", random_state=random_seed)
sc.pp.neighbors(adata_exh_hvg, n_neighbors=15, n_pcs=20, random_state=random_seed)
sc.tl.umap(adata_exh_hvg, min_dist=0.3, spread=1.0, random_state=random_seed)

exh_leiden_key = "leiden_res_0.4_exhausted"
sc.tl.leiden(adata_exh_hvg, resolution=0.4, key_added=exh_leiden_key)

copy_embedding_results(
    adata_source=adata_exh_hvg,
    adata_target=adata_exh,
    obs_cols=[exh_leiden_key],
    uns_keys=["neighbors", "pca", "umap", exh_leiden_key, f"{exh_leiden_key}_colors"],
)

corrected_marker_set = set(exhaustion_markers)
module_gene_sets = {
    "progenitor_module": [gene for gene in ["CCR7", "SELL", "LEF1", "TCF7", "CD28", "SLAMF6"] if gene in adata_exh.var_names],
    "cytotoxic_module": [gene for gene in ["IFNG", "PRF1", "GZMA", "GZMB", "GNLY", "TBX21", "CX3CR1", "KLRG1", "ZEB2", "PRDM1"] if gene in adata_exh.var_names],
    "terminal_exhaustion_module": [
        gene
        for gene in ["TOX", "TOX2", "ENTPD1", "BATF", "LAG3", "TIGIT", "PDCD1", "CTLA4", "TRIB1", "HAVCR2"]
        if gene in adata_exh.var_names
    ],
}

for module_name, module_genes in module_gene_sets.items():
    score_col = f"{module_name}_score"
    if len(module_genes) == 0:
        adata_exh.obs[score_col] = np.nan
        continue
    sc.tl.score_genes(
        adata_exh,
        gene_list=module_genes,
        score_name=score_col,
        random_state=random_seed,
        use_raw=False,
    )

score_cols = [f"{module_name}_score" for module_name in module_gene_sets]
cluster_order = sorted(adata_exh.obs[exh_leiden_key].astype(str).unique(), key=lambda value: int(value))
cluster_name_map = {}
cluster_base_name_counts = defaultdict(int)
cluster_annotation_rows = []

for cluster_id in cluster_order:
    cluster_obs = adata_exh.obs.loc[adata_exh.obs[exh_leiden_key].astype(str) == cluster_id].copy()
    dominant_reference_state = cluster_obs["cell_subtype_3"].astype(str).value_counts().idxmax()
    cluster_score_means = cluster_obs[score_cols].mean()
    top_module = cluster_score_means.idxmax()

    if dominant_reference_state == "Tpex":
        base_cluster_name = "Tpex-like progenitor"
    elif dominant_reference_state == "GZMK+ Tex":
        base_cluster_name = "GZMK-Tex cytotoxic"
    elif dominant_reference_state == "XBP1+ Tex":
        base_cluster_name = "XBP1-Tex stress-adapted"
    else:
        base_cluster_name = "ITGAE-Tex terminal"

    cluster_base_name_counts[base_cluster_name] += 1
    if cluster_base_name_counts[base_cluster_name] == 1:
        cluster_name = base_cluster_name
    else:
        cluster_name = f"{base_cluster_name}_{cluster_base_name_counts[base_cluster_name]}"

    cluster_name_map[cluster_id] = cluster_name
    cluster_annotation_rows.append(
        {
            "cluster_id": cluster_id,
            "cluster_name": cluster_name,
            "dominant_reference_state": dominant_reference_state,
            "top_module": top_module,
            "n_cells": int(cluster_obs.shape[0]),
            **{score_col: cluster_score_means.get(score_col, np.nan) for score_col in score_cols},
        }
    )

adata_exh.obs["exhausted_cluster_name"] = adata_exh.obs[exh_leiden_key].astype(str).map(cluster_name_map)
named_cluster_order = [cluster_name_map[cluster_id] for cluster_id in cluster_order]
adata_exh.obs["exhausted_cluster_name"] = pd.Categorical(
    adata_exh.obs["exhausted_cluster_name"],
    categories=named_cluster_order,
    ordered=True,
)

cluster_annotation_df = pd.DataFrame(cluster_annotation_rows)
cluster_annotation_df.to_csv(data_dir / "exhausted_cluster_annotations.csv", index=False)

sc.pl.umap(
    adata_exh,
    color=["cell_subtype_3", exh_leiden_key, "disease_type", "meta_tissue_type", "exhausted_cluster_name"],
    ncols=2,
    wspace=0.35,
    show=False,
)
plt.savefig(figures_dir / "umap_exhausted_only.png")
plt.close()


# 8. Receptor expression summaries in exhausted-like cells
receptor_expression_exh = adata_exh[:, receptor_genes].to_df() if receptor_genes else pd.DataFrame(index=adata_exh.obs_names)
receptor_detection = receptor_expression_exh.gt(0).mean(axis=0).sort_values(ascending=False) if receptor_genes else pd.Series(dtype=float)
top_receptor_genes = receptor_detection.head(20).index.tolist()

disease_type_order = adata_exh.obs["disease_type"].astype(str).value_counts().index.tolist()
meta_tissue_type_order = adata_exh.obs["meta_tissue_type"].astype(str).value_counts().index.tolist()

adata_exh.obs["disease_type"] = pd.Categorical(adata_exh.obs["disease_type"].astype(str), categories=disease_type_order, ordered=True)
adata_exh.obs["meta_tissue_type"] = pd.Categorical(
    adata_exh.obs["meta_tissue_type"].astype(str),
    categories=meta_tissue_type_order,
    ordered=True,
)

save_matrixplot(
    adata_plot=adata_exh,
    genes=top_receptor_genes,
    groupby="disease_type",
    output_path=figures_dir / "matrixplot_top_inhibitory_receptors_exhausted_by_disease_type.png",
    categories_order=disease_type_order,
    title="Top inhibitory receptors in exhausted-like cells by disease type",
)

save_matrixplot(
    adata_plot=adata_exh,
    genes=top_receptor_genes,
    groupby="meta_tissue_type",
    output_path=figures_dir / "matrixplot_top_inhibitory_receptors_exhausted_by_meta_tissue_type.png",
    categories_order=meta_tissue_type_order,
    title="Top inhibitory receptors in exhausted-like cells by tissue context",
)


# 9. Differential expression and heatmaps
rank_genes_key = "rank_genes_groups_exhausted"
sc.tl.rank_genes_groups(adata_exh, groupby="exhausted_cluster_name", method="wilcoxon", use_raw=False, key_added=rank_genes_key)
rank_genes_df = sc.get.rank_genes_groups_df(adata_exh, group=None, key=rank_genes_key)
rank_genes_df.to_csv(data_dir / "exhausted_cluster_markers.csv", index=False)

top_marker_genes = []
for cluster_name in named_cluster_order:
    cluster_marker_df = rank_genes_df.loc[rank_genes_df["group"] == cluster_name].copy()
    if "logfoldchanges" in cluster_marker_df.columns:
        cluster_marker_df = cluster_marker_df.loc[cluster_marker_df["logfoldchanges"] > 0].copy()
    top_marker_genes.extend(cluster_marker_df["names"].head(5).tolist())

heatmap_genes = unique_preserve_order(top_marker_genes + exhaustion_markers + top_receptor_genes)
save_heatmap(
    adata_plot=adata_exh,
    genes=heatmap_genes,
    groupby="exhausted_cluster_name",
    output_path=figures_dir / "heatmap_exhausted_cluster_markers.png",
    title="Top exhausted-cluster markers, exhaustion markers, and inhibitory receptors",
)

coexpression_genes = receptor_detection.loc[receptor_detection >= 0.01].index.tolist()
if len(coexpression_genes) >= 2:
    receptor_binary_df = receptor_expression_exh[coexpression_genes].gt(0).astype(int)
    phi_corr_df = receptor_binary_df.corr(method="pearson")
    cluster_grid = sns.clustermap(phi_corr_df, cmap="vlag", center=0, linewidths=0.1, figsize=(10, 10))
    cluster_grid.fig.suptitle("Receptor co-expression in exhausted-like cells")
    cluster_grid.savefig(figures_dir / "heatmap_receptor_coexpression_exhausted.png")
    plt.close("all")
else:
    save_placeholder_figure(
        figures_dir / "heatmap_receptor_coexpression_exhausted.png",
        "Fewer than two inhibitory receptors passed the 1% detection threshold.",
    )


# 10. Module summaries and receptor-family summaries
cluster_module_summary_df = (
    adata_exh.obs.groupby("exhausted_cluster_name", observed=False)[score_cols]
    .mean()
    .reset_index()
    .merge(cluster_annotation_df, left_on="exhausted_cluster_name", right_on="cluster_name", how="left")
)
cluster_module_summary_df.to_csv(data_dir / "exhaustion_module_scores_by_cluster.csv", index=False)

family_to_genes = {}
for receptor_family, family_df in receptor_mapping_df.loc[receptor_mapping_df["used_for_plots"]].groupby("receptor_family", observed=False):
    family_genes = unique_preserve_order(family_df["mapped_gene"].dropna().tolist())
    if len(family_genes) > 0:
        family_to_genes[receptor_family] = family_genes

receptor_family_summary_rows = []
for context_variable in ["disease_type", "meta_tissue_type"]:
    if receptor_expression_exh.empty:
        continue

    family_score_df = pd.DataFrame(index=adata_exh.obs_names)
    for receptor_family, family_genes in family_to_genes.items():
        family_score_df[receptor_family] = receptor_expression_exh[family_genes].mean(axis=1)

    family_score_df[context_variable] = adata_exh.obs[context_variable].astype(str).values
    group_counts = family_score_df[context_variable].value_counts()
    grouped_family_scores = family_score_df.groupby(context_variable, observed=False).mean(numeric_only=True)

    for context_group, grouped_row in grouped_family_scores.iterrows():
        for receptor_family in family_to_genes:
            receptor_family_summary_rows.append(
                {
                    "context_variable": context_variable,
                    "context_group": context_group,
                    "receptor_family": receptor_family,
                    "mean_family_score": grouped_row.get(receptor_family, np.nan),
                    "n_cells": int(group_counts.get(context_group, 0)),
                }
            )

receptor_family_summary_df = pd.DataFrame(receptor_family_summary_rows)
receptor_family_summary_df.to_csv(data_dir / "receptor_family_scores_by_context.csv", index=False)


# 11. Donor-level context models
context_results = []
age_results = []
sex_results = []
skipped_rows = []

context_subsets = [("all_exhausted", "All exhausted-like", adata_exh.obs_names)]

for reference_state in exhausted_reference_states:
    state_index = adata_exh.obs_names[adata_exh.obs["cell_subtype_3"].astype(str) == reference_state]
    if len(state_index) > 0:
        context_subsets.append(("reference_state", reference_state, state_index))

for cluster_name in named_cluster_order:
    cluster_index = adata_exh.obs_names[adata_exh.obs["exhausted_cluster_name"].astype(str) == cluster_name]
    cluster_donors = adata_exh.obs.loc[cluster_index, "individual_id"].nunique()
    if cluster_donors >= 5:
        context_subsets.append(("cluster", cluster_name, cluster_index))

for subset_type, subset_label, subset_index in context_subsets:
    subset_obs = adata_exh.obs.loc[subset_index].copy()
    subset_expr = receptor_expression_exh.loc[subset_index].copy()

    disease_results_df, disease_skipped_rows = run_binary_context_models(
        expression_df=subset_expr,
        obs_df=subset_obs,
        subset_type=subset_type,
        subset_label=subset_label,
        context_variable="disease_type",
    )
    tissue_results_df, tissue_skipped_rows = run_binary_context_models(
        expression_df=subset_expr,
        obs_df=subset_obs,
        subset_type=subset_type,
        subset_label=subset_label,
        context_variable="meta_tissue_type",
    )

    context_results.append(disease_results_df)
    context_results.append(tissue_results_df)
    skipped_rows.extend(disease_skipped_rows)
    skipped_rows.extend(tissue_skipped_rows)

all_exhausted_age_df, age_skipped_rows = run_age_models(
    expression_df=receptor_expression_exh,
    obs_df=adata_exh.obs.copy(),
    subset_type="all_exhausted",
    subset_label="All exhausted-like",
)
all_exhausted_sex_df, sex_skipped_rows = run_sex_models(
    expression_df=receptor_expression_exh,
    obs_df=adata_exh.obs.copy(),
    subset_type="all_exhausted",
    subset_label="All exhausted-like",
)

age_results.append(all_exhausted_age_df)
sex_results.append(all_exhausted_sex_df)
skipped_rows.extend(age_skipped_rows)
skipped_rows.extend(sex_skipped_rows)

context_results_df = pd.concat(context_results, ignore_index=True) if len(context_results) > 0 else pd.DataFrame()
age_results_df = pd.concat(age_results, ignore_index=True) if len(age_results) > 0 else pd.DataFrame()
sex_results_df = pd.concat(sex_results, ignore_index=True) if len(sex_results) > 0 else pd.DataFrame()
skipped_results_df = pd.DataFrame(skipped_rows)

context_results_df.to_csv(data_dir / "receptor_context_associations.csv", index=False)
age_results_df.to_csv(data_dir / "receptor_age_associations.csv", index=False)
sex_results_df.to_csv(data_dir / "receptor_sex_associations.csv", index=False)
skipped_results_df.to_csv(data_dir / "skipped_context_comparisons.csv", index=False)


# 12. save processed data
adata_hvg.write_h5ad(data_dir / "adata_processed.h5ad")
adata_full.write_h5ad(data_dir / "adata_full_processed.h5ad")
adata_exh.write_h5ad(data_dir / "adata_exhausted_processed.h5ad")
