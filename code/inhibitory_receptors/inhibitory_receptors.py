from cd8atlas.code import pipeline_elements as pipeline
from cd8atlas.code.pipeline_elements import bioreport
import re



def load_inhibitory_receptor_genes(filename="inhibitory_receptor_list.csv") -> list[str]:
    receptor_file = bioreport.file_load(
        filename,
        data_dir=pipeline.root_dir / "code" / "inhibitory_receptors",
    )
    receptor_df = pipeline.pd.read_csv(receptor_file, sep=";", dtype=str).fillna("")
    receptor_genes: set[str] = set()
    manual_name_map = {
        "CD300a": "CD300A",
        "CD300f": "CD300F",
        "Fc receptor-like protein 1": "FCRL1",
        "Fc receptor-like protein 2": "FCRL2",
        "Fc receptor-like protein 3": "FCRL3",
        "Fc receptor-like protein 4": "FCRL4",
        "Fc receptor-like protein 5": "FCRL5",
        "Fc receptor-like protein 6": "FCRL6",
        "Fc?RIIa (CD32a)": "FCGR2A",
        "Fc?RIIb (CD32b)": "FCGR2B",
        "Fc?RIIc (CD32c)": "FCGR2C",
        "PILR A": "PILRA",
        "CD28H (Transmembrane and immunoglobulin domain-containing protein 2)": "TMIGD2",
        "SIRP?": "SIRPA",
        "Megakaryocyte and platelet inhibitory receptor G6b": "MPIG6B",
    }

    for _, row in receptor_df.iterrows():
        gene_name = row["Gene_name"].strip()
        if gene_name:
            receptor_genes.add(gene_name.replace("-", ""))
            continue

        name = row["Name"].strip()
        if not name:
            continue
        if name in manual_name_map:
            receptor_genes.add(manual_name_map[name])
            continue

        candidates: list[str] = []
        paren_matches = re.findall(r"\(([^)]+)\)", name)
        for match in paren_matches:
            candidates.extend(part.strip() for part in match.split(","))
        if not candidates:
            candidates = [part.strip() for part in name.split(",")]
        candidates.append(name)

        for candidate in candidates:
            normalized = re.sub(r"[^A-Za-z0-9]", "", candidate).upper()
            if re.fullmatch(r"[A-Z0-9]+", normalized) and any(char.isalpha() for char in normalized):
                receptor_genes.add(normalized)
                break

    return sorted(receptor_genes)


def label_inhibitory_receptor_score(
    adata,
    receptor_genes: list[str],
    score_col: str = "inhibitory_receptor_score",
):
    feature_symbols = pipeline.pd.Series(
        adata.var["gene_symbol"].astype(str).to_numpy(),
        index=adata.var_names,
    )
    matched_features = feature_symbols.index[feature_symbols.isin(receptor_genes)].tolist()
    pipeline.sc.tl.score_genes(
        adata,
        gene_list=matched_features,
        score_name=score_col,
        ctrl_size=max(1, min(50, len(matched_features))),
        random_state=pipeline.random_seed,
        use_raw=False,
        layer="log1p",
    )
    bioreport.log(f"Calculated inhibitory receptor score from {len(matched_features)} genes", heading=3)
    return adata


def inhibitory_receptor_heatmap(
    adata,
    receptor_genes: list[str],
    groupby_col: str = "cluster_annotation",
    figure_name: str = "heatmap_cluster_annotation_X_inhibitory_receptors",
):
    feature_symbols = pipeline.pd.Series(
        adata.var["gene_symbol"].astype(str).to_numpy(),
        index=adata.var_names,
    )
    matched_features = feature_symbols.index[feature_symbols.isin(receptor_genes)].tolist()
    adata.var["inhibitory_receptors_selected"] = adata.var_names.isin(matched_features)
    adata_plot = adata[:, matched_features].copy()
    adata_plot.var_names = pipeline.plot_var_names(adata_plot)
    heatmap_layer = "log1p" if "log1p" in adata.layers else None

    expr = adata_plot.layers[heatmap_layer] if heatmap_layer is not None else adata_plot.X
    if pipeline.sp.issparse(expr):
        expr = expr.toarray()
    expr_df = pipeline.pd.DataFrame(
        expr,
        index=adata_plot.obs[groupby_col].astype(str).to_numpy(),
        columns=adata_plot.var_names,
    )
    mean_expr = expr_df.groupby(level=0, sort=False).mean()
    mean_expr = mean_expr.loc[mean_expr.index != ""]
    heatmap_df = mean_expr.T
    heatmap_df = heatmap_df.sub(heatmap_df.mean(axis=1), axis=0)
    row_std = heatmap_df.std(axis=1).replace(0, 1)
    heatmap_df = heatmap_df.div(row_std, axis=0)

    fig_width = max(8, 0.8 * heatmap_df.shape[1] + 4)
    fig_height = max(14, 0.28 * heatmap_df.shape[0] + 2)
    pipeline.plt.figure(figsize=(fig_width, fig_height))
    pipeline.sns.heatmap(
        heatmap_df,
        cmap="vlag",
        center=0,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "Mean scaled log1p expression"},
    )
    pipeline.plt.xlabel(groupby_col)
    pipeline.plt.ylabel("Inhibitory receptor genes")
    pipeline.plt.xticks(rotation=45, ha="right")
    pipeline.plt.yticks(rotation=0)
    pipeline.plt.tight_layout()
    bioreport.figure(figure_name)
    return adata



pipeline.init("inhibitory_receptors")

adata = pipeline.load_data("qc+subsampled_100000.h5ad")

adata = pipeline.label_genes(adata)
receptor_genes = load_inhibitory_receptor_genes()
adata = pipeline.feature_selection_full(
    adata,
    include_markers="state_markers.csv",
    whitelist_genes=receptor_genes,
)
dataset_genes = set(map(str, adata.var_names))
dataset_genes.update(gene for gene in adata.var["gene_symbol"].astype(str) if gene)
present_receptors = sorted(gene for gene in receptor_genes if gene in dataset_genes)
missing_receptors = sorted(gene for gene in receptor_genes if gene not in dataset_genes)

bioreport.log(f"Missing inhibitory receptor genes: {len(missing_receptors)}, present genes: {len(present_receptors)}", heading=3)

adata = pipeline.umap_basic(
    adata,
    n_pcs=25,
    n_neighbors=12,
    min_dist=0.1,
    leiden_resolution=0.4,
    save=True,
)

adata = label_inhibitory_receptor_score(adata, receptor_genes)
adata = pipeline.label_states( adata, "state_markers.csv")  # from https://www.biorxiv.org/content/10.1101/2024.07.17.603780v1.full
adata = pipeline.load_cluster_state_summary_csv(
    adata,
    filename="cluster_state_summary_names.csv",
    cluster_col="leiden",
    annotation_col="cluster_annotation",
    source_col="manual_cluster_name",
    data_dir=pipeline.root_dir / "data",
)
pipeline.umap_plot(adata,
    umap_color_cols=["leiden", "state_label", "state_score_max", "state_score_margin", "inhibitory_receptor_score"],
    name="umap_states")

adata.var["state_markers_selected"] = adata.var_names.isin(pipeline.marker_genes_from_csv("state_markers.csv"))


pipeline.dotplot(adata, cluster_col="state_label", genes_bool_col="state_markers_selected", figure_name="dotplot_state_label_X_state_markers")

pipeline.dotplot(adata, cluster_col="leiden",
    genes_bool_col="state_markers_selected",
    gene_groups=pipeline.marker_gene_groups_from_csv("state_markers.csv"),
    swap_axes=False,
    figure_name=f"dotplot_leiden_X_state_markers")

adata = pipeline.dotplot_allmarkers(adata, n_markers_per_cluster=3, cluster_col="leiden")
pipeline.export_cluster_state_summary_csv(adata, cluster_col="leiden", state_col="state_label", filename="cluster_state_summary.csv")
adata = inhibitory_receptor_heatmap(adata, receptor_genes, groupby_col="cluster_annotation")



pipeline.export()
