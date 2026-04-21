from cd8atlas.code import pipeline_elements as pipeline
from cd8atlas.code.pipeline_elements import bioreport


cluster_col = "leiden"

pipeline.init("inhibitory_receptors")

adata = pipeline.load_data("qc+subsampled_100000.h5ad")

adata = pipeline.label_genes(adata)
adata = pipeline.feature_selection_full(adata, include_markers="state_markers.csv")

adata = pipeline.umap_basic(
    adata,
    n_pcs=25,
    n_neighbors=10,
    min_dist=0.2,
    leiden_resolution=0.6,
    save=True,
)

adata = pipeline.label_states(
    adata,
    "state_markers.csv",
)  # from https://www.biorxiv.org/content/10.1101/2024.07.17.603780v1.full
pipeline.umap_plot(
    adata,
    umap_color_cols=[cluster_col, "state_label", "state_score_max", "state_score_margin"],
    name="states",
)

adata.var["state_markers_selected"] = adata.var_names.isin(
    pipeline.marker_genes_from_csv("state_markers.csv")
)
pipeline.dotplot(
    adata,
    cluster_col=cluster_col,
    genes_bool_col="state_markers_selected",
    swap_axes=True,
    figure_name=f"dotplot_{cluster_col}_X_state_markers",
)

adata = pipeline.dotplot_allmarkers(
    adata,
    n_markers_per_cluster=3,
    cluster_col=cluster_col,
)

pipeline.export()
