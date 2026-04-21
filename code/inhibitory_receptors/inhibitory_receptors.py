from cd8atlas.code import pipeline_elements as pipeline
from cd8atlas.code.pipeline_elements import bioreport

pipeline.init("inhibitory_receptors")

adata = pipeline.load_data("qc+subsampled_100000.h5ad")

adata = pipeline.label_genes(adata)
adata = pipeline.feature_selection_full(adata, include_markers="state_markers.csv")

adata = pipeline.umap_basic(
    adata,
    n_pcs=25,
    n_neighbors=12,
    min_dist=0.1,
    leiden_resolution=0.4,
    save=True,
)

adata = pipeline.label_states( adata, "state_markers.csv")  # from https://www.biorxiv.org/content/10.1101/2024.07.17.603780v1.full
pipeline.umap_plot(adata,
    umap_color_cols=["leiden", "state_label", "state_score_max", "state_score_margin"],
    name="umap_states")

adata.var["state_markers_selected"] = adata.var_names.isin(pipeline.marker_genes_from_csv("state_markers.csv"))


pipeline.dotplot(adata, cluster_col="state_label", genes_bool_col="state_markers_selected", figure_name="dotplot_state_label_X_state_markers")

pipeline.dotplot(adata, cluster_col="leiden",
    genes_bool_col="state_markers_selected",
    gene_groups=pipeline.marker_gene_groups_from_csv("state_markers.csv"),
    swap_axes=False,
    figure_name=f"dotplot_leiden_X_state_markers")

adata = pipeline.dotplot_allmarkers(adata, n_markers_per_cluster=3, cluster_col="leiden")

pipeline.export()
