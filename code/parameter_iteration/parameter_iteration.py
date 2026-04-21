from cd8atlas.code import pipeline_elements as pipeline
from cd8atlas.code.pipeline_elements import bioreport


HVG_N_TOP_GENES = [4000]
N_PCS = [25]
N_NEIGHBORS = [8, 10, 12]
MIN_DIST = [0.1, 0.2]
LEIDEN_RESOLUTIONS = [0.4, 0.5, 0.6]
REPORT_FIG_DPI = 60


pipeline.init("parameter_iteration")

UMAP_METADATA_COLS = ["cell_subtype_3", "disease_type", "meta_tissue_type"]
UMAP_STATE_COLS = ["state_label", "state_score_max", "state_score_margin"]

bioreport.log("Loading and preprocessing data", heading=1)

adata = pipeline.load_data("subsampled_100000.h5ad")
adata = pipeline.qc_full(adata)
adata = pipeline.normalise(adata)

bioreport.log("Hardcoded parameter grid", heading=1)
bioreport.table(
    [
        {"parameter": "hvg_n_top_genes", "values": HVG_N_TOP_GENES},
        {"parameter": "n_pcs", "values": N_PCS},
        {"parameter": "n_neighbors", "values": N_NEIGHBORS},
        {"parameter": "min_dist", "values": MIN_DIST},
        {"parameter": "leiden_resolution", "values": LEIDEN_RESOLUTIONS},
    ]
)

statistics = []
clustering_run_count = len(HVG_N_TOP_GENES) * len(N_PCS) * len(N_NEIGHBORS) * len(LEIDEN_RESOLUTIONS)
umap_plot_count = len(HVG_N_TOP_GENES) * len(N_PCS) * len(N_NEIGHBORS) * len(MIN_DIST)

bioreport.log(
    f"Running {clustering_run_count} clustering evaluations "
    f"({len(HVG_N_TOP_GENES)} HVG settings x {len(N_PCS)} PC settings x {len(N_NEIGHBORS)} neighbor settings x "
    f"{len(LEIDEN_RESOLUTIONS)} Leiden resolutions)"
)
bioreport.log(
    f"Rendering {umap_plot_count} UMAP panels for visualization "
    f"({len(MIN_DIST)} min_dist settings for each clustering graph)"
)

bioreport.log("UMAP summary grid", heading=1)

leiden_run_number = 0
for hvg_n_top_genes in HVG_N_TOP_GENES:
    bioreport.silence(True)
    adata_hvg = adata.copy()
    adata_hvg = pipeline.feature_selection_full(
        adata_hvg,
        n_top_genes=hvg_n_top_genes,
        include_markers="state_markers.csv",
    )
    top_genes = adata_hvg.var_names[
        pipeline.np.asarray(adata_hvg.layers["log1p"].mean(axis=0)).ravel().argsort()[::-1][:20]
    ].tolist()
    bioreport.silence(False)

    for n_pcs in N_PCS:
        bioreport.log(f"HVG {hvg_n_top_genes}, PCs {n_pcs}", heading=2)
        bioreport.log(
            f"Columns are n_neighbors = {', '.join(str(n) for n in N_NEIGHBORS)}. "
            f"Each panel contains all Leiden resolutions: {', '.join(str(r) for r in LEIDEN_RESOLUTIONS)}."
        )
        bioreport.silence(True)
        adata_pca = adata_hvg.copy()
        adata_pca = pipeline.pca(adata_pca, n_pcs=n_pcs)
        bioreport.silence(False)

        graph_runs: list[tuple[int, str, object, list[str]]] = []
        for n_neighbors in N_NEIGHBORS:
            base_name = f"hvg{hvg_n_top_genes}_pcs{n_pcs}_nn{n_neighbors}"
            bioreport.silence(True)
            adata_graph = adata_pca.copy()
            pipeline.sc.pp.neighbors(
                adata_graph,
                n_neighbors=n_neighbors,
                n_pcs=n_pcs,
                random_state=pipeline.random_seed,
            )

            leiden_keys = []
            for leiden_resolution in LEIDEN_RESOLUTIONS:
                leiden_run_number += 1
                leiden_key = f"leiden_res_{leiden_resolution}"
                pipeline.sc.tl.leiden(adata_graph, resolution=leiden_resolution, key_added=leiden_key)
                cluster_counts = adata_graph.obs[leiden_key].value_counts()
                leiden_keys.append(leiden_key)

                statistics.append(
                    {
                        "run_no": leiden_run_number,
                        "name": f"{base_name}_res{leiden_resolution}",
                        "base_name": base_name,
                        "n_cells": adata_graph.n_obs,
                        "n_genes": adata_graph.n_vars,
                        "hvg_n_top_genes": hvg_n_top_genes,
                        "n_pcs": n_pcs,
                        "n_neighbors": n_neighbors,
                        "leiden_resolution": leiden_resolution,
                        "n_clusters": int(cluster_counts.shape[0]),
                        "smallest_cluster_n": int(cluster_counts.min()),
                        "largest_cluster_fraction": float(cluster_counts.max() / adata_graph.n_obs),
                        "top20_highly_expressed_genes": top_genes,
                    }
                )
            bioreport.silence(False)
            graph_runs.append((n_neighbors, base_name, adata_graph, leiden_keys))

        for min_dist in MIN_DIST:
            bioreport.log(f"min_dist = {min_dist}", heading=3)
            for n_neighbors, base_name, adata_graph, leiden_keys in graph_runs:
                adata_umap = adata_graph.copy()
                pipeline.sc.tl.umap(
                    adata_umap,
                    min_dist=min_dist,
                    spread=1.0,
                    random_state=pipeline.random_seed,
                )
                bioreport.silence(True)
                adata_umap = pipeline.label_states(adata_umap, "state_markers.csv")
                bioreport.silence(False)
                umap_color_cols = leiden_keys + [
                    col for col in UMAP_METADATA_COLS if col in adata_umap.obs.columns
                ] + [
                    col for col in UMAP_STATE_COLS if col in adata_umap.obs.columns
                ]
                pipeline.sc.pl.umap(
                    adata_umap,
                    color=umap_color_cols,
                    ncols=min(3, len(umap_color_cols)),
                    wspace=0.35,
                    show=False,
                )
                bioreport.figure(
                    f"{base_name}_md{min_dist}_leiden_resolutions",
                    sameline=True,
                    dpi=REPORT_FIG_DPI,
                )

# Display stats

statistics_df = pipeline.pd.DataFrame(statistics)
csv_path = pipeline.bioreport.file_save("parameter_iteration_summary_statistics.csv")
statistics_df.to_csv(csv_path, index=False)
csv_link = csv_path.relative_to(pipeline.bioreport.report_path.parent).as_posix()

bioreport.log("Summary statistics", heading=1)
bioreport.log("These statistics summarize clustering runs only. UMAP `min_dist` is treated as a visualization parameter and is not included in the CSV or summary plots.")
bioreport.log(f"Full summary statistics CSV: [{csv_link}]({csv_link})")
bioreport.table(statistics_df.drop(columns=["top20_highly_expressed_genes"]), maxrows=50)

bioreport.log("Parameter iteration summary plots", heading=2)

pipeline.plt.figure(figsize=(7, 5))
pipeline.sns.scatterplot(
    data=statistics_df,
    x="n_clusters",
    y="largest_cluster_fraction",
    hue="leiden_resolution",
    size="n_neighbors",
)
pipeline.plt.title("Cluster count vs largest cluster fraction")
bioreport.figure("parameter_iteration_clusters_vs_dominance", dpi=REPORT_FIG_DPI)

pipeline.plt.figure(figsize=(7, 5))
cluster_heatmap = statistics_df.pivot_table(
    index="leiden_resolution",
    columns="n_neighbors",
    values="n_clusters",
    aggfunc="mean",
)
pipeline.sns.heatmap(cluster_heatmap, annot=True, fmt=".1f", cmap="viridis")
pipeline.plt.title("Mean number of clusters")
bioreport.figure("parameter_iteration_heatmap_clusters", dpi=REPORT_FIG_DPI)

pipeline.plt.figure(figsize=(7, 5))
smallest_cluster_heatmap = statistics_df.pivot_table(
    index="n_pcs",
    columns="hvg_n_top_genes",
    values="smallest_cluster_n",
    aggfunc="mean",
)
pipeline.sns.heatmap(smallest_cluster_heatmap, annot=True, fmt=".1f", cmap="magma")
pipeline.plt.title("Mean smallest cluster size")
bioreport.figure("parameter_iteration_heatmap_smallest_cluster", dpi=REPORT_FIG_DPI)

top_gene_counts = (
    statistics_df["top20_highly_expressed_genes"]
    .explode()
    .value_counts()
    .head(20)
    .reset_index()
)
top_gene_counts.columns = ["gene", "count"]

pipeline.plt.figure(figsize=(8, 6))
pipeline.sns.barplot(data=top_gene_counts, x="count", y="gene")
pipeline.plt.title("Most frequent top expressed genes across runs")
bioreport.figure("parameter_iteration_top_gene_frequency", dpi=REPORT_FIG_DPI)

pipeline.export()
