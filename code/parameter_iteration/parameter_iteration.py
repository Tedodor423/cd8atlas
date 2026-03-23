from cd8atlas.code import pipeline_elements as pipeline_elements
from pipeline_elements import Path, pd, bioreport

pipeline_elements.init("parameter_iteration")

SPREAD = 1.0

bioreport.log("Loading and preprocessing data", heading=1)

adata = pipeline_elements.load_data("subsampled_100000.h5ad")
# pipeline_elements.load_data("huARdb_v2_GEX.CD8.all_genes.h5ad")
# pipeline_elements.subsample_data(n_obs=100000)

adata = pipeline_elements.qc_full(adata)
adata = pipeline_elements.normalise(adata)
pipeline_elements.save_data(adata, "parameter-iteration_preprocessed.h5ad", root=True)

# get parameter combinations
param_csv = Path(pipeline_elements.root_dir / "code/parameter_iteration/parameters.csv")
params_table = pd.read_csv(param_csv)
param_combinations = [
    {
            "hvg_n_top_genes": int(hvg_n_top_genes),
            "n_pcs": int(n_pcs),
            "n_neighbors": int(n_neighbors),
            "leiden_resolution": float(leiden_resolution),
            "min_dist": float(min_dist),
        }
        for hvg_n_top_genes in params_table["hvg_n_top_genes"].dropna().unique()
        for n_pcs in params_table["n_pcs"].dropna().unique()
        for n_neighbors in params_table["n_neighbors"].dropna().unique()
        for leiden_resolution in params_table["leiden_resolution"].dropna().unique()
        for min_dist in params_table["min_dist"].dropna().unique()
    ]

bioreport.log(f"Loaded parameters from {param_csv}:", heading=1)
bioreport.table(params_table)
statistics = []

# main run
bioreport.log(f"Running with {len(param_combinations)} parameter combinations")
for run_number, combination in enumerate(param_combinations, start=1):
    hvg_n_top_genes = combination["hvg_n_top_genes"]
    n_pcs = combination["n_pcs"]
    n_neighbors = combination["n_neighbors"]
    leiden_resolution = combination["leiden_resolution"]
    min_dist = combination["min_dist"]

    name = f"hvg{hvg_n_top_genes}_pcs{n_pcs}_nn{n_neighbors}_res{leiden_resolution}_md{min_dist}"

    pipeline_elements.bioreport.log(f"Run {run_number}/{len(param_combinations)}: {name}", heading=2)

    adata_run = adata.copy()
    adata_run = pipeline_elements.feature_selection_full(adata_run, n_top_genes=hvg_n_top_genes)
    adata_run = pipeline_elements.umap_clustering(
        adata_run,
        n_neighbors=n_neighbors
        min_dist=min_dist,
        spread=SPREAD,
        leiden_resolution=leiden_resolution,
        figure_name=f"umap_{name}",
        output_file=f"adata_{name}.h5ad",
    )

    statistics.append(
        {
            "run": run_number,
            "name": name,
            **combination,
            "n_obs": adata_run.n_obs,
            "n_vars": adata_run.n_vars,
        }
    )

pipeline_elements.bioreport.table(statistics)
pipeline_elements.export()