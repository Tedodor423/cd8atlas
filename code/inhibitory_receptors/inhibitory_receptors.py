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
pipeline.export_cluster_state_summary_csv(adata, cluster_col="leiden", state_col="state_label", filename="cluster_state_summary.csv")

pipeline.export()
