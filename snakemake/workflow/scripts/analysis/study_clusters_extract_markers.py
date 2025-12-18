# script to cluster mesenchymal cell populations separately per study and extract their marker genes

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from harmony import harmonize

file_path = snakemake.input.adata
marker_output = snakemake.output.marker
cluster_output = snakemake.output.cluster


# Read in the AnnData object
adata = sc.read_h5ad(file_path)


def filter_patients_by_cell_count(adata, min_cells=10, sample_col="sample"):
    """
    Function to remove patients from adata object that have less than min cells.
    Returns filtered adata object.
    """
    # Assuming 'patient' column exists in adata.obs
    if sample_col not in adata.obs.columns:
        raise ValueError(f"The sample column {sample_col} is not found in adata.obs.")

    # Count the number of cells per patient
    patient_counts = adata.obs[sample_col].value_counts()

    # Identify patients with at least `min_cells` cells
    patients_to_keep = patient_counts[patient_counts >= min_cells].index

    # Filter the AnnData object to keep only the patients with at least `min_cells` cells
    adata_filtered = adata[adata.obs[sample_col].isin(patients_to_keep)].copy()

    return adata_filtered


# Extract mesenchymal cells
fibroblast_cells = adata[adata.obs["annotation_MOFA"] == "mesenchymal"]
del adata

# filter for samples that have at least 10 fibroblast cells
fibroblast_cells = filter_patients_by_cell_count(
    fibroblast_cells, min_cells=10, sample_col="sample"
)

# Further cell & gene filtering
sc.pp.filter_cells(fibroblast_cells, min_genes=200)
sc.pp.filter_genes(fibroblast_cells, min_cells=10)

# normalization
fibroblast_cells.layers["counts"] = fibroblast_cells.X.copy()
sc.pp.normalize_total(fibroblast_cells, target_sum=1e4)
sc.pp.log1p(fibroblast_cells)

# calculation of HVGs
sc.pp.highly_variable_genes(
    fibroblast_cells, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key="sample"
)

# calculate the proportion of samples each gene was named highly variable
fibroblast_cells.var["hvg_prop"] = (
    fibroblast_cells.var["highly_variable_nbatches"]
    / fibroblast_cells.obs["sample"].nunique()
)

num_hvg_study = fibroblast_cells.var["highly_variable"].sum()
# sort by proportion
sorted_df = fibroblast_cells.var.sort_values(by="hvg_prop", ascending=False)
top_n_values = sorted_df.head(num_hvg_study)
min_value = sorted_df.head(num_hvg_study)["hvg_prop"].min()
print(f"At least {min_value} of samples reported {num_hvg_study} hv genes")
print(f"Using {fibroblast_cells.var['highly_variable'].sum()} HVGs")


print("regressing out total counts & mitochondrial counts...")
sc.pp.regress_out(fibroblast_cells, ["nCount_RNA", "pct_counts_mt"])
# scale & center for PCA
sc.pp.scale(fibroblast_cells, max_value=10)
# PCA
sc.tl.pca(fibroblast_cells, svd_solver="arpack")

# run harmony to mix samples
Z = harmonize(fibroblast_cells.obsm["X_pca"], fibroblast_cells.obs, batch_key="sample")
fibroblast_cells.obsm[f"X_harmony_single_study"] = Z

sc.pp.neighbors(
    fibroblast_cells,
    use_rep=f"X_harmony_single_study",
    n_neighbors=10,
    n_pcs=40,
    key_added=f"neighbors_harmony_single_study",
)
sc.tl.umap(
    fibroblast_cells,
    neighbors_key=f"neighbors_harmony_single_study",
    random_state=33,
    min_dist=0.7,
)

# run leiden clusterin with several resolutions
marker_dict = {}
for resolution in [0.3, 0.5, 0.7]:
    sc.tl.leiden(
        fibroblast_cells,
        neighbors_key=f"neighbors_harmony_single_study",
        resolution=resolution,
        key_added=f"leiden_{str(resolution)}",
    )

    sc.tl.rank_genes_groups(
        fibroblast_cells, f"leiden_{str(resolution)}", method="wilcoxon", use_raw=False
    )
    marker_genes = pd.DataFrame(
        fibroblast_cells.uns["rank_genes_groups"]["names"]
    ).iloc[:100, :]
    new_colnames = [f"{col}_leiden_{resolution}" for col in marker_genes.columns]
    marker_genes.columns = new_colnames
    marker_dict[resolution] = marker_genes.T

# save marker genes per cluster
all_markers = pd.concat(marker_dict.values()).T
all_markers.to_csv(marker_output)

# save cluster results
harm = pd.DataFrame(fibroblast_cells.obsm[f"X_harmony_single_study"])
harm["id"] = fibroblast_cells.obs.index
harm[["UMAP1", "UMAP2"]] = fibroblast_cells.obsm["X_umap"]
harm[["leiden_0.3", "leiden_0.5", "leiden_0.7"]] = fibroblast_cells.obs[
    ["leiden_0.3", "leiden_0.5", "leiden_0.7"]
].to_numpy()
harm.to_csv(cluster_output)
