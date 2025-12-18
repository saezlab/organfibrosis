import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc


# Get input and output paths from snakemake
input = snakemake.input[0]
metadata = snakemake.output[0]
pb_data = snakemake.output[1]
col_data = snakemake.output[2]
organ = snakemake.wildcards["organ"]
study = snakemake.wildcards["study"]
days_post_treatment = snakemake.params[0]
metacolumns = snakemake.params[1]

# Read AnnData object
adata = sc.read_h5ad(input)


# Subset AnnData to cells with annotation
adata = adata[adata.obs["annotation_MOFA"].notna()].copy()

# Create meta_data DataFrame with selected columns
meta_columns_sel = list(set(adata.obs.columns) & set(metacolumns))
meta_data = adata.obs[meta_columns_sel].drop_duplicates()
print(adata.obs[["sample"]].drop_duplicates().shape)

# Add grouping if missing
if "grouping" not in meta_data.columns:
    days = days_post_treatment[meta_data["study"][0]]
    meta_data["grouping"] = days

# Add organ column and save meta_data
meta_data["organ"] = organ
meta_data.to_csv(metadata)

# Compute cell-type numbers needed for MOFAcell in coldata
cell_type_numbers = adata.obs.groupby(["sample", "annotation_MOFA"])[
    "annotation_MOFA"
].count()
cell_type_numbers = (
    cell_type_numbers.to_frame()
    .rename({"annotation_MOFA": "counts"}, axis=1)
    .reset_index()
)

# OCEAN data was SoupX corrected, following steps cannot all handle floats
if study == "McCown_2025_sn":
    adata.layers["counts"] = np.round(adata.X)

# Create pseudobulk using decoupler
padata = dc.get_pseudobulk(
    adata,
    sample_col="sample",
    groups_col="annotation_MOFA",
    layer="counts",
    min_prop=0,
    min_smpls=0,
)

# Save pseudobulk data to CSV
pb_dat = pd.DataFrame(padata.X)
pb_dat.columns = padata.var.index.values
pb_dat.index = padata.obs.index.values
pb_dat.to_csv(pb_data)

# Prepare coldata for DESeq2 and save
pb_coldata = padata.obs.copy()
pb_coldata["colname"] = pb_coldata.index.values
pb_coldata = pb_coldata.merge(
    cell_type_numbers, on=["annotation_MOFA", "sample"], how="left"
)
pb_coldata.to_csv(col_data)
