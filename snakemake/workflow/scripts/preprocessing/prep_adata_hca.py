# Preprocessing script for HCA lung single-cell data
# Author: Leonie Küchenhoff

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages

# Snakemake input/output/params
input_file = snakemake.input[0]
cmap_cells = snakemake.params[0]
cmap_conditions = snakemake.params[1]
meta_columns = snakemake.params[2]
outputs = snakemake.output.adata
output_pdfs = snakemake.output.umap
qc_outputs = snakemake.output.qc_df
ctype_counts = snakemake.output.ctype_count
ctype_count_patients = snakemake.output.ctype_count_patients

# Studies will be renamed to fit to rest
study_renaming = {
    "Misharin_Budinger_2018":"Reyfman_2019",
    "Banovich_Kropski_2020":"Habermann_2020",
    "Kaminski_2020":"Adams_2020",
    "Sheppard_2020":"Tsukui_2020",
    "Schiller_2020":"Strunz_2020",
}

# Print output file paths for debugging
for output, output_pdf in zip(outputs, output_pdfs):
    study = output.split("/")[-1][:-5]
    print(output)
    print(output_pdf)

# Load and preprocess data
adata = sc.read_h5ad(input_file)
adata.obs["annotation_MOFA"] = None
adata.obs["cond_test"] = "control"

# Set raw data to count matrix
adata.X = adata.raw.X

# Subset to diseases of interest
adata_sub = adata[
    adata.obs["disease"].isin(
        [
            "normal",
            "cystic fibrosis",
            "pulmonary fibrosis",
            "hypersensitivity pneumonitis",
            "interstitial lung disease",
            "non-specific interstitial pneumonia",
            "pulmonary sarcoidosis",
        ]
    )
]
del adata

# Keep only reference samples that belong to disease studies
disease_studies = (
    adata_sub.obs[
        adata_sub.obs["disease"].isin(
            [
                "cystic fibrosis",
                "pulmonary fibrosis",
                "hypersensitivity pneumonitis",
                "interstitial lung disease",
                "non-specific interstitial pneumonia",
                "pulmonary sarcoidosis",
            ]
        )
    ]["study"]
    .unique()
    .tolist()
)
adata_sub = adata_sub[adata_sub.obs["study"].isin(disease_studies)].copy()

# Set gene names as var index instead of ensembl ids
adata_sub.var = adata_sub.var.reset_index().rename(columns={"index": "ensembl_gene_id"})
adata_sub.var.index = adata_sub.var["feature_name"].tolist()

# Calculate QC metrics
adata_sub.var["mt"] = adata_sub.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata_sub, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

# Rename columns for uniform nomenclature across atlases
adata_sub.obs["nCount_RNA"] = adata_sub.X.sum(axis=1)
adata_sub.obs["nFeature_RNA"] = np.count_nonzero(
    adata_sub.X.astype("int16").toarray(), axis=1
)
adata_sub.obs["library"] = [
    i + "_" + j
    for i, j in zip(
        adata_sub.obs["dataset"].tolist(),
        [entry[0] for entry in adata_sub.obs["disease"]],
    )
]
adata_sub.obs = adata_sub.obs.rename(
    columns={
        "ann_level_2": "cell_type1",
        "ann_level_3": "cell_type2",
        "assay": "tech",
        "age_or_mean_of_age_range": "age",
    }
)
adata_sub.obs["sex"] = adata_sub.obs["sex"].replace(
    {"male": "Male", "female": "Female"}
)
adata_sub.obs["study"] = adata_sub.obs["study"].replace(study_renaming)

print(f"adata_sub start: {adata_sub.shape}")

# Extract QC values for plotting
qc_df = adata_sub.obs[["pct_counts_mt", "nFeature_RNA", "nCount_RNA", "study"]]

# Do QC filtering
adata_sub = adata_sub[adata_sub.obs["pct_counts_mt"] < 25, :]
sc.pp.filter_cells(adata_sub, min_genes=200)
sc.pp.filter_cells(adata_sub, min_counts=500)
sc.pp.filter_cells(adata_sub, max_counts=40000)
sc.pp.filter_genes(adata_sub, min_cells=3)

# Split single nucleus from single cell Seq data
adata_sub = adata_sub[adata_sub.obs["suspension_type"] == "cell"].copy()

# Homogenize cell type labels
adata_sub.obs.loc[
    adata_sub.obs["cell_type1"].str.contains("Fibroblast", case=False),
    "annotation_MOFA",
] = "fibroblast"
adata_sub.obs.loc[
    adata_sub.obs["cell_type1"].str.contains("Smooth muscle", case=False),
    "annotation_MOFA",
] = "pericytesSMCs"
adata_sub.obs.loc[adata_sub.obs["ann_level_4"] == "Pericytes", "annotation_MOFA"] = (
    "pericytesSMCs"
)
adata_sub.obs.loc[
    adata_sub.obs["cell_type1"].str.contains("Myeloid", case=False), "annotation_MOFA"
] = "myeloid"
adata_sub.obs.loc[
    adata_sub.obs["cell_type1"].str.contains("Lymphoid", case=False), "annotation_MOFA"
] = "lymphoid"
adata_sub.obs.loc[
    adata_sub.obs["ann_level_1"].str.contains("Endo", case=False), "annotation_MOFA"
] = "endothelial"
adata_sub.obs.loc[
    adata_sub.obs["ann_level_1"].isin(["Mesothelium", "Epithelial"]), "annotation_MOFA"
] = "epithelial"
adata_sub.obs.loc[
    adata_sub.obs["cell_type1"].isin(["Mesothelium"]), "annotation_MOFA"
] = "epithelial"


print(f"adata_sub after cell types: {adata_sub.shape}")

# Merge info from ann_level_3 and 4 into one column and replace NaNs
adata_sub.obs["ann_level_4"] = adata_sub.obs["ann_level_4"].fillna("Unknown")
adata_sub[adata_sub.obs["ann_level_4"] == "None"].obs["ann_level_4"] = "Unknown"
adata_sub.obs["cell_type2"] = (
    adata_sub.obs["cell_type2"].astype(str)
    + "_"
    + adata_sub.obs["ann_level_4"].astype(str)
)

# Rename columns for uniform nomenclature across atlases
adata_sub.obs["cond_test"] = "control"
adata_sub.obs.loc[adata_sub.obs["disease"] != "normal", "cond_test"] = "fibrosis"
adata_sub.obs["grouping"] = adata_sub.obs["disease"]
adata_sub.obs["grouping"] = adata_sub.obs["grouping"].cat.rename_categories(
    {"normal": "control"}
)
adata_sub.layers["counts"] = adata_sub.X.copy()
adata_sub.obs["modality"] = "sc"

adata_sub.obs["ctype"] = adata_sub.obs["annotation_MOFA"]
adata_sub.obs.loc[
    adata_sub.obs["annotation_MOFA"].isin(["pericytesSMCs", "fibroblast"]),
    "annotation_MOFA",
] = "mesenchymal"

adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()

print(f"adata_new: {adata_new.shape}")

# --- Save outputs ---
# Save each study separately
for output, output_pdf, qc_output, ctype_count, ctype_count_patient in zip(
    outputs, output_pdfs, qc_outputs, ctype_counts, ctype_count_patients
):
    study = output.split("/")[-1][:-5]
    print("study")
    adata_study = adata_new[adata_new.obs["study"] == study].copy()
    print(f"in save loop adata_study: {adata_study.shape}")

    adata_full_study = adata_sub[adata_sub.obs["study"] == study].copy()
    print(f"in save loop adata_full_study: {adata_full_study.shape}")
    # Keep only relevant metadata columns
    cols_to_keep = list(set(adata_full_study.obs.columns) & set(meta_columns))
    adata_full_study.obs = adata_full_study.obs[cols_to_keep]

    adata_full_study.obs["study"] = study

    adata_full_study.write(output)

    # Summarize cell type counts
    counts = {
        "endothelial": [],
        "myeloid": [],
        "lymphoid": [],
        "mesenchymal": [],
        "epithelial": [],
    }
    ctypes = ["endothelial", "myeloid", "lymphoid", "epithelial", "mesenchymal"]

    for i in ctypes:
        count = adata_full_study.obs[
            adata_full_study.obs["annotation_MOFA"] == i
        ].shape[0]
        counts[i] = count
    counts["other"] = adata_full_study.obs[
        ~adata_full_study.obs["annotation_MOFA"].isin(ctypes)
    ].shape[0]
    counts["pericytesSMCs"] = adata_full_study.obs[
        adata_full_study.obs["ctype"] == "pericytesSMCs"
    ].shape[0]
    counts["fibroblast"] = adata_full_study.obs[
        adata_full_study.obs["ctype"] == "fibroblast"
    ].shape[0]
    # Creating a DataFrame from counts
    counts_df = pd.DataFrame({study: counts})
    counts_df.to_csv(ctype_count)

    # Also get counts per patient
    if pd.api.types.is_categorical_dtype(adata_full_study.obs["annotation_MOFA"]):
        adata_full_study.obs["annotation_MOFA"] = adata_full_study.obs[
            "annotation_MOFA"
        ].cat.add_categories(["other"])
    obs_data = adata_full_study.obs.copy()
    obs_data["annotation_MOFA"].fillna("other", inplace=True)
    counts = obs_data.groupby(["sample", "annotation_MOFA"]).size().reset_index()
    pivoted_df = counts.pivot(index="sample", columns="annotation_MOFA", values=0)
    counts_withcond = pivoted_df.merge(
        obs_data[["sample", "grouping", "cond_test"]], on="sample", how="left"
    ).drop_duplicates()
    counts_withcond.to_csv(ctype_count_patient)

    with PdfPages(output_pdf) as output_pdf:
        for i in [adata_full_study, adata_study]:
            fig, axs = plt.subplots(1, 3, tight_layout=True, figsize=(15, 4))
            sc.pl.umap(
                i, color=["annotation_MOFA"], ax=axs[0], show=False, palette=cmap_cells
            )
            sc.pl.umap(
                i, color=["cond_test"], ax=axs[1], show=False, palette=cmap_conditions
            )
            sc.pl.umap(i, color=["sample"], ax=axs[2], show=False, legend_loc=None)
            output_pdf.savefig(fig)

    qc_study = qc_df[qc_df["study"] == study]
    qc_study.to_csv(qc_output)
