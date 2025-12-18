# This script reads in all anndata objects and filters for mesenchymal cells.
# It annotates HVGs and finally results in a merged fibroblast anndata object,
# as well as 4 lists of common HVGs - one per organ and one joint list.


import anndata as ad
import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

# Snakemake inputs
file_paths = snakemake.input
save_path = snakemake.output.adata
organs = snakemake.params["organs"]
num_hvg_study = snakemake.params["hvg"]
study_organ = snakemake.params["study_organ"]





# Functions
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


def extract_and_merge_fibroblasts(file_paths):
    """
    Function to extract fibroblasts and HVGs from each adata object.
    Returns merged adata object and one list of HVGs per study in the shape of a dictionary
    """
    # Initialize an empty list to hold the AnnData objects
    adata_list = []
    hvg_list = {}

    for file_path in file_paths:
        print(file_path)
        # Read in the AnnData object
        adata = sc.read_h5ad(file_path)

        # extract study name (later used in HVG dictionary)
        study = adata.obs["study"].iloc[0]

        # Extract fibroblast cells
        fibroblast_cells = adata[adata.obs["annotation_MOFA"].isin(["mesenchymal"])]

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
            fibroblast_cells,
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5,
            batch_key="sample",
        )

        # calculate the proportion of samples each gene was named highly variable
        fibroblast_cells.var["hvg_prop"] = (
            fibroblast_cells.var["highly_variable_nbatches"]
            / fibroblast_cells.obs["sample"].nunique()
        )

        # sort by proportion
        sorted_df = fibroblast_cells.var.sort_values(by="hvg_prop", ascending=False)
        top_n_values = sorted_df.head(num_hvg_study)
        min_value = sorted_df.head(num_hvg_study)["hvg_prop"].min()
        print(f"at least {min_value} of samples reported {num_hvg_study} hv genes")

        if "logcounts" in fibroblast_cells.layers:
            del fibroblast_cells.layers["logcounts"]

        for col in fibroblast_cells.obs.columns:
            if isinstance(fibroblast_cells.obs[col], pd.DataFrame):
                print(f"  -> .obs['{col}'] is a DataFrame!")
        for col in fibroblast_cells.var.columns:
            if isinstance(fibroblast_cells.var[col], pd.DataFrame):
                print(f"  -> .var['{col}'] is a DataFrame!")

        # Append the fibroblast cells to the list
        adata_list.append(fibroblast_cells)
        # Add HVGs to dictionary
        hvg_list[study] = top_n_values.index

        # Clear the memory
        del adata
        del fibroblast_cells

        test = ad.concat(adata_list, join="outer")
        del test

    # Concatenate all the extracted fibroblast cells into one AnnData object
    merged_adata = ad.concat(adata_list, join="outer")

    return merged_adata, hvg_list


def joint_hvg(gene_list, save_path_plot):
    """
    Function to extract shared HVGs. Extracts those genes that were named HVG most often.
    Returns list of HVGs.
    """
    # Count the frequency of each gene in 'hvg_list' and sort in decreasing order
    gene_counts = pd.Series(gene_list).value_counts().sort_values(ascending=False)

    # Create a DataFrame and group by gene count
    gene_selection_plt = pd.DataFrame(
        {"gene": gene_counts.index, "ngenes": gene_counts.values}
    )
    gene_selection_plt = gene_selection_plt.groupby("gene").sum().reset_index()

    # Create a bar plot using matplotlib
    fig, ax = plt.subplots(1)
    plt.hist(gene_selection_plt.ngenes, bins=10, edgecolor="black")
    plt.xlabel("Number of Datasets")
    plt.ylabel("Number of Genes")
    plt.savefig(save_path_plot, bbox_inches="tight")

    # make sure that no gene is selected that was reported by only one study
    gene_selection_plt[gene_selection_plt["ngenes"] > 1]
    # Select the top 'num_hvg' genes
    gene_selection = pd.DataFrame(gene_counts.index[:num_hvg_study])

    return gene_selection




print("reading in anndata objects...")
merged_adata, hvg_list = extract_and_merge_fibroblasts(file_paths)


# We have to make different HVG lists, one for all organs combined...
all_studies = merged_adata.obs["study"].unique()
genes = np.array([hvg_list[name] for name in all_studies]).flatten()
joint_gene_list = joint_hvg(genes, "plots/preprocessing/fibroblasts/hvg_all_organs.pdf")
joint_gene_list.to_csv("results/preprocessing/fibroblasts/all_organs_hvg.csv")
# ... and one per organ
for organ in organs:
    organ_studies = study_organ[study_organ["organ"] == organ]["study"].unique()
    genes = np.array([hvg_list[name] for name in organ_studies]).flatten()
    joint_gene_list = joint_hvg(
        genes, f"plots/preprocessing/fibroblasts/hvg_{organ}.pdf"
    )
    joint_gene_list.to_csv(f"results/preprocessing/fibroblasts/{organ}.csv")


# extract relevant columns in adata object
meta_columns = list(
    set(merged_adata.obs.columns)
    & set(
        [
            "identifier",
            "grouping",
            "study",
            "cond_test",
            "sample",
            "region",
            "study",
            "sex",
            "batch",
            "tech",
            "age",
            "organ",
            "fibrosis score (interstitial fibrosis) in %",
            "ischemia time in sec",
            "LVEF",
            "BMI",
            "Trichrome % fibrotic",
            "modality",
            "Fibrosis",
            "eGFR",
            "nCount_RNA",
            "pct_counts_mt",
            "annotation_MOFA",
            "cell_type1",
            "cell_type2",
        ]
    )
)


merged_adata.obs = merged_adata.obs[meta_columns]
# add info to which organ sample belong
merged_adata.obs = (
    merged_adata.obs.reset_index()
    .merge(study_organ, how="left", on="study")
    .set_index("index")
)

# change dtype of columns
merged_adata.obs["age"] = merged_adata.obs["age"].astype("float64")
merged_adata.obs["eGFR"] = merged_adata.obs["eGFR"].astype("float64")
merged_adata.obs["BMI"] = merged_adata.obs["BMI"].astype("float64")

print(f"fibroblast object will be saved to {save_path}...")
merged_adata.write(save_path)
