# Script for meta-analysis of differential expression results across organs and studies
# for disease fibroblasts vs other mesenchymal cells
# This script loads per-study DEG results, reformats them, and performs random-effects meta-analysis
# using statsmodels to combine effect sizes within organs and across organs. Results are saved as pickle files
# for downstream analysis. D

import pickle

import pandas as pd
from statsmodels.stats.meta_analysis import combine_effects

# snakemake inputs
deg_files = snakemake.input["deg_file"]

# snakemake parameters
org_colors = snakemake.params["organ_colors"]
organs = org_colors.keys()
real_names = snakemake.params["organ_names"]

# snakemake outputs
organ_spec_path = snakemake.output["organ_spec"]
cross_organ_path = snakemake.output["cross_organ"]

input_df = pd.DataFrame(columns=["organ", "study", "deg_file"])
for count, deg_file in enumerate(deg_files):
    organ = deg_file.split("/")[-2]
    study = deg_file.split("/")[-1][0:-4]
    input_df.loc[count, "organ"] = organ
    input_df.loc[count, "study"] = study
    input_df.loc[count, "deg_file"] = deg_file

all_studies = input_df["study"].unique()


def prep_deg_files_allgenes(input_df):
    """
    function retruns diff. gene expression results for a cell type of interest for all genes.
    input:
    output:
    reformated dfs:
    stat_df - pandas df with logFC results per study and gene
    se_df - pandas df with se value results per study and gene
    """
    stat_col = f"log2FoldChange"
    se_col = f"lfcSE"

    stat_df = pd.DataFrame()
    se_df = pd.DataFrame()

    for index, row in input_df.iterrows():
        study = row["study"]
        organ = row["organ"]
        deg = pd.read_csv(row["deg_file"], index_col=0)

        cols = deg.columns
        if len(cols) > 0:
            deg_subset = deg.loc[:, cols].rename(columns={stat_col: study})
            stat = deg_subset[study]
            stat_df = stat_df.merge(
                stat, left_index=True, right_index=True, how="outer"
            )

            se = deg_subset.loc[:, deg_subset.columns != study].rename(
                columns={se_col: study}
            )
            se = se[study]
            se_df = se_df.merge(se, left_index=True, right_index=True, how="outer")

        else:
            stat_df[study] = -1
            se_df[study] = -1

    stat_df = stat_df.fillna(-1)
    se_df = se_df.fillna(-1)

    return (stat_df, se_df)


def combine_gene_summary(gene, studies):
    """
    For a given gene (by gene name/index), this function:
    - Extracts the effect sizes (logFC) and standard errors (se) from the
        corresponding DataFrames.
    - Filters out studies with missing data (coded as -1).
    - If at least two studies have valid data, it performs a meta-analysis.
    - Returns the full summary_frame DataFrame from combine_effects.
    - If fewer than two studies are available, returns a DataFrame with the
        same summary_frame index and columns filled with NaNs.
    """
    # Get the gene's row from both DataFrames
    logfc_row = deg_dict["logFC"].loc[gene, studies]
    se_row = deg_dict["se"].loc[gene, studies]

    # Create a boolean mask for valid data (value != -1)
    valid_mask = (logfc_row != -1) & (se_row != -1)

    if valid_mask.sum() < 2:
        return None

    # Extract valid effect sizes and compute variances
    effect_sizes = logfc_row[valid_mask]
    variances = (se_row[valid_mask]) ** 2
    study_names = effect_sizes.index.tolist()

    # Perform the random-effects meta-analysis
    result = combine_effects(
        effect_sizes, variances, method_re="iterated", row_names=study_names
    )

    # Return the full summary frame (which includes the details for each study and overall summaries)
    return result.summary_frame()


def get_organ_info(df, organ_name):
    """
    get info from one organ from organ specific meta analysis results
    """
    org = df.reset_index().set_index("gene")
    org_effect = (
        org[org["summary_row"] == "random effect"]["eff"]
        .rename(organ_name)
        .astype("float64")
    )
    org_sd = (
        org[org["summary_row"] == "random effect"]["sd_eff"]
        .rename(organ_name)
        .astype("float64")
    )

    return org_effect, org_sd


def cross_organ_dl(gene, effect_df, se_df):
    """
    For the given gene, extract the effect sizes and standard errors,
    filter out missing values (NaN), and if at least two studies have valid data,
    perform a random-effects meta-analysis.

    Returns:
        The full summary_frame DataFrame from combine_effects.
        Returns None if fewer than two studies have data.
    """
    # Extract gene-specific rows from both dataframes
    eff_row = effect_df.loc[gene]
    se_row = se_df.loc[gene]

    # Create a boolean mask for studies with valid data.
    valid_mask = (~eff_row.isna()) & (~se_row.isna())

    # If fewer than two studies have data, skip this gene.
    if valid_mask.sum() < 2:
        return None

    # Extract valid effect sizes and compute variances from SEs.
    valid_effects = eff_row[valid_mask]
    valid_variances = (se_row[valid_mask]) ** 2
    study_names = valid_effects.index.tolist()

    # Perform the meta-analysis.
    result = combine_effects(
        valid_effects, valid_variances, method_re="iterated", row_names=study_names
    )

    # Return the full summary_frame.
    return result.summary_frame()


def flatten_concatenation(matrix):
    flat_list = []
    for row in matrix:
        flat_list += row
    return flat_list


# load data and run models

deg_dict = {}

deg_count = pd.DataFrame(index=flatten_concatenation(all_studies))

stat_df, se_df = prep_deg_files_allgenes(input_df)
deg_dict["logFC"] = stat_df
deg_dict["se"] = se_df


# run model per organ
organ_results = {}
for organ in organs:
    print(organ)
    studies = input_df[input_df["organ"] == organ]["study"].unique()
    print(studies)

    # We'll store the results in a dictionary keyed by gene name.
    results_dict = {
        gene: combine_gene_summary(gene, studies) for gene in deg_dict["se"].index
    }

    # combine all summary frames into one multi-index DataFrame.
    # The outer index will be the gene, and the inner index will be the summary row labels.
    combined_results = pd.concat(results_dict, names=["gene", "summary_row"])
    organ_results[organ] = combined_results

# save results s
with open(organ_spec_path, "wb") as fp:
    pickle.dump(organ_results, fp)

# combine organ-specific results into one score
co_results = {}

lung_effect, lung_se = get_organ_info(organ_results["HCAlung"], "lung")
heart_effect, heart_se = get_organ_info(organ_results["reheatHeart"], "heart")
kidney_effect, kidney_se = get_organ_info(organ_results["kidney"], "kidney")
liver_effect, liver_se = get_organ_info(organ_results["liver"], "liver")

effect_df = pd.concat([lung_effect, heart_effect, kidney_effect, liver_effect], axis=1)
se_df = pd.concat([lung_se, heart_se, kidney_se, liver_se], axis=1)

results_dict_cross_org = {
    gene: cross_organ_dl(gene, effect_df, se_df) for gene in effect_df.index
}
# Combine the summary frames into one DataFrame with a multi-index:
# The outer index is the gene and the inner index is the summary_frame row label.
combined_results_cross_org = pd.concat(
    results_dict_cross_org, names=["gene", "summary_row"]
)

print(combined_results_cross_org.head())

co_results["myofib"] = combined_results_cross_org

# save results
with open(cross_organ_path, "wb") as fp:
    pickle.dump(combined_results_cross_org, fp)
