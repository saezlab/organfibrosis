# Author: Leonie Küchenhoff
# Script that predicts disease status of patients
# With diff. expressed genes of another study
# to assess the similarity of regulated genes
# between two studies
# Option A: Make assessment with top 500 regulated genes
# Option B: Make assessment with unique genes (unique for organ)

import platform as platform
import os as os
import yaml as yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import numpy as np
import scanpy as sc
from sklearn.metrics import roc_auc_score
import pickle
import upsetplot


# plot paraneters
plt.rcParams.update({"font.size": 22})

# snakemake inputs
deg_files = snakemake.input["deg_file"]
pb_files = snakemake.input["pb_file"]
metadata_files = snakemake.input["metadata_file"]
gene_dict_up_output = snakemake.input["gene_dict_up_output"]
gene_dict_down_output = snakemake.input["gene_dict_down_output"]
# snakemake parameters
org_colors = snakemake.params["organ_colors"]
organs = org_colors.keys()

category_colors = snakemake.params["celltype_colors"]
views = snakemake.params["views"]

study_colors = snakemake.params["study_colors"]
studies = study_colors.keys()

real_names = snakemake.params["organ_names"]

org_color_real = {
    real_names_key: org_colors[real_names_dict_key]
    for real_names_dict_key, real_names_key in real_names.items()
}

# snakemake outputs
predicitions_auroc = snakemake.output["predicitions_auroc"]
predicitions_clustered = snakemake.output["predicitions_clustered"]
predicitions_non_clustered = snakemake.output["predicitions_non_clustered"]
predicitions_box_withinorgan = snakemake.output["predicitions_box_withinorgan"]
predicitions_organ_agg = snakemake.output["predicitions_organ_agg"]
heatmap_output = snakemake.output["heatmap_output"]
stripplot_output = snakemake.output["stripplot_output"]

# snakemake wildcard
mode = snakemake.wildcards["mode"]


print(mode)


def prep_deg_files_allgenes(ctype):
    stat_col = f"log2FoldChange_{ctype}"
    pval_col = f"padj_{ctype}"
    tval_col = f"stat_{ctype}"

    stat_df = pd.DataFrame()
    pval_df = pd.DataFrame()
    tval_df = pd.DataFrame()

    for count, deg_file in enumerate(deg_files):
        organ = deg_file.split("/")[-2]
        study = "_".join(deg_file.split("/")[-1].split("_")[0:-3])
        deg = pd.read_csv(deg_file, index_col=0)
        cols = [col for col in deg.columns if ctype in col]
        if len(cols) > 0:
            deg_subset = deg.loc[:, cols].rename(columns={stat_col: study})
            stat = deg_subset[study]
            stat_df = stat_df.merge(
                stat, left_index=True, right_index=True, how="outer"
            )

            pval = deg_subset.loc[:, deg_subset.columns != study].rename(
                columns={pval_col: study}
            )
            pval = pval[study]
            pval_df = pval_df.merge(
                pval, left_index=True, right_index=True, how="outer"
            )

            tval = deg_subset.loc[:, deg_subset.columns != study].rename(
                columns={tval_col: study}
            )
            tval = tval[study]
            tval_df = tval_df.merge(
                tval, left_index=True, right_index=True, how="outer"
            )
        else:
            pval_df[study] = 1.0
            stat_df[study] = 0.0
            tval_df[study] = 0.0

    stat_df = stat_df.fillna(0)
    pval_df = pval_df.fillna(1)
    tval_df = tval_df.fillna(0)

    return (stat_df, pval_df, tval_df)


def read_pdata(pb_path, meta_path):
    """
    function to read in pseudobulk and metadata
    input: paths to pseudobulk and metadata
    output: anndata object with metadata in obs
    """

    metadata = pd.read_csv(meta_path)
    pdata = sc.read_csv(pb_path)
    index_names = pdata.obs.index.str.split("_")[:].tolist()
    pdata.obs["ctype"] = [item[-1] for item in index_names]
    pdata.obs["sample"] = [
        "_".join(item[:-1]) if len(item) > 1 else item[1] for item in index_names
    ]
    pdata.obs = pdata.obs.merge(metadata, on="sample", how="left")

    pdata.obs.index = pdata.obs["sample"]

    sc.pp.normalize_total(pdata, target_sum=1e4)
    sc.pp.log1p(pdata)
    sc.pp.scale(pdata)

    return pdata


input_df = pd.DataFrame(columns=["organ", "study", "deg_file"])


for count, deg_file in enumerate(deg_files):
    organ = deg_file.split("/")[-2]
    study = "_".join(deg_file.split("/")[-1].split("_")[0:-3])
    input_df.loc[count, "organ"] = organ
    input_df.loc[count, "study"] = study
    input_df.loc[count, "deg_file"] = deg_file

all_studies = input_df["study"].unique()

print(all_studies)


def get_organ_for_study(study):
    """
    function to receive organ name when study name is given
    """
    organ = input_df[input_df["study"] == study]["organ"].iloc[0]
    return organ


deg_dict = {}
for ctype in views:
    deg_dict[ctype] = {}
    stat_df, pval_df, tval_df = prep_deg_files_allgenes(ctype)

    deg_dict[ctype] = tval_df


if mode == "unique":
    with open(gene_dict_up_output, "rb") as fp:
        up = pickle.load(fp)
    with open(gene_dict_down_output, "rb") as fp:
        down = pickle.load(fp)


results_dict_A = {}


# iterate through all cell types to make predictions
for ctype in views:
    # empty df to store results
    results_df = pd.DataFrame(columns=all_studies, index=all_studies)

    deg_ctype = deg_dict[ctype]

    for train_study in all_studies:
        train_organ = get_organ_for_study(train_study)
        for metafile, pb_file in zip(metadata_files, pb_files):
            organ = pb_file.split("/")[-3]
            test_study = pb_file.split("/")[-1].split(".")[0]

            pdata = read_pdata(pb_file, metafile)

            if (mode == "unique") & (train_organ != organ):
                # pick only non-unique genes

                # get common up- and down-regulated genes into "shared" list
                up_df = upsetplot.from_contents(up[ctype])
                overlap = (
                    up_df.reset_index().set_index("id").loc[:, [train_organ, organ]]
                )
                shared = list(overlap[overlap.sum(axis=1) > 1].index)
                down_df = upsetplot.from_contents(down[ctype])
                overlap = (
                    down_df.reset_index().set_index("id").loc[:, [train_organ, organ]]
                )
                shared_down = list(overlap[overlap.sum(axis=1) > 1].index)
                shared = shared + shared_down

                # exclude shared genes from dataframe
                excluded_genes = deg_ctype.loc[~deg_ctype.index.isin(shared)]
                # pick top 500 most up or down regulated genes
                top_genes = (
                    np.abs(excluded_genes.loc[:, train_study]).nlargest(500).index
                )
            else:
                top_genes = np.abs(deg_ctype.loc[:, train_study]).nlargest(500).index

            # get overlapping genes to deg analysis
            gene_overlap = list(set(top_genes) & set(pdata.var_names))

            print(len(gene_overlap))

            # get train data (only t values from one study with relevant genes)
            deg_sub = deg_ctype.loc[gene_overlap, train_study]

            # get test data (only pseudobulks from one study with relevant genes)
            pdata_sub = pdata[pdata.obs["ctype"] == ctype][:, gene_overlap]
            pdata_vals = pdata_sub.X
            # make condition array (either 0 or 1, depending on condition)
            condition_array = np.where(pdata_sub.obs["cond_test"] == "fibrosis", 1, 0)

            # calculate dot_product
            dot_product = pdata_vals.dot(deg_sub.to_numpy())

            if dot_product.shape[0] == 0:
                print(f"test study: {test_study}, train_study: {train_study}")
                print(deg_sub.shape)
                print(pdata_sub.shape)
            else:
                if np.unique(condition_array).shape[0] > 1:
                    # calculate AUROC
                    auroc = roc_auc_score(condition_array, dot_product)

                    # save AUROC in df
                    results_df.loc[train_study, test_study] = auroc
                else:
                    print(
                        f"only on condition in test study: {test_study}, train_study: {train_study}"
                    )
                    print(deg_sub.shape)
                    print(pdata_sub.shape)

        results_df = results_df.fillna(0)
    # save results for one cell type
    results_dict_A[ctype] = results_df


# export results
with open(predicitions_auroc, "wb") as file:
    pickle.dump(results_dict_A, file)

# plot results
with PdfPages(predicitions_clustered) as output_pdf:
    for ctype in views:
        df = results_dict_A[ctype].astype(float)

        studies = df.columns
        organs_plot = [get_organ_for_study(i) for i in studies]
        organ_col = [org_colors[organ] for organ in organs_plot]

        # Create heatmaps for positive and negative Jaccard indices
        fig = plt.figure(figsize=(30, 14), tight_layout=True)
        g = sns.clustermap(
            df,
            # Turn off the clustering
            row_cluster=True,
            col_cluster=True,
            row_colors=organ_col,
            col_colors=organ_col,
            linewidths=0,
            cmap="Blues",
            vmin=0.5,
            vmax=1,
        )

        # set labels and titltes
        g.ax_heatmap.set_ylabel("train study")
        g.ax_heatmap.set_xlabel("test study")
        g.fig.suptitle(f"AUROC in {ctype} cells", size=20)

        x0, _y0, _w, _h = g.cbar_pos
        g.ax_cbar.set_position([x0 - 0.15, _y0 - 0.3, _w / 2, _h * 1])
        plt.tight_layout()

        g.ax_cbar.set_ylabel("AUROC", rotation=90, labelpad=10)
        output_pdf.savefig(bbox_inches="tight")


with PdfPages(predicitions_non_clustered) as output_pdf:
    for ctype in views:
        df = results_dict_A[ctype].astype(float)

        studies = df.columns
        organs_plot = [get_organ_for_study(i) for i in studies]
        organ_col = [org_colors[organ] for organ in organs_plot]

        # Create heatmaps for positive and negative Jaccard indices
        fig = plt.figure(figsize=(30, 14), tight_layout=True)
        g = sns.clustermap(
            df,
            # Turn off the clustering
            row_cluster=False,
            col_cluster=False,
            row_colors=organ_col,
            col_colors=organ_col,
            linewidths=0,
            cmap="Blues",
            vmin=0.5,
            vmax=1,
        )

        # set labels and titltes
        g.ax_heatmap.set_ylabel("train study")
        g.ax_heatmap.set_xlabel("test study")
        g.fig.suptitle(f"AUROC in {ctype} cells", size=20)

        x0, _y0, _w, _h = g.cbar_pos
        g.ax_cbar.set_position([x0 - 0.01, _y0 - 0.3, _w / 2, _h * 1])
        g.ax_cbar.set_ylabel("AUROC", rotation=90, labelpad=10)
        output_pdf.savefig(bbox_inches="tight")


with PdfPages(predicitions_box_withinorgan) as output_pdf:
    stripplot_dict = {}
    for ctype in views:
        # boxplots
        test = results_dict_A[ctype].melt(ignore_index=False)
        # exclude values from studies that could not be predited / predict other studies (due to not enough samples etc.)
        test_filtered = test[
            (test.index != test["variable"])
            & (test["value"] != 0)
            & (test["value"] != 0.5)
        ]
        test_filtered.index = [
            real_names[get_organ_for_study(study)] for study in test_filtered.index
        ]
        test_filtered.loc[:, "variable"] = [
            real_names[get_organ_for_study(study)]
            for study in test_filtered["variable"]
        ]
        test_filtered_twice = test_filtered[
            test_filtered.index == test_filtered["variable"]
        ]
        print(ctype)
        print(test_filtered_twice.groupby("variable").median())
        fig, ax = plt.subplots(1, figsize = (5,5))

        sns.boxplot(
            y=test_filtered_twice["value"],
            x=test_filtered_twice["variable"],
            palette=org_color_real,
            fliersize=0,
        )
        sns.stripplot(
            y=test_filtered_twice["value"],
            x=test_filtered_twice["variable"],
            palette=org_color_real,
            dodge=False,
            edgecolor="black",
            linewidth=0.5,
        )
        ax.set_title(ctype)
        ax.set_ylabel("AUROC")
        ax.set_xlabel("organ")
        ax.set_ylim(0, 1.1)
        output_pdf.savefig(bbox_inches="tight")
        stripplot_dict[ctype] = test_filtered_twice

# export results
with open(stripplot_output, "wb") as file:
    pickle.dump(stripplot_dict, file)


with PdfPages(predicitions_organ_agg) as output_pdf:
    data_heatmap_dict = {}
    for ctype in views:
        # boxplots
        test = results_dict_A[ctype].melt(ignore_index=False)
        # exclude values from studies that could not be predited / predict other studies (due to not enough samples etc.)
        test_filtered = test[
            (test.index != test["variable"])
            & (test["value"] != 0)
            & (test["value"] != 0.5)
        ]
        test_filtered.index = [
            real_names[get_organ_for_study(study)] for study in test_filtered.index
        ]
        test_filtered.loc[:, "variable"] = [
            real_names[get_organ_for_study(study)]
            for study in test_filtered["variable"]
        ]
        data_heatmap = (
            test_filtered.reset_index()
            .groupby(["variable", "index"])
            .median()
            .reset_index()
        )
        data_heatmap = data_heatmap.pivot(
            columns="variable", index="index", values="value"
        )
        data_heatmap = data_heatmap.astype(float)
        fig, ax = plt.subplots(1, figsize=(5, 4))
        sns.heatmap(data=data_heatmap, cmap="Blues", vmin=0.5, vmax=1, annot=True)
        ax.set_ylabel("train")
        ax.set_xlabel("test")
        ax.set_title(f"{ctype} \n")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right')
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha='right')
        output_pdf.savefig(bbox_inches="tight")
        data_heatmap_dict[ctype] = data_heatmap


# export results
with open(heatmap_output, "wb") as file:
    pickle.dump(data_heatmap_dict, file)
