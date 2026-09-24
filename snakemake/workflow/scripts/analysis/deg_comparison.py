# This script analyzes differential gene expression (DEG) results across multiple studies and organs,
# identifies consensus and organ-specific DEGs for various cell types, and generates summary plots
# (dotplots, barplots, UpSet plots, heatmaps) and JSON outputs for downstream interpretation.
# Author: Leonie Küchenhoff


import json
import os
import platform
import warnings

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import upsetplot
import yaml
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
from scipy import stats

# Set plot parameters and suppress warnings
plt.rcParams.update({"font.size": 18})
warnings.filterwarnings("ignore")

# Snakemake inputs

deg_files = snakemake.input["deg_file"]

# Snakemake parameters
org_colors = snakemake.params["organ_colors"]  
organs = org_colors.keys()  
category_colors = snakemake.params["celltype_colors"]  
views = snakemake.params["views"]
study_colors = snakemake.params["study_colors"] 
real_names = snakemake.params["organ_names"] 
# Snakemake outputs
organ_spec_dot = snakemake.output["organ_spec_dot"]
no_organ_consensus = snakemake.output["no_organ_consensus"]
overlap_organ_consensus = snakemake.output["overlap_organ_consensus"]
overall_consensus = snakemake.output["overall_consensus"]
deg_organ_file_up = snakemake.output["organ_consenus_up"]
deg_organ_file_down = snakemake.output["organ_consenus_down"]
common_deg_file = snakemake.output["common_deg_file"]
common_deg_file_up = snakemake.output["common_deg_file_up"]
common_deg_file_down = snakemake.output["common_deg_file_down"]
combined_pval_organ = snakemake.output["combined_pval_organ"]
deg_count_path = snakemake.output["deg_count"]
deg_sim_clust = snakemake.output["deg_sim_clust"]
deg_sim_noclust = snakemake.output["deg_sim_noclust"]

org_color_real = {
    real_names_key: org_colors[real_names_dict_key] 
    for real_names_dict_key, real_names_key in real_names.items()
}


# Build input dataframe describing all studies/organs
input_df = pd.DataFrame(columns=["organ", "study", "deg_file"])
for count, deg_file in enumerate(deg_files):
    organ = deg_file.split("/")[-2]
    study = "_".join(deg_file.split("/")[-1].split("_")[0:-3])
    input_df.loc[count, "organ"] = organ
    input_df.loc[count, "study"] = study
    input_df.loc[count, "deg_file"] = deg_file

all_studies = input_df["study"].unique()

# Utility functions
def get_organ_for_study(study):
    """
    Return organ name for a given study name.
    """
    organ = input_df[input_df["study"] == study]["organ"].iloc[0]
    return organ

ctype = "not known yet"  # Placeholder for cell type name


def plot_dotplot(
    df,
    x="study",
    y="index",
    color="FC",
    size="pval_log",
    color_range=(-1.5, 1.5),
    size_range=(0, 3),
    size_title="-log10(adj. pval)",
    color_title="log2(FC)",
    title=ctype,
):
    """
    Plot the differential expression results of a list of genes in a certain cell type.
    Args:
        df: pandas DataFrame with DEG results (must have columns for x, y, color, size)
        x: column for x-axis (default: 'study')
        y: column for y-axis (default: 'index', i.e. gene)
        color: column for dot color (default: 'FC')
        size: column for dot size (default: 'pval_log')
        color_range: tuple for color normalization
        size_range: tuple for size normalization
        size_title: legend title for size
        color_title: legend title for color
        title: plot title
    """
    if df[y].nunique() > 0:
        fig, ax = plt.subplots(
            1,
            figsize=(
                round(len(df[x].unique()) / 5) + 0.5,
                (0.15 * len(df[y].unique()) + 1.5),
            ),
            tight_layout=True,
        )
        # Scatterplot with limited size range
        sns.scatterplot(
            ax=ax,
            data=df,
            x=x,
            y=y,
            size=size,
            hue=color,
            sizes=((size_range[0] * 60) + 20, (size_range[1] * 60) + 20),
            palette="RdBu_r",
            size_norm=size_range,
            hue_norm=color_range,
        )
        yrange = [-0.5, df[y].nunique() + 2]
        ax.set_ylim(yrange)
        xrange = [-1, df[x].nunique()]
        ax.set_xlim(xrange)
        norm = plt.Normalize(color_range[0], color_range[1])
        sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm)
        sm.set_array([])
        # Custom size legend
        handles, labels = ax.get_legend_handles_labels()
        size_legend_labels = np.linspace(size_range[0], size_range[1], num=5)
        size_legend_handles = [
            plt.scatter([], [], s=(val * 60) + 20, color="black")
            for val in size_legend_labels
        ]
        ax.legend(
            size_legend_handles,
            [f"{val:.1f}" for val in size_legend_labels],
            loc="center left",
            bbox_to_anchor=(1.8, 0.5),
            title=size_title,
        )
        # Add colorbar
        fig.subplots_adjust(right=0.95)
        sub_ax = plt.axes([1, 0.35, 0.1, 0.25])
        ax.figure.colorbar(sm, label=color_title, cax=sub_ax)
        ax.tick_params(axis="x", labelsize=10)
        ax.tick_params(axis="y", labelsize=10)
        labels = ax.get_xticklabels()
        ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor")
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("")
        # Add organ color bar/markers
        organ_colors_bar = [
            org_colors.get(get_organ_for_study(study.get_text())) for study in labels
        ]
        y_max = df[y].nunique()
        ax.scatter(
            x=np.arange(len(df[x].unique())),
            y=[y_max + 1] * len(df[x].unique()),
            c=organ_colors_bar,
            s=100,
            marker="s",
            label="Organ",
        )


def deg_results_for_plot(genes, ctype):
    """
    Prepare gene expression results for plotting for a given cell type and gene list.
    Args:
        genes: list of gene names
        ctype: cell type string
    Returns:
        DataFrame with columns for plotting (FC, pval, etc.)
    """
    fc = (
        deg_dict[ctype]["logFC"]
        .loc[genes, all_studies]
        .melt(var_name="study", value_name="FC", ignore_index=False)
        .reset_index()
    )
    pval = (
        deg_dict[ctype]["p_adj"]
        .loc[genes, all_studies]
        .melt(var_name="study", value_name="pval", ignore_index=False)
        .reset_index()
    )
    fc = fc.fillna(0)
    pval = pval.fillna(1)
    plot = fc.merge(pval, on=["index", "study"], how="outer")
    plot["pval_log"] = np.log10(plot["pval"]) * -1
    plot.replace([np.inf, -np.inf], 0, inplace=True)
    return plot


def prep_deg_files_allgenes(ctype):
    """
    Return DEG results for a cell type of interest for all genes across studies.
    Args:
        ctype: cell type string
    Returns:
        stat_df: logFC per study/gene
        pval_df: adj. pval per study/gene
        tval_df: t value per study/gene
    """
    stat_col = f"log2FoldChange_{ctype}"
    pval_col = f"padj_{ctype}"
    tval_col = f"stat_{ctype}"
    stat_df = pd.DataFrame()
    pval_df = pd.DataFrame()
    tval_df = pd.DataFrame()
    for index, row in input_df.iterrows():
        study = row["study"]
        deg = pd.read_csv(row["deg_file"], index_col=0)
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


def get_summed_weights_perorgan(factor_multiplied):
    """
    For each organ, sum the weights for each gene across studies.
    Returns a dict of significant genes per organ.
    """
    test = factor_multiplied.T
    test["organ"] = [get_organ_for_study(i) for i in factor_multiplied.columns]
    sig_organ = {}
    for organ in organs:
        all_genes = test[test["organ"] == organ].iloc[:, :-1].sum(axis=0)
        sig_genes = list(all_genes[all_genes > 0.6].index)
        sig_organ[organ] = sig_genes
    return sig_organ


def get_nr_studies_inversed(df):
    """
    For each study, return 1/(number of studies in its organ).
    Used for weighting studies in consensus calculations.
    """
    nr_studies_inversed_in = []
    studies = df.columns
    organs_list = [get_organ_for_study(study) for study in studies]
    organ_nr = pd.DataFrame(organs_list).groupby(0).size()
    for organ in organs_list:
        nr_studies_inversed_in.append(1 / organ_nr[organ])
    return nr_studies_inversed_in


def check_if_gene_common(
    stat_df,
    pval_df,
    direction,
    val_cutoff_pos=0.5,
    val_cutoff_neg=-0.5,
    pval_cutoff=0.05,
):
    """
    Check if a gene is common (per organ): present in at least 60% of studies within the organ.
    Args:
        stat_df: DataFrame with test statistics
        pval_df: DataFrame with p values
        direction: 'up' or 'down' (regulation)
        val_cutoff_pos: cutoff for upregulation
        val_cutoff_neg: cutoff for downregulation
        pval_cutoff: significance threshold
    Returns:
        dict with genes per organ
    """
    nr_studies_inversed_in = get_nr_studies_inversed(stat_df)
    if direction == "up":
        significant_df = pd.DataFrame(
            (stat_df > val_cutoff_pos) & (pval_df < pval_cutoff)
        )
    elif direction == "down":
        significant_df = pd.DataFrame(
            (stat_df < val_cutoff_neg) & (pval_df < pval_cutoff)
        )
    else:
        print("Please put a valid value for direction. Either 'up' or 'down'.")
    factor_multiplied = significant_df * nr_studies_inversed_in
    summed_weights = factor_multiplied.sum(axis=1)
    gene_list = get_summed_weights_perorgan(factor_multiplied)
    return gene_list


# Main DEG comparison and plotting logic

# Reformat and save DEG results in a dictionary

deg_dict = {}
sig_genes_ctype_organ_up = {}
sig_genes_ctype_organ_down = {}
deg_count = pd.DataFrame(index=all_studies)
for ctype in views:
    deg_dict[ctype] = {}
    stat_df, pval_df, tval_df = prep_deg_files_allgenes(ctype)
    deg_dict[ctype]["logFC"] = stat_df
    deg_dict[ctype]["p_adj"] = pval_df
    deg_dict[ctype]["t_val"] = tval_df
    # Upregulated genes
    sig_genes_ctype_organ_up[ctype] = check_if_gene_common(stat_df, pval_df, "up")
    # Downregulated genes
    sig_genes_ctype_organ_down[ctype] = check_if_gene_common(stat_df, pval_df, "down")
    deg_count_ctype = pd.DataFrame(
        ((np.abs(stat_df) > 0.5) & (pval_df < 0.05)).sum(axis=0)
    ).rename(columns={0: ctype})
    deg_count = deg_count.merge(deg_count_ctype, left_index=True, right_index=True)

# Custom plotting order for studies
custom_dict = {study: nr for nr, study in enumerate(study_colors.keys())}


# Plot organ-specific genes (dotplots)
with PdfPages(organ_spec_dot) as output_pdf:
    for organ in organs:
        for ctype in views:
            genes = sig_genes_ctype_organ_up[ctype][organ]
            plot = deg_results_for_plot(genes, ctype)
            plot = plot.sort_values(
                by=["study", "index"], key=lambda x: x.map(custom_dict)
            )
            if plot.shape[0] > 0:
                fig = plot_dotplot(plot, title=f"{ctype} upregulated")
                output_pdf.savefig(bbox_inches="tight")
            genes = sig_genes_ctype_organ_down[ctype][organ]
            plot = deg_results_for_plot(genes, ctype)
            plot = plot.sort_values(
                by=["study", "index"], key=lambda x: x.map(custom_dict)
            )
            if plot.shape[0] > 0:
                fig = plot_dotplot(plot, title=f"{ctype} downregulated")
                output_pdf.savefig(bbox_inches="tight")


# Calculate and plot number of organ consensus genes
deg_count_perorgan = pd.DataFrame(columns=organs, index=views)
with PdfPages(no_organ_consensus) as output_pdf:
    for name, reg_genes in zip(
        ["up", "down"], [sig_genes_ctype_organ_up, sig_genes_ctype_organ_down]
    ):
        deg_count_perorgan = pd.DataFrame(columns=organs, index=views)
        for organ in organs:
            for ctype in views:
                deg_count_perorgan.loc[ctype, organ] = len(reg_genes[ctype][organ])
        deg_count_perorgan = deg_count_perorgan.melt(ignore_index=False).reset_index()
        fig, ax = plt.subplots(1, figsize=(5, 4))
        sns.barplot(
            data=deg_count_perorgan,
            y="value",
            x="index",
            hue="variable",
            palette=org_colors,
        )
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
        ax.set_ylabel("gene count")
        ax.set_ylabel("cell type")
        ax.set_title(f"Common {name}regulated genes per organ and cell type")
        output_pdf.savefig(bbox_inches="tight")


# UpSet plot of overlap of organ consensus genes per cell type
with PdfPages(overlap_organ_consensus) as output_pdf:
    for ctype in views:
        for name, reg_genes in zip(
            ["up", "down"], [sig_genes_ctype_organ_up, sig_genes_ctype_organ_down]
        ):
            upset_df = upsetplot.from_contents(reg_genes[ctype])
            upset_df.index.names = [
                real_names[old_name] for old_name in upset_df.index.names
            ]
            if upset_df.shape[0] > 0:
                upsetplot.plot(upset_df, show_counts=True, sort_categories_by="input", element_size = 20)
                plt.suptitle(f"{name}regulated genes in {ctype} cells", size=16)
                output_pdf.savefig(bbox_inches="tight")


# Save organ consensus genes (up/down) as JSON
with open(deg_organ_file_up, "w") as outfile:
    json.dump(sig_genes_ctype_organ_up, outfile)
with open(deg_organ_file_down, "w") as outfile:
    json.dump(sig_genes_ctype_organ_down, outfile)


# Find and plot cross-organ regulated genes (up/down/any)
for ctype in views:
    nr_studies_inversed = get_nr_studies_inversed(deg_dict[ctype]["logFC"])
    # Either up or down
    significant_df = pd.DataFrame(
        (abs(deg_dict[ctype]["logFC"]) > 0.5) & (deg_dict[ctype]["p_adj"] < 0.05)
    )[all_studies]
    factor_multiplied = significant_df * nr_studies_inversed
    summed_weights = factor_multiplied.sum(axis=1)
    # Up regulated
    significant_df_up = pd.DataFrame(
        (deg_dict[ctype]["logFC"] > 0.5) & (deg_dict[ctype]["p_adj"] < 0.05)
    )[all_studies]
    factor_multiplied_up = significant_df_up * nr_studies_inversed
    summed_weights_up = factor_multiplied_up.sum(axis=1)
    # Down regulated
    significant_df_down = pd.DataFrame(
        (deg_dict[ctype]["logFC"] < -0.5) & (deg_dict[ctype]["p_adj"] < 0.05)
    )[all_studies]
    factor_multiplied_down = significant_df_down * nr_studies_inversed
    summed_weights_down = factor_multiplied_down.sum(axis=1)
    nr_below = np.sum(significant_df, axis=1)
    deg_dict[ctype]["p_adj"]["nr_below"] = nr_below
    deg_dict[ctype]["p_adj"]["summed_weights"] = summed_weights
    deg_dict[ctype]["p_adj"]["summed_weights_up"] = summed_weights_up
    deg_dict[ctype]["p_adj"]["summed_weights_down"] = summed_weights_down
    plot_hist = deg_dict[ctype]["p_adj"].groupby("nr_below").size()
    plot_hist = pd.DataFrame(plot_hist).rename(columns={0: "count"})

# Save cross-organ regulated, up-regulated, and down-regulated genes in dictionaries
common_genes_dict = {}
common_genes_dict_up = {}
common_genes_dict_down = {}
with PdfPages(overall_consensus) as output_pdf:
    for ctype in views:
        fig, ax = plt.subplots(1)
        sns.barplot(x=plot_hist.index, y=plot_hist["count"], ax=ax, color="darkgrey")
        ax.set_yscale("log")
        ax.bar_label(ax.containers[0])
        ax.yaxis.grid(color="gray")
        ax.set_title(ctype)
        output_pdf.savefig(bbox_inches="tight")
        fig, ax = plt.subplots(1)
        sns.histplot(summed_weights, color="darkgrey")
        ax.set_yscale("log")
        ax.set_title(ctype)
        ax.set_xlabel("score")
        output_pdf.savefig(bbox_inches="tight")
        common_genes_dict[ctype] = list(
            deg_dict[ctype]["p_adj"][
                deg_dict[ctype]["p_adj"]["summed_weights"] > 1.5
            ].index
        )
        common_genes_dict_up[ctype] = list(
            deg_dict[ctype]["p_adj"][
                deg_dict[ctype]["p_adj"]["summed_weights_up"] > 1.5
            ].index
        )
        common_genes_dict_down[ctype] = list(
            deg_dict[ctype]["p_adj"][
                deg_dict[ctype]["p_adj"]["summed_weights_down"] > 1.5
            ].index
        )
        plot = deg_results_for_plot(common_genes_dict[ctype], ctype)
        plot = plot.sort_values(by=["study", "index"], key=lambda x: x.map(custom_dict))
        if plot.shape[0] > 0:
            fig = plot_dotplot(plot, title=ctype)
            output_pdf.savefig(bbox_inches="tight")
        plot = deg_results_for_plot(common_genes_dict_up[ctype], ctype)
        plot = plot.sort_values(by=["study", "index"], key=lambda x: x.map(custom_dict))
        if plot.shape[0] > 0:
            fig = plot_dotplot(plot, title=ctype)
            output_pdf.savefig(bbox_inches="tight")
        plot = deg_results_for_plot(common_genes_dict_down[ctype], ctype)
        plot = plot.sort_values(by=["study", "index"], key=lambda x: x.map(custom_dict))
        if plot.shape[0] > 0:
            fig = plot_dotplot(plot, title=ctype)
            output_pdf.savefig(bbox_inches="tight")

# Save common genes as JSON
with open(common_deg_file, "w") as outfile:
    json.dump(common_genes_dict, outfile)
with open(common_deg_file_up, "w") as outfile:
    json.dump(common_genes_dict_up, outfile)
with open(common_deg_file_down, "w") as outfile:
    json.dump(common_genes_dict_down, outfile)

# Combine p-values for each organ individually and plot
with PdfPages(combined_pval_organ) as output_pdf:
    for ctype in views:
        for organ in input_df["organ"].unique():
            organ_studies = input_df[input_df["organ"] == organ]["study"].unique()
            organ_subp = deg_dict[ctype]["p_adj"][organ_studies]
            organ_subf = deg_dict[ctype]["logFC"][organ_studies]
            # Combine p-values using Fisher's method
            sig_results = deg_dict[ctype]["p_adj"][organ_studies].apply(
                stats.combine_pvalues, axis=1
            )
            comb_pvals = [i[1] for i in sig_results]
            organ_subp["comb_pval"] = comb_pvals
            # Data for upper plot
            sorted_pvallist = organ_subp.sort_values(by=["comb_pval"])
            ranks = np.arange(0, len(sorted_pvallist), 1)
            # Extract top 500 genes and get their fold change values
            sig_genes = sorted_pvallist[organ_studies].index
            plot_FC = organ_subf.loc[sig_genes, :]
            plot_FC_top500 = plot_FC.iloc[:500, :]
            plot_FC_sorted = plot_FC_top500.loc[
                plot_FC_top500.mean(axis=1).sort_values(ascending=False).index, :
            ]
            # Plot
            fig, ax = plt.subplots(2, figsize=(13, 9), tight_layout=True, sharex=False)
            sns.scatterplot(
                x=ranks,
                y=-1 * np.log10(sorted_pvallist["comb_pval"]),
                color="black",
                linewidth=0,
                ax=ax[0],
            )
            heatmap = sns.heatmap(
                plot_FC_sorted.T,
                cmap="coolwarm",
                vmin=-2,
                vmax=2,
                ax=ax[1],
                cbar_kws={"shrink": 0.4, "aspect": 10, "label": "log2(FC)"},
            )
            ax[0].set_ylabel("-log10(BH p-value)")
            ax[0].set_xlabel("rank")
            ax[0].axvline(500, color="black")
            ax[0].set_title(ctype)
            output_pdf.savefig(bbox_inches="tight")


# Plot DEG count per organ/cell type
deg_count["organ"] = [get_organ_for_study(i) for i in deg_count.index]
plot_deg_count = deg_count.reset_index(names="study").melt(
    id_vars=["organ", "study"], value_name="count", var_name="cell type"
)
fig, ax = plt.subplots(
    1, len(views), figsize=(16, 6), sharey=True, sharex=True, tight_layout=True
)
for count, ctype in enumerate(views):
    plot_deg_count_ctype = plot_deg_count[plot_deg_count["cell type"] == ctype]
    plot_deg_count_ctype['tissue'] = [
        real_names[i] for i in plot_deg_count_ctype['organ']
        ]
    sns.boxplot(
        ax=ax[count],
        data=plot_deg_count_ctype,
        x="tissue",
        y="count",
        boxprops={"alpha": 0.4},
        palette=org_color_real,
    )
    sns.swarmplot(
        ax=ax[count],
        data=plot_deg_count_ctype,
        x="tissue",
        y="count",
        hue="study",
        palette=study_colors,
        s=8,
    )
    ax[count].legend().remove()
    ax[count].set_xticklabels(ax[count].get_xticklabels(), rotation=45, ha="right")
    ax[count].set_title(ctype)
ax[-1].legend(bbox_to_anchor=(1, 1), loc="upper left", title="study", fontsize=14)
plt.suptitle("Number of diff. expressed genes")
plt.savefig(deg_count_path, bbox_inches="tight")

plot_deg_count.to_csv(deg_count_path.replace(".pdf", ".csv"))


# Plot DEG similarity (clustered and non-clustered heatmaps)
plt.rcParams.update({"font.size": 18})
with PdfPages(deg_sim_noclust) as output_pdf_noclust:
    with PdfPages(deg_sim_clust) as output_pdf_clust:
        for ctype in views:
            ctype_tval = deg_dict[ctype]["t_val"]
            dot_ctype = ctype_tval.T.dot(ctype_tval)
            betrag = ctype_tval.apply(np.linalg.norm, axis=0)
            norm_factor = np.outer(betrag, betrag)
            norm_dot_products = dot_ctype / norm_factor
            # Custom color map for the categories
            studies = ctype_tval.columns
            organs_plot = [get_organ_for_study(i) for i in studies]
            organ_col = [org_colors[organ] for organ in organs_plot]
            # Clustered heatmap
            fig = plt.figure(figsize=(20, 10), tight_layout=True)
            g = sns.clustermap(
                norm_dot_products.fillna(0),
                row_cluster=True,
                col_cluster=True,
                row_colors=organ_col,
                col_colors=organ_col,
                linewidths=0,
                cmap="coolwarm",
                vmin=-0.8,
                vmax=0.8,
            )
            g.fig.suptitle(f"Normed dot product DGE t-values {ctype}", size=20)
            x0, _y0, _w, _h = g.cbar_pos
            g.ax_cbar.set_position([x0 - 0.1, _y0 - 0.4, _w / 2, _h * 1])
            plt.tight_layout()
            # Custom legend for organs
            legend_org_patches = [
                Patch(color=color, label=real_names[organ])
                for organ, color in org_colors.items()
            ]
            legend_org = plt.legend(
                handles=legend_org_patches,
                title="Organs",
                bbox_to_anchor=(-0.4, 1.1),
                loc="lower left",
            )
            plt.gca().add_artist(legend_org)
            output_pdf_clust.savefig(bbox_inches="tight")
            # Non-clustered heatmap
            fig = plt.figure(figsize=(30, 14), tight_layout=True)
            g = sns.clustermap(
                norm_dot_products,
                row_cluster=False,
                col_cluster=False,
                row_colors=organ_col,
                col_colors=organ_col,
                linewidths=0,
                cmap="coolwarm",
                vmin=-0.8,
                vmax=0.8,
            )
            g.fig.suptitle(f"Normed dot product DGE t-values {ctype}", size=20)
            x0, _y0, _w, _h = g.cbar_pos
            g.ax_cbar.set_position([x0 - 0.05, _y0 - 0.3, _w / 2, _h * 1])
            plt.tight_layout()
            legend_org_patches = [
                Patch(color=color, label=real_names[organ])
                for organ, color in org_colors.items()
            ]
            legend_org = plt.legend(
                handles=legend_org_patches,
                title="Organs",
                bbox_to_anchor=(-0.4, 1.1),
                loc="lower left",
            )
            plt.gca().add_artist(legend_org)
            output_pdf_noclust.savefig(bbox_inches="tight")
