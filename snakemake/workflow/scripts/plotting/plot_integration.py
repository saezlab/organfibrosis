# Script to plot label transfer results from label transfer

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import os
import pandas as pd
import decoupler as dc
import plotnine as p9
import PyComplexHeatmap as pch
from matplotlib.backends.backend_pdf import PdfPages


ref_path = snakemake.input["ref"]
query_path = snakemake.input["query"]


wildcards = snakemake.wildcards

confusion_matrix = snakemake.output["confusion_matrix"]
#scArches_ecdf = snakemake.output["scArches_ecdf"]
scArches_umap = snakemake.output["scArches_umap"]

if snakemake.wildcards['study'] == 'Wilson_2022':
    original_label = "author_cell_type"
elif snakemake.wildcards['study'] == 'Muto_2022':
    original_label = "celltype"

reference_broad = 'cell_type1'
original_sublabel = "label_before"
reference_sublabel = "cell_type2"
label = "cell_type1_pred"
sub_label = "cell_type2_pred"


def plot_confusion_matrix(df):
    col_norm_df = df.div(df.sum(axis=0), axis=1).fillna(0)
    row_norm_df = df.div(df.sum(axis=1), axis=0).fillna(0)

    titles = [
        "Origin of cells in predicted labels",
        "Destination of cells from original labels",
    ]

    figures = []

    for norm_df, title in zip([col_norm_df, row_norm_df], titles):

        # Create annotations inside the loop to avoid layout reuse issues
        orig_ha = pch.HeatmapAnnotation(
            Orig_labels=pch.anno_barplot(df.sum(axis=1)),
            legend=False,
            axis=0,
            label_side="bottom",
            label_kws={
                "rotation": -30,
                "horizontalalignment": "left",
                "verticalalignment": "top",
            },
        )

        pred_ha = pch.HeatmapAnnotation(
            Predicted=pch.anno_barplot(df.sum(axis=0)),
            legend=False,
            label_side="left",
        )

        # Make figure larger and do NOT use tight_layout=True here
        f = plt.figure(figsize=(10, 9), dpi=300)

        cm = pch.ClusterMapPlotter(
            data=norm_df,
            top_annotation=pred_ha,
            right_annotation=orig_ha,
            col_cluster=False,
            row_cluster=False,
            show_rownames=True,
            row_names_side="left",
            col_names_side="bottom",
            show_colnames=True,
            row_dendrogram=False,
            col_dendrogram=False,
            cmap="Blues",
            legend_gap=5,
            label="fraction of cells",
            plot=True,
            edgecolors="black",
        )

        f.suptitle(title, y=0.98)

        # Manually leave room for bottom x labels
        f.subplots_adjust(
            bottom=0.28,
            top=0.90,
            left=0.18,
            right=0.88,
        )

        figures.append(f)

    return figures


def pref_ecdf(data, wildcards):
    return (
        p9.ggplot(data[data["Score"] > 0.1], p9.aes(x="Score", color="celltype"))
        + p9.stat_ecdf(geom="step")
        + p9.theme_bw()
        + p9.theme(figure_size=(10, 10), dpi=300)
        + p9.facet_wrap(sub_label)
        + p9.labs(
            x="Prediction score",
            y="Cumulative distribution",
            color="reference" + "\ncell type\n",
            title="Predictions scores for predicted cell types from " + "reference",
        )
    )


def orig_ecdf(data, wildcards):
    return (
        p9.ggplot(data[data["Score"] > 0.1], p9.aes(x="Score", color="celltype"))
        + p9.stat_ecdf(geom="step")
        + p9.theme_bw()
        + p9.theme(
            figure_size=(10, 10), axis_text_x=p9.element_text(rotation=90), dpi=300
        )
        + p9.facet_wrap(sub_label)
        + p9.labs(
            x="Prediction score",
            y="Cumulative distribution",
            color="reference" + "\ncell type\n",
            title="Prediction scores for " + wildcards["study"] + " cell types",
        )
    )


def soft_ecdf(data, wildcards):
    yield orig_ecdf(data, wildcards)
    yield pref_ecdf(data, wildcards)


# Load data
query = sc.read_h5ad(query_path)
query.obs["dataset"] = wildcards["study"]

ref = sc.read_h5ad(ref_path)
ref.obs["dataset"] = "reference"

adata_full = sc.concat(
    [query, ref],
).copy()

# map sub cell types to main cell type
cell_type_unique_dict = (
    ref.obs[[reference_broad, reference_sublabel]]
    .drop_duplicates()
    .groupby(reference_broad)[reference_sublabel]
    .apply(list)
    .to_dict()
)
cell_type2_to_cell_type1_map = {
    v: k for k, values in cell_type_unique_dict.items() for v in values
}
query.obs[label] = query.obs[sub_label].map(cell_type2_to_cell_type1_map)

# Run neighbors and umap
sc.pp.neighbors(adata_full, use_rep="X_scANVI")
adata_full.obsm["X_umap_scANVI"] = sc.tl.umap(adata_full, copy=True).obsm["X_umap"]

# put back in original objects
query.obsm["X_umap_scANVI"] = adata_full[
    adata_full.obs["dataset"] == wildcards["study"]
][query.obs.index, :].obsm["X_umap_scANVI"]
ref.obsm["X_umap_scANVI"] = adata_full[adata_full.obs["dataset"] == "reference"][
    ref.obs.index, :
].obsm["X_umap_scANVI"]
#soft_adata = dc.get_acts(query, obsm_key="scArches_soft_" + sub_label)

# Plot embeddings
datasets = {"query": wildcards["study"], "ref": "reference"}
meta_column = {"reference": reference_broad, wildcards["study"]: label}
titles = {
    "reference": "Labels from" + "reference",
    wildcards["study"]: "Predictions" + "for " + wildcards["study"],
}
adatas = {"reference": ref, wildcards["study"]: query}
fig, axes = plt.subplots(1, 2, figsize=(10, 5), dpi=300)
axes = axes.flatten()
plt.suptitle("scArches embeddings")
for ax, ds in zip(axes, list(datasets.keys())[::-1]):
    background_ds = datasets["ref"] if ds == "query" else datasets["query"]

    sc.pl.embedding(
        adatas[background_ds], basis="X_umap_scANVI", frameon=False, show=False, ax=ax
    )

    sc.pl.embedding(
        adatas[datasets[ds]],
        basis="X_umap_scANVI",
        color=[meta_column[datasets[ds]]],
        frameon=False,
        title=titles[datasets[ds]],
        ax=ax,
        show=False,
        alpha=0.7,
        legend_loc="on data" if ds == "ref" else "right margin",
    )

# add title
axes[0].set_title("Reference labels from " + "reference")
axes[1].set_title("Predictions for " + wildcards["study"])

plt.savefig(scArches_umap, dpi=300)

'''
# get prediction scores
soft = pd.concat(
    [
        query.obsm["scArches_soft_" + sub_label],
        query.obs.filter([sub_label, original_sublabel], axis=1),
    ],
    axis=1,
)

# go from wide to long format
long_soft = pd.melt(
    soft.reset_index(),
    id_vars=["index", sub_label, original_sublabel],
    var_name="celltype",
    value_name="Score",
)


p9.save_as_pdf_pages(soft_ecdf(long_soft, wildcards), scArches_ecdf)
'''

# Plot confusion matrix

combinations = [
    [original_label, label],
    [original_label, sub_label],
    [original_sublabel, label],
    [original_sublabel, sub_label],
]

with PdfPages(confusion_matrix) as pdf:
    for combination in combinations:
        # compute confusion matrix between old and new labels
        df = query.obs.groupby(combination).size().unstack(fill_value=0)

        figures = plot_confusion_matrix(df)

        for f in figures:
            pdf.savefig(f, bbox_inches="tight")
            plt.close(f)
