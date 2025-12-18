import seaborn as sns
import numpy as np
import pandas as pd
import mofax as mfx
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch

# plot paraneters
plt.rcParams.update({"font.size": 12})


markers = snakemake.input["markers"]
models = snakemake.input["models"]

category_colors = snakemake.params["celltype_colors"]
org_colors = snakemake.params["organ_colors"]
organs = snakemake.params["organs"]
pairing = snakemake.params["pairing"]

jaccard_marker = snakemake.output["jaccard_marker"]
abs_marker = snakemake.output["abs_marker"]

jaccard_marker_csv = snakemake.output["jaccard_marker_csv"]
abs_marker_csv = snakemake.output["abs_marker_csv"]


jaccard_feature = snakemake.output["jaccard_feat"]
abs_feature = snakemake.output["abs_feat"]

jaccard_feat_csv = snakemake.output["jaccard_feat_csv"]
abs_feat_csv = snakemake.output["abs_feat_csv"]


study_organ_df = pd.DataFrame(pairing, columns=["organ", "study"])


def flatten_concatenation(matrix):
    # function to flatten a concatenated list
    flat_list = []
    for row in matrix:
        flat_list += row
    return flat_list


def get_organ_for_study(study, study_organ_df):
    # function to get organ from study name
    organ = study_organ_df[study_organ_df["study"] == study]["organ"].item()
    return organ


# Create heatmaps for Jaccard indices
def plot_heatmap(
    matrix, title, savepath, colors_x, colors_y, organ_col_x, organ_col_y, type
):
    plt.figure(figsize=(15, 15), tight_layout=True)
    if type == "jaccard":
        g = sns.clustermap(
            matrix,
            # Turn off the clustering
            row_cluster=False,
            col_cluster=False,
            row_colors=[organ_col_x, colors_x],
            col_colors=[organ_col_y, colors_y],
            linewidths=0,
            cmap="viridis",
            vmin=0,
            vmax=0.7,
            yticklabels=False,
            xticklabels=False,
        )
    if type == "absolute":
        g = sns.clustermap(
            matrix,
            # Turn off the clustering
            row_cluster=False,
            col_cluster=False,
            row_colors=[organ_col_x, colors_x],
            col_colors=[organ_col_y, colors_y],
            linewidths=0,
            cmap="viridis",
            vmin=0,
            norm=LogNorm(),
            yticklabels=False,
            xticklabels=False,
        )
    # Create custom legend for organs
    legend_org_patches = [
        Patch(color=color, label=organ) for organ, color in org_colors.items()
    ]
    legend_org = plt.legend(
        handles=legend_org_patches,
        title="Organs",
        bbox_to_anchor=(-0.2, 1.1),
        loc="lower left",
    )

    # Create custom legend for cell types
    legend_ctype_patches = [
        Patch(color=color, label=ctype) for ctype, color in category_colors.items()
    ]
    legend_ctype = plt.legend(
        handles=legend_ctype_patches,
        title="cell types",
        bbox_to_anchor=(-0.2, -0.1),
        loc="upper left",
    )

    # Add the legends to the plot
    plt.gca().add_artist(legend_org)
    plt.gca().add_artist(legend_ctype)

    # Change cbar position
    x0, _y0, _w, _h = g.cbar_pos
    g.ax_cbar.set_position([x0, _y0 - 0.5, _w / 2, _h * 1.5])
    g.ax_cbar.set_ylabel(title, rotation=90, labelpad=10)

    plt.tight_layout()
    plt.savefig(savepath, bbox_inches="tight")


# extract marker genes
marker_dict = {}
studies = []
for markerfile in markers:
    df = pd.read_csv(markerfile)
    organ = markerfile.split("/")[-3]
    study = markerfile.split("/")[-1][:-4]
    df["organ"] = organ
    df["study"] = study
    studies.append(study)
    marker_dict[study] = df


# extract features per model
feature_dict = {}
for model in models:
    study = model.split("/")[-1][:-10]
    mofa = mfx.mofa_model(model)
    features = mofa.features

    # fix annotation of features (all include celltype_ before gene name)
    for key, item in features.items():
        fixed = np.array(
            [gene.split("_", 1)[1] if "_" in gene else "" for gene in item]
        )
        features[key] = fixed

    feature_dict[study] = features


full_markers = pd.concat(marker_dict.values(), ignore_index=True)
full_markers_filtered = full_markers[
    (full_markers["FDR"] < 0.01) & (full_markers["logFC"] > 2)
]
ctypes = full_markers["name"].unique().tolist()

studies_ctypes = flatten_concatenation(
    [[i + " " + ctype for i in studies] for ctype in ctypes]
)

# Calculate Jaccard index & absolute overlap for sig marker genes
jaccard_matrix = pd.DataFrame(0, index=studies_ctypes, columns=studies_ctypes)
absolute_matrix = jaccard_matrix.copy()
categories_x = []
categories_y = []

for study1 in studies:
    # get two sets of marker genes from a certain study + cell type
    for ctype1 in ctypes:
        genes1 = full_markers_filtered[
            (full_markers_filtered["name"] == ctype1)
            & (full_markers_filtered["study"] == study1)
        ]["gene"].tolist()
        categories_x.append(ctype1)
        for study2 in studies:
            for ctype2 in ctypes:
                genes2 = full_markers_filtered[
                    (full_markers_filtered["name"] == ctype2)
                    & (full_markers_filtered["study"] == study2)
                ]["gene"].tolist()
                categories_y.append(ctype2)

                # calculate intersection & union size of to gene sets
                intersection_size = len(set(genes1) & set(genes2))
                union_size = len(genes1) + len(genes2) - intersection_size
                if union_size != 0:
                    # calculate jaccard index
                    jaccard_index = intersection_size / union_size
                else:
                    jaccard_index = 0
                # insert index into matrix
                jaccard_matrix.at[study1 + " " + ctype1, study2 + " " + ctype2] = (
                    jaccard_index
                )
                absolute_matrix.at[study1 + " " + ctype1, study2 + " " + ctype2] = (
                    intersection_size
                )

# plot results
ctypes_y = [i[-1] for i in jaccard_matrix.index.str.split(" ")]
ctypes_x = [i[-1] for i in jaccard_matrix.columns.str.split(" ")]
organs_y = [
    get_organ_for_study(i[0], study_organ_df)
    for i in jaccard_matrix.index.str.split(" ")
]
organs_x = [
    get_organ_for_study(i[0], study_organ_df)
    for i in jaccard_matrix.columns.str.split(" ")
]

# Create a custom color map for the categories
colors_x = [category_colors[category] for category in ctypes_x]
colors_y = [category_colors[category] for category in ctypes_y]

organ_col_x = [org_colors[organ] for organ in organs_x]
organ_col_y = [org_colors[organ] for organ in organs_y]


plot_heatmap(
    jaccard_matrix,
    "jaccard index",
    jaccard_marker,
    colors_x,
    colors_y,
    organ_col_x,
    organ_col_y,
    "jaccard",
)
plot_heatmap(
    absolute_matrix,
    "absolute overlap",
    abs_marker,
    colors_x,
    colors_y,
    organ_col_x,
    organ_col_y,
    "absolute",
)


jaccard_matrix.to_csv(jaccard_marker_csv)
absolute_matrix.to_csv(abs_marker_csv)

# Do the same for input genes per model

# Calculate Jaccard index for input genes
jaccard_matrix = pd.DataFrame(0, index=studies_ctypes, columns=studies_ctypes)
absolute_matrix = jaccard_matrix.copy()

print("keys:")
print(feature_dict.keys())

for study1 in studies:
    for ctype1 in ctypes:
        if ctype1 not in feature_dict[study1].keys():
            jaccard_index = 0
        else:
            genes1 = feature_dict[study1][ctype1]
            for study2 in studies:
                for ctype2 in ctypes:
                    if ctype2 in feature_dict[study2].keys():
                        genes2 = feature_dict[study2][ctype2]
                        intersection_size = len(set(genes1) & set(genes2))
                        union_size = len(genes1) + len(genes2) - intersection_size
                        jaccard_index = intersection_size / union_size
                    else:
                        jaccard_index = 0
                    jaccard_matrix.at[study1 + " " + ctype1, study2 + " " + ctype2] = (
                        jaccard_index
                    )
                    absolute_matrix.at[study1 + " " + ctype1, study2 + " " + ctype2] = (
                        intersection_size
                    )


# plot results
ctypes_y = [i[-1] for i in jaccard_matrix.index.str.split(" ")]
ctypes_x = [i[-1] for i in jaccard_matrix.columns.str.split(" ")]
organs_y = [
    get_organ_for_study(i[0], study_organ_df)
    for i in jaccard_matrix.index.str.split(" ")
]
organs_x = [
    get_organ_for_study(i[0], study_organ_df)
    for i in jaccard_matrix.columns.str.split(" ")
]


# Create a custom color map for the categories
colors_x = [category_colors[category] for category in ctypes_x]
colors_y = [category_colors[category] for category in ctypes_y]
organ_col_x = [org_colors[organ] for organ in organs_x]
organ_col_y = [org_colors[organ] for organ in organs_y]


plot_heatmap(
    jaccard_matrix,
    "jaccard index",
    jaccard_feature,
    colors_x,
    colors_y,
    organ_col_x,
    organ_col_y,
    "jaccard",
)
plot_heatmap(
    absolute_matrix,
    "absolute overlap",
    abs_feature,
    colors_x,
    colors_y,
    organ_col_x,
    organ_col_y,
    "absolute",
)

jaccard_matrix.to_csv(jaccard_feat_csv)
absolute_matrix.to_csv(abs_feat_csv)
