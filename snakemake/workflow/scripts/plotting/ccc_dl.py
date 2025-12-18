# Author: Leonie Küchenhoff
# This script plots the results of a mixed effects model
# cell-cell communication (CCC) data    

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
from matplotlib.backends.backend_pdf import PdfPages
import decoupler as dc
import numpy as np
import scanpy as sc
import json
import yaml
from statsmodels.stats.meta_analysis import combine_effects
import pickle
import liana as li
from itertools import product
import networkx as nx

# plot paraneters
plt.rcParams.update({'font.size':20})

config = snakemake.config  

organs = config["meta_organs"]
real_names = config["general_plotting"]["organ_names"]
org_colors = config["general_plotting"]["organ_colors"]
views = config['general_plotting']['views']

org_color_real = {
    real_names_key: org_colors[real_names_dict_key] 
    for real_names_dict_key, real_names_key in real_names.items()
}


ccc_result_path = snakemake.input['ccc_results']
eff_cutoff_pos = snakemake.params['eff_cutoff_pos']
eff_cutoff_neg = snakemake.params['eff_cutoff_neg']

output_pdf = snakemake.output['plots']
os.makedirs(os.path.dirname(output_pdf), exist_ok=True)
pdf = PdfPages(output_pdf)


ccc_results = pd.read_pickle(ccc_result_path)

interaction_no_dict = {}
interaction_no_dict["pos"] = {}
interaction_no_dict["neg"] = {}
for organ in organs:
    ccc_result_organ = ccc_results[organ]

    no_above_cutoff_single = ccc_result_organ[
        (ccc_result_organ["interaction_eff"] > eff_cutoff_pos)
        & (
            ccc_result_organ["interaction_eff"] - ccc_result_organ["interaction_sd_eff"]
            > 0
        )
    ]

    grouped = (
        no_above_cutoff_single.groupby(["source", "target"])
        .size()
        .reset_index()
        .rename(columns={0: "no interactions"})
    )
    grouped["organ"] = organ
    grouped["communicators"] = grouped["source"] + " -> " + grouped["target"]
    interaction_no_dict["pos"][organ] = grouped

    no_above_cutoff_single = ccc_result_organ[
        ccc_result_organ["interaction_eff"] < eff_cutoff_neg
    ]

    grouped = (
        no_above_cutoff_single.groupby(["source", "target"])
        .size()
        .reset_index()
        .rename(columns={0: "no interactions"})
    )
    grouped["organ"] = organ
    grouped["communicators"] = grouped["source"] + " -> " + grouped["target"]
    interaction_no_dict["neg"][organ] = grouped


all_interactions_pos = pd.concat(interaction_no_dict["pos"], ignore_index=True)
all_interactions_neg = pd.concat(interaction_no_dict["neg"], ignore_index=True)


fig, axs = plt.subplots(
    1, 2, figsize=(9, 5), sharex=True, sharey=True, tight_layout=True
)

sns.barplot(
    data=all_interactions_pos,
    hue="organ",
    x="source",
    y="no interactions",
    errorbar=None,
    palette=org_colors,
    ax=axs[0],
)
sns.barplot(
    data=all_interactions_neg,
    hue="organ",
    x="source",
    y="no interactions",
    errorbar=None,
    palette=org_colors,
    ax=axs[1],
)

for ax in axs:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.get_legend().remove()


axs[0].set_title("upregulated CCC events")
axs[1].set_title("downregulated CCC events")

pdf.savefig(fig, bbox_inches="tight")
plt.close(fig)


# Create barplots for the number of interactions per target cell type, split by up/down regulation
fig, axs = plt.subplots(
    1, 2, figsize=(9, 5), sharex=True, sharey=True, tight_layout=True
)

sns.barplot(
    data=all_interactions_pos,
    hue="organ",
    x="target",
    y="no interactions",
    errorbar=None,
    palette=org_colors,
    ax=axs[0],
)
sns.barplot(
    data=all_interactions_neg,
    hue="organ",
    x="target",
    y="no interactions",
    errorbar=None,
    palette=org_colors,
    ax=axs[1],
)
# Rotate x-axis labels for better readability and remove legends from both subplots
for ax in axs:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.get_legend().remove()

# Set titles for the new barplots
axs[0].set_title("upregulated CCC events")
axs[1].set_title("downregulated CCC events")

pdf.savefig(fig, bbox_inches="tight")
plt.close(fig)


# Calculate the percentage of interactions per organ for upregulated events
all_interactions_pos = all_interactions_pos.merge(
    all_interactions_pos.groupby("organ").sum()["no interactions"],
    left_on="organ",
    right_index=True,
    suffixes=["", "_sum"],
)
all_interactions_pos["percentage"] = (
    all_interactions_pos["no interactions"]
    / all_interactions_pos["no interactions_sum"]
)

# Calculate the percentage of interactions per organ for downregulated events
all_interactions_neg = all_interactions_neg.merge(
    all_interactions_neg.groupby("organ").sum()["no interactions"],
    left_on="organ",
    right_index=True,
    suffixes=["", "_sum"],
)
all_interactions_neg["percentage"] = (
    all_interactions_neg["no interactions"]
    / all_interactions_neg["no interactions_sum"]
)

layout = {
    "endothelial": np.array([1.00000000e00, 0.34202023]),
    "mesenchymal": np.array([0.12364823, 0.98480774]),
    "myeloid": np.array([-0.9396926, 0.34202023]),
    "lymphoid": np.array([-0.4999999, -0.64278773]),
    "epithelial": np.array([0.76604432, -0.64278773]),
}


for organ in organs:
    sub = all_interactions_pos[all_interactions_pos["organ"] == organ]
    fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150, tight_layout=True)
    G = nx.from_pandas_edgelist(
        sub, source="source", target="target", edge_attr=["percentage"]
    )

    G.add_nodes_from([k for k in layout if k not in G.nodes()])
    edges = G.edges()
    nodes = G.nodes()
    weights = [G[u][v]["percentage"] for u, v in edges]
    weights = np.log2(np.array(weights) + 1) * 50

    nx.draw(
        G,
        pos=layout,
        with_labels=True,
        edge_color=org_colors[organ],
        width=weights,
        ax=ax,
        alpha=1,
        node_size=1000,
        node_color="white",
        edgecolors="gray",
        font_size=18,
    )
    ax.margins(0.2)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


fig, axs = plt.subplots(
    1, 2, figsize=(9, 5), sharex=True, sharey=False, tight_layout=True
)

sns.barplot(
    data=all_interactions_pos,
    hue="organ",
    x="source",
    y="percentage",
    errorbar=None,
    palette=org_colors,
    ax=axs[0],
)
sns.barplot(
    data=all_interactions_neg,
    hue="organ",
    x="source",
    y="percentage",
    errorbar=None,
    palette=org_colors,
    ax=axs[1],
)

# Rotate x-axis labels for better readability and remove legends from both subplots
for ax in axs:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.get_legend().remove()


axs[0].set_title("upregulated CCC events")
axs[1].set_title("downregulated CCC events")

pdf.savefig(fig, bbox_inches="tight")
plt.close(fig)


fig, axs = plt.subplots(
    1, 2, figsize=(9, 5), sharex=True, sharey=False, tight_layout=True
)

sns.barplot(
    data=all_interactions_pos,
    hue="organ",
    x="target",
    y="percentage",
    errorbar=None,
    palette=org_colors,
    ax=axs[0],
)
sns.barplot(
    data=all_interactions_neg,
    hue="organ",
    x="target",
    y="percentage",
    errorbar=None,
    palette=org_colors,
    ax=axs[1],
)

# Rotate x-axis labels for better readability and remove legends from both subplots
for ax in axs:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.get_legend().remove()


axs[0].set_title("upregulated CCC events")
axs[1].set_title("downregulated CCC events")

pdf.savefig(fig, bbox_inches="tight")
plt.close(fig)

interactions = np.sort(all_interactions_pos['communicators'].unique())

fig, axs = plt.subplots(
    1, 4, figsize=(35, 7), sharex=True, sharey=True, tight_layout=True
)
for count, organ in enumerate(organs):
    sub = all_interactions_pos[all_interactions_pos["organ"] == organ]
    sns.barplot(
        data=sub,
        x="communicators",
        y="percentage",
        errorbar=None,
        color=org_colors[organ],
        ax=axs[count],
        order = interactions
    )
    axs[count].set_xticklabels(axs[count].get_xticklabels(), rotation=45, ha="right")

pdf.savefig(fig, bbox_inches="tight")
plt.close(fig)


counts_dict_s = {}
counts_dict_t = {}
for organ in organs:
    ccc_result_organ = ccc_results[organ]
    sig = ccc_result_organ[(ccc_result_organ['interaction_eff'] > eff_cutoff_pos) & (ccc_result_organ['interaction_eff'] - ccc_result_organ['interaction_sd_eff']> 0)]
    heatmap_target = sig.groupby(['interaction','target']).count()['ligand'].reset_index().pivot(columns = 'target', values = 'ligand', index = 'interaction')
    heatmap_source = sig.groupby(['interaction','source']).count()['ligand'].reset_index().pivot(columns = 'source', values = 'ligand', index = 'interaction')

    counts_dict_s[organ] = heatmap_source
    counts_dict_t[organ] = heatmap_target


# Ensure all DataFrames have the same structure
dfs_s = [counts_dict_s[organ] for organ in organs]
dfs_t = [counts_dict_t[organ] for organ in organs]

# Outer join to align indices across all DataFrames
merged_df = pd.concat(dfs_s, axis=0, join='outer').groupby(level=0).sum(min_count=1)

# Create count DataFrames to track contributions
count_df1 = counts_dict_s['reheatHeart'].sum(axis = 1)
count_df2 = counts_dict_s['HCAlung'].sum(axis = 1)
count_df3 = counts_dict_s['kidney'].sum(axis = 1)
count_df4 = counts_dict_s['liver'].sum(axis = 1)
# Fill NaN values with 0 to ensure proper summation
sum_df = pd.concat([df.fillna(0) for df in dfs_s], axis=0, join='outer').groupby(level=0).sum()

sum_df_t = pd.concat([df.fillna(0) for df in dfs_t], axis=0, join='outer').groupby(level=0).sum()

# Add contribution count columns
sum_df['heart'] = count_df1
sum_df['lung'] = count_df2
sum_df['kidney'] = count_df3
sum_df['liver'] = count_df4

sum_df = sum_df.fillna(0)
sum_df['organ_contibutions'] = np.count_nonzero(sum_df.loc[:, ['heart','lung','kidney', 'liver']], axis = 1)

subset = sum_df[sum_df['organ_contibutions'] > 3]
interactions_to_plot = subset.index 


fig, ax = plt.subplots(1, 3, sharey = True, figsize = (9, subset.shape[0] * 0.5))


sns.heatmap(sum_df_t.loc[interactions_to_plot, views], ax = ax[1], cbar=None, annot=True)
sns.heatmap(subset.loc[:, views], ax = ax[0], cbar=None, annot=True)
ax[2] = plt.barh(
    y = subset.index, width = subset['heart'], 
    color = org_color_real['heart'], align = 'edge')
ax[2] = plt.barh(
    y = subset.index, width = subset['lung'], left = subset['heart'], 
    color = org_color_real['lung'], align = 'edge')
ax[2] = plt.barh(
    y = subset.index, width = subset['kidney'], left = subset['heart'] +  subset['lung'], 
    color = org_color_real['kidney'], align = 'edge')
ax[2] = plt.barh(
    y = subset.index, width = subset['liver'], left = subset['heart'] +  subset['lung'] + subset['kidney'], 
    color = org_color_real['liver'], align = 'edge')

pdf.savefig(fig, bbox_inches="tight")
plt.close(fig)

pdf.close()
plt.close("all")

