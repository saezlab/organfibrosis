"""Compare cross-organ LMM gene effects with MOFA factor weights and visualize their correlations."""

import pickle
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import pearsonr


# Retrieve the workflow's declared inputs, outputs, and plotting parameters.
dl_results_path = snakemake.input["dl_results"]
mofa_weights_paths = snakemake.input["mofa_weights"]

scatterplots_path = snakemake.output["scatterplots"]
heatmap_path = snakemake.output["heatmap"]
corr_df_path = snakemake.output["corr_df"]

organs = snakemake.params["organs"]
views = snakemake.params["views"]
real_names = snakemake.params["organ_names"]

plt.rcParams.update({"font.size": 18})

# Load the cross-organ LMM results and the MOFA gene weights for each organ.
with open(dl_results_path, "rb") as handle:
    dl_dict = pickle.load(handle)

mofa_weights = {
    organ: pd.read_csv(path, index_col=0)
    for organ, path in zip(organs, mofa_weights_paths)
}

corr_results = []

# Calculate correlations and collect every organ/cell-type scatterplot in one PDF.
with PdfPages(scatterplots_path) as output_pdf:
    for organ in organs:
        for ctype in views:
            df_oneview = dl_dict[ctype].reset_index()
            df_oneorgan = df_oneview[
                df_oneview["summary_row"] == real_names[organ]
            ]

            mofa_ctype = mofa_weights[organ][
                mofa_weights[organ]["celltype"] == ctype
            ]
            merged_df = df_oneorgan.merge(mofa_ctype, on="gene", how="outer")

            # The final two columns contain the two organ-specific MOFA factors.
            fig, axes = plt.subplots(2, figsize=(6, 9))

            for index, y_col in enumerate(merged_df.columns[-2:]):
                valid = merged_df[["eff", y_col]].dropna()

                if len(valid) >= 2:
                    pearson_r, p_value = pearsonr(valid["eff"], valid[y_col])
                else:
                    pearson_r, p_value = np.nan, np.nan

                corr_results.append(
                    {
                        "organ": organ,
                        "celltype": ctype,
                        "factor": y_col,
                        "pearson_r": pearson_r,
                        "p_value": p_value,
                        "n": len(valid),
                    }
                )

                sns.scatterplot(
                    ax=axes[index],
                    data=merged_df,
                    x="eff",
                    y=y_col,
                )
                axes[index].set_title(
                    f"{real_names[organ]} | {ctype} | "
                    f"{y_col.replace('_', ' ')}\n"
                    f"r = {pearson_r:.2f}, n = {len(valid)}"
                )

            fig.tight_layout()
            output_pdf.savefig(fig)
            plt.close(fig)

# Save the statistics underlying both the scatterplots and correlation heatmap.
corr_df = pd.DataFrame(corr_results)
corr_df.to_csv(corr_df_path, index=False)


def prettify_factor(name):
    """Remove the organ suffix from a factor name for the heatmap label."""
    match = re.match(r"[Ff]actor\s*(\d+)", name)
    return f"Factor {match.group(1)}" if match else name.replace("_", " ").strip()


# Summarize correlations across cell types in one heatmap panel per organ.
sns.set_style("white")
plt.rcParams.update({"font.size": 20})

fig, axes = plt.subplots(
    1,
    len(organs),
    figsize=(8, 4),
    sharey=True,
    gridspec_kw={"wspace": 0.18, "right": 0.9},
)
axes = np.atleast_1d(axes)
cbar_ax = fig.add_axes([0.925, 0.30, 0.017, 0.45])

for count, organ in enumerate(organs):
    ax = axes[count]
    heatmap_df = (
        corr_df[corr_df["organ"] == organ]
        .pivot_table(index="celltype", columns="factor", values="pearson_r")
        .reindex(views)
    )
    heatmap_df.columns = [prettify_factor(column) for column in heatmap_df.columns]

    sns.heatmap(
        heatmap_df,
        annot=True,
        fmt=".2f",
        annot_kws={"size": 16},
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        linewidths=1.2,
        linecolor="white",
        square=True,
        cbar=(count == 0),
        cbar_ax=(cbar_ax if count == 0 else None),
        ax=ax,
    )

    ax.set_title(real_names[organ], fontsize=20, fontweight="bold", pad=10)
    ax.set_xlabel("")
    ax.set_ylabel("cell type" if count == 0 else "", fontsize=20)
    ax.tick_params(axis="y", labelsize=20, rotation=0)
    ax.tick_params(axis="x", labelsize=20)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

cbar_ax.tick_params(labelsize=18)
cbar_ax.set_ylabel("pearson corr.", rotation=270, labelpad=26, fontsize=22)
fig.savefig(heatmap_path, bbox_inches="tight")
plt.close(fig)
