import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


cross_organ_dl_path = snakemake.input["cross_organ_dl"]
plot_paths = list(snakemake.output["plots"])
data_paths = list(snakemake.output["data"])

comparisons = list(snakemake.params["comparisons"])
n_top = snakemake.params["n_top"]
effect_threshold = snakemake.params["effect_threshold"]
meta_organs = list(snakemake.params["organs"])
organ_colors = dict(snakemake.params["organ_colors"])
organ_names = dict(snakemake.params["organ_names"])

organ_colors_by_name = {
    organ_names[organ]: organ_colors[organ]
    for organ in organ_names
}

plt.rcParams.update({"font.size": 20})

with open(cross_organ_dl_path, "rb") as file_handle:
    cross_organ_dl = pickle.load(file_handle)


def build_summary(celltype, organs=meta_organs):
    """Build a random-effect table with effect sizes and confidence intervals."""
    results = cross_organ_dl[celltype].reset_index().set_index("gene")
    summary = results.loc[results["summary_row"] == "random effect"].copy()

    for organ in organs:
        organ_results = results.loc[results["summary_row"] == organ]
        summary[organ] = organ_results["eff"].reindex(summary.index)
        summary[f"{organ}_ci_low"] = organ_results["ci_low"].reindex(summary.index)
        summary[f"{organ}_ci_upp"] = organ_results["ci_upp"].reindex(summary.index)

    return summary


def select_plot_data(celltype, organ1, organ2):
    """Select the top shared genes and return the values displayed in the plot."""
    summary = build_summary(celltype)
    selection_columns = [
        organ1,
        organ2,
        f"{organ1}_ci_low",
        f"{organ2}_ci_low",
    ]
    shared = summary[selection_columns].dropna()
    shared = shared[
        (shared[organ1] > effect_threshold)
        & (shared[organ2] > effect_threshold)
        & (shared[f'{organ1}_ci_low'] > 0) 
        & (shared[f'{organ2}_ci_low'] > 0)
    ]
    genes = (
        shared.mean(axis=1)
        .sort_values(ascending=False)
        .index[:n_top]
        .tolist()[::-1]
    )

    output_columns = [
        organ1,
        f"{organ1}_ci_low",
        f"{organ1}_ci_upp",
        organ2,
        f"{organ2}_ci_low",
        f"{organ2}_ci_upp",
    ]
    plot_data = summary.loc[genes, output_columns].copy()
    plot_data.index.name = "gene"
    return plot_data


def forest_plot(celltype, organ1, organ2, plot_data):
    """Plot effect sizes and confidence intervals for two organs."""
    figure_height = max(2.5, 0.5 * len(plot_data) + 1.5)
    figure, axis = plt.subplots(figsize=(7, figure_height))

    if plot_data.empty:
        axis.text(
            0.5,
            0.5,
            (
                f"No genes with effect size > {effect_threshold}\n"
                f"in both {organ1} and {organ2}"
            ),
            ha="center",
            va="center",
            transform=axis.transAxes,
        )
        axis.set_axis_off()
        axis.set_title(celltype)
        return figure

    genes = plot_data.index.tolist()
    y_positions = np.arange(len(genes))
    offsets = {organ1: 0.18, organ2: -0.18}

    for position in y_positions[::2]:
        axis.axhspan(
            position - 0.5,
            position + 0.5,
            color="0.92",
            linewidth=0,
            zorder=-2,
        )

    for organ in [organ1, organ2]:
        effects = plot_data[organ].to_numpy()
        lower = plot_data[f"{organ}_ci_low"].to_numpy()
        upper = plot_data[f"{organ}_ci_upp"].to_numpy()
        axis.errorbar(
            effects,
            y_positions + offsets[organ],
            xerr=np.vstack([effects - lower, upper - effects]),
            fmt="o",
            markersize=8,
            color=organ_colors_by_name[organ],
            ecolor=organ_colors_by_name[organ],
            elinewidth=2,
            capsize=4,
            linestyle="none",
            label=organ,
        )

    axis.axvline(0, color="grey", linestyle="--", linewidth=1, zorder=-1)
    axis.set_yticks(y_positions)
    axis.set_yticklabels(genes)
    axis.set_ylim(-0.5, len(genes) - 0.5)
    axis.set_xlabel("organ effect size")
    axis.set_title(celltype)
    axis.legend(frameon=False, fontsize=16, loc="best", bbox_to_anchor=(1, 1))
    sns.despine(ax=axis)
    return figure


if not len(comparisons) == len(plot_paths) == len(data_paths):
    raise ValueError("Each comparison must have one PDF and one CSV output.")

for comparison, plot_path, data_path in zip(comparisons, plot_paths, data_paths):
    celltype, organ1, organ2 = comparison
    plot_data = select_plot_data(celltype, organ1, organ2)

    Path(plot_path).parent.mkdir(parents=True, exist_ok=True)
    Path(data_path).parent.mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(data_path)

    figure = forest_plot(celltype, organ1, organ2, plot_data)
    figure.tight_layout()
    figure.savefig(plot_path, bbox_inches="tight")
    plt.close(figure)
