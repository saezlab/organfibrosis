# This script identifies the factors that
# associate to the diaseas condition most signifiantly. 
# It then plots the loadings of these facotrs for each organ and study.


import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import decoupler as dc
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mofax as mfx
import pandas as pd
import seaborn as sns
import numpy as np


plt.rcParams.update({"font.size": 20})
sns.set_style("whitegrid")


def group_paths_by_organ(paths: Iterable[str], organs: Sequence[str]) -> Dict[str, List[str]]:
    """Return a dictionary {organ: [paths]} filtered by the supplied organs."""
    organ_set = set(organs)
    grouped: Dict[str, List[str]] = defaultdict(list)

    for raw_path in paths:
        path = Path(raw_path)
        parts = path.parts

        if "metadata" in parts:
            organ = parts[-3]
        else:
            organ = path.parent.name

        if organ in organ_set:
            grouped[organ].append(str(path))

    return {organ: sorted(files) for organ, files in grouped.items()}


def find_best_factors(
    loadings: pd.DataFrame, metadata: pd.DataFrame, variable: str
) -> Tuple[str, str, pd.DataFrame, pd.DataFrame]:
    """Identify the top two factors associated with a metadata variable."""
    if variable not in metadata.columns:
        raise KeyError(f"Column '{variable}' missing from metadata.")

    # Build a metadata frame aligned with the MOFA loadings index.
    aligned_metadata = metadata.loc[loadings.index, [variable]]
    associations = dc.get_metadata_associations(
        [loadings, aligned_metadata], obs_keys=[variable]
    )

    assoc_subset = (
        associations.loc[associations["variable"] == variable]
        .sort_values("p_adj")
        .copy()
    )
    if assoc_subset.empty:
        raise ValueError(f"No associations found for '{variable}'.")

    top_factors = assoc_subset["factor"].tolist()
    factor1 = top_factors[0]
    factor2 = top_factors[1] if len(top_factors) > 1 else top_factors[0]

    per_study: List[pd.DataFrame] = []
    for study, study_meta in metadata.groupby("study"):
        study_samples = study_meta.index.intersection(loadings.index)
        if len(study_samples) < 2:
            continue

        sub_loadings = loadings.loc[study_samples]
        sub_metadata = study_meta.loc[study_samples, [variable]]
        study_assoc = dc.get_metadata_associations(
            [sub_loadings, sub_metadata], obs_keys=[variable]
        )

        filtered = study_assoc[study_assoc["factor"].isin([factor1, factor2])].copy()
        if filtered.empty:
            continue
        filtered["study"] = study
        per_study.append(filtered)

    study_assoc_df = (
        pd.concat(per_study, ignore_index=True)
        if per_study
        else pd.DataFrame(columns=list(associations.columns) + ["study"])
    )

    return factor1, factor2, assoc_subset, study_assoc_df


def summarise_r2(r2_df: pd.DataFrame, factors: Sequence[str]) -> pd.DataFrame:
    """Aggregate R2 per group/view for a set of factors."""
    if not factors:
        return pd.DataFrame(columns=["Group", "View", "R2"])

    filtered = r2_df[r2_df["Factor"].isin(factors)].copy()
    if filtered.empty:
        return pd.DataFrame(columns=["Group", "View", "R2"])

    return (
        filtered.groupby(["Group", "View"], as_index=False)["R2"]
        .sum()
        .sort_values(["Group", "View"])
    )


def annotate_significance(
    ax: plt.Axes,
    data: pd.DataFrame,
    factor: str,
    study_assoc: pd.DataFrame,
    study_order: Sequence[str],
    alpha: float = 0.05,
) -> None:
    """Annotate significant study-level differences with asterisks.

    study_order ensures we respect the categorical order used by seaborn even when
    tick labels are hidden on shared axes.
    """
    if study_assoc.empty:
        return
    pval_column = "p_adj" if "p_adj" in study_assoc.columns else "pval"
    if pval_column not in study_assoc.columns:
        return

    # Seaborn shares x-axes across rows, so derive positions from the explicit order
    # rather than relying on tick labels that might be hidden.
    pos_by_label = {label: idx for idx, label in enumerate(study_order)}

    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min if y_max > y_min else 1.0
    pad = 0.05 * y_range

    lookup = (
        study_assoc[study_assoc["factor"] == factor]
        .set_index("study")[pval_column]
        .to_dict()
    )

    for study, p_value in lookup.items():
        if pd.isna(p_value) or p_value >= alpha:
            continue
        if study not in pos_by_label:
            continue

        x = pos_by_label[study]
        y = data.loc[data["study"] == study, factor].max()
        y = y + pad

        if y > ax.get_ylim()[1]:
            ax.set_ylim(ax.get_ylim()[0], y + pad)

        ax.text(x, y, "*", ha="center", va="bottom", fontsize=16, fontweight="bold")


def build_palette(
    condition_labels: Iterable[str],
    display_name: str,
    organ_color: str,
    condition_colors: Mapping[str, str],
) -> Dict[str, str]:
    """Create a color palette for scatter/box plots."""
    palette = {"control": "lightgrey"}
    palette[display_name] = organ_color

    for label in condition_labels:
        if label in palette:
            continue
        palette[label] = condition_colors.get(label, organ_color)

    return palette


def main(snakemake) -> None:  # noqa: D401 - Snakemake entrypoint
    organs: Sequence[str] = snakemake.params["organs"]
    organ_names: Mapping[str, str] = snakemake.params["organ_names"]
    organ_colors: Mapping[str, str] = snakemake.params["organ_colors"]
    condition_colors: Mapping[str, str] = snakemake.params["condition_colors"]

    model_paths = group_paths_by_organ(snakemake.input["models"], organs)
    metadata_paths = group_paths_by_organ(snakemake.input["metadata"], organs)

    # Snakemake defines per-organ R2 output targets.
    r2_outputs = {
        organ: output_path
        for organ, output_path in zip(organs, snakemake.output["r2"])
    }
    scatter_output = snakemake.output["scatter"]
    boxplot_output = snakemake.output["boxplot"]

    organ_results = {}

    for organ in organs:
        display_name = organ_names[organ]
        if organ not in model_paths:
            raise ValueError(f"No MOFA model path supplied for organ '{organ}'.")
        if organ not in metadata_paths:
            raise ValueError(f"No metadata paths supplied for organ '{organ}'.")

        model_path = model_paths[organ][0]
        model = mfx.mofa_model(model_path)

        loadings = model.get_factors(df=True)
        r2 = model.get_r2()

        # Metadata CSVs can be split per study; concatenate and deduplicate samples.
        metadata_frames = [pd.read_csv(path) for path in metadata_paths[organ]]
        metadata = (
            pd.concat(metadata_frames, ignore_index=True)
            .drop_duplicates(subset=["sample", "study"])
            .set_index("sample")
        )

        # Align samples present in both loadings and metadata
        common_samples = loadings.index.intersection(metadata.index)
        metadata_subset = (
            metadata.loc[common_samples, ["study", "cond_test"]]
            .dropna(subset=["cond_test"])
        )
        loadings = loadings.loc[metadata_subset.index]

        factor1, factor2, assoc, study_assoc = find_best_factors(
            loadings, metadata_subset, "cond_test"
        )

        # Collect significantly associated factors (fallback to the best pair).
        significant_factors = assoc.loc[assoc["p_adj"] < 0.05, "factor"].unique().tolist()
        if not significant_factors:
            significant_factors = [factor1]
            if factor2 != factor1:
                significant_factors.append(factor2)

        # Summarise explained variance for the selected factors.
        r2_summary = summarise_r2(r2, significant_factors)
        r2_summary.to_csv(r2_outputs[organ], index=False)

        # Prepare tidy dataframe for plotting downstream.
        plot_df = (
            loadings.join(metadata_subset, how="left")
            .rename(columns={"cond_test": "condition"})
            .reset_index()
            .rename(columns={"index": "sample"})
        )

        plot_df["condition_label"] = plot_df["condition"].replace(
            {"fibrosis": display_name}
        )

        palette = build_palette(
            plot_df["condition_label"].unique(), display_name, organ_colors[organ], condition_colors
        )
        study_order = plot_df["study"].drop_duplicates().tolist()

        # Cache per-organ artefacts for reuse in subsequent plots.
        organ_results[organ] = {
            "display_name": display_name,
            "loadings": plot_df,
            "factor1": factor1,
            "factor2": factor2,
            "palette": palette,
            "study_order": study_order,
            "study_assoc": study_assoc,
        }

    # Scatter plot of the two leading fibrosis-associated factors
    n_organs = len(organs)
    ncols = 2 if n_organs > 1 else 1
    nrows = math.ceil(n_organs / ncols)

    scatter_fig, scatter_axes = plt.subplots(
        nrows,
        ncols,
        figsize=(6 * ncols, 4 * nrows),
        tight_layout=True,
    )
    scatter_flat = np.array(scatter_axes).reshape(-1)

    for idx, organ in enumerate(organs):
        # Each organ panel reuses the precomputed summaries.
        result = organ_results[organ]
        ax = scatter_flat[idx]

        sns.scatterplot(
            data=result["loadings"],
            x=result["factor1"],
            y=result["factor2"],
            hue="condition_label",
            palette=result["palette"],
            s=60,
            ax=ax,
            edgecolor="black",
            linewidth=0.8,
            facecolors="none",
        )

        ax.set_title(result["display_name"])
        ax.set_xlabel(result["factor1"])
        ax.set_ylabel(result["factor2"])
        ax.get_legend().remove()

    for extra_ax in scatter_flat[n_organs:]:
        extra_ax.axis("off")

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="black",
            label="Control",
            markerfacecolor=condition_colors.get("control","lightgrey"),
            markeredgecolor="black",
            markersize=8,
            linestyle="None",
            markeredgewidth=0.8,
        )
    ]
    for organ in organs:
        display_name = organ_results[organ]["display_name"]
        legend_elements.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="black",
                label=f"Fibrosis ({display_name})",
                markerfacecolor=organ_colors[organ],
                markeredgecolor="black",
                markersize=8,
                linestyle="None",
                markeredgewidth=0.8,
            )
        )

    legend_ax = scatter_flat[min(n_organs, scatter_flat.size) - 1]
    legend_ax.legend(
        handles=legend_elements,
        title="Sample type",
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )
    scatter_fig.savefig(scatter_output, bbox_inches="tight")
    plt.close(scatter_fig)

    # Study-level boxplots for the top two fibrosis-associated factors.
    box_fig, box_axes = plt.subplots(
        2,
        n_organs,
        figsize=(4 * n_organs, 8),
        tight_layout=True,
        sharex="col",
    )
    box_axes = np.array(box_axes).reshape(2, n_organs)

    for idx, organ in enumerate(organs):
        result = organ_results[organ]
        factors = [result["factor1"], result["factor2"]]
        # Ensure a stable hue order (control first, fibrosis second).
        condition_order = [
            label for label in ["control", result["display_name"]] if label in result["palette"]
        ]
        condition_order.extend(
            [label for label in result["palette"] if label not in condition_order]
        )

        for row, factor in enumerate(factors):
            ax = box_axes[row, idx]
            sns.boxplot(
                data=result["loadings"],
                x="study",
                y=factor,
                hue="condition_label",
                palette=result["palette"],
                order=result["study_order"],
                hue_order=condition_order,
                ax=ax,
            )

            if row == 0:
                ax.set_title(result["display_name"])
            ax.set_ylabel(factor)

            if row == 1:  # bottom row
                ax.set_xlabel("Study")
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
            else:
                ax.set_xlabel("")
                ax.tick_params(labelbottom=False)  # explicitly hide labels on top row

            # Overlay study-wise significance markers for the current factor.
            annotate_significance(
                ax,
                result["loadings"],
                factor,
                result["study_assoc"],
                result["study_order"],
            )
            if ax.get_legend():
                ax.get_legend().remove()

    box_fig.savefig(boxplot_output, bbox_inches="tight")
    plt.close(box_fig)


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]
    except NameError as exc:  # pragma: no cover - Snakemake injects this variable
        raise RuntimeError(
            "This script is intended to be executed via Snakemake."
        ) from exc
    main(snakemake)  # type: ignore[name-defined]
