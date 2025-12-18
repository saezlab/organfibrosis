
import json
from collections.abc import Mapping
from itertools import product
from math import ceil
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from typing import Dict, List, Optional


config = snakemake.config  # type: ignore[name-defined]

organs = config["meta_organs"]
real_names = config["general_plotting"]["organ_names"]
org_colors = config["general_plotting"]["organ_colors"]

ccc_path = Path(snakemake.input.ccc)  # type: ignore[name-defined]
cosine_sim_input = snakemake.input.cosine_sim  # type: ignore[name-defined]

plt.rcParams.update({"font.size": 20})

if isinstance(cosine_sim_input, Mapping):
    cosine_sim_paths = {organ: Path(path) for organ, path in cosine_sim_input.items()}
else:
    cosine_sim_paths = {
        organ: Path(path) for organ, path in zip(organs, cosine_sim_input)
    }

per_interaction_output_path = Path(
    snakemake.output.per_interaction  # type: ignore[name-defined]
)
merged_output_path = Path(snakemake.output.merged_matrix)  # type: ignore[name-defined]
top_genes_output_path = Path(snakemake.output.top_genes)  # type: ignore[name-defined]
top_interactors_output_path = Path(
    snakemake.output.top_interactors  # type: ignore[name-defined]
)
top_genes_heatmap_path = Path(
    snakemake.output.top_genes_heatmap  # type: ignore[name-defined]
)
top_interactors_plot_path = Path(
    snakemake.output.top_interactors_plot  # type: ignore[name-defined]
)
top_genes_boxplot_path = Path(
    snakemake.output.top_genes_boxplot  # type: ignore[name-defined]
)
top_genes_heatmap_csv_path = Path(
    snakemake.output.top_genes_heatmap_csv  # type: ignore[name-defined]
)
top_interactors_plot_csv_path = Path(
    snakemake.output.top_interactors_plot_csv  # type: ignore[name-defined]
)
top_genes_boxplot_csv_path = Path(
    snakemake.output.top_genes_boxplot_csv  # type: ignore[name-defined]
)

def split_complex(value):
    if isinstance(value, str) and value:
        return value.split("_")
    return [value]


def expand_complex_interactions(df: pd.DataFrame) -> pd.DataFrame:
    expanded_rows = []
    for _, row in df.iterrows():
        ligand_candidates = (
            split_complex(row["ligand_complex"])
            if "ligand_complex" in row.index and pd.notna(row["ligand_complex"])
            else split_complex(row["ligand"])
        )
        receptor_candidates = (
            split_complex(row["receptor_complex"])
            if "receptor_complex" in row.index and pd.notna(row["receptor_complex"])
            else split_complex(row["receptor"])
        )
        if len(ligand_candidates) > 1 or len(receptor_candidates) > 1:
            for ligand, receptor in product(ligand_candidates, receptor_candidates):
                new_row = row.copy()
                new_row["interaction"] = f"{ligand}^{receptor}"
                if "ligand" in new_row.index:
                    new_row["ligand"] = ligand
                if "receptor" in new_row.index:
                    new_row["receptor"] = receptor
                expanded_rows.append(new_row)
        else:
            expanded_rows.append(row.copy())
    return pd.DataFrame(expanded_rows).reset_index(drop=True)


def test_significance_with_fdr(
    cosine_sim_dict: Dict[str, pd.DataFrame],
    organ_list: List[str],
    interactions_to_plot: List[str],
) -> Dict[str, Dict[str, Dict[str, Optional[float]]]]:
    raw_results: List[Optional[float]] = []
    test_keys: List[tuple[str, str]] = []

    for organ in organ_list:
        data = cosine_sim_dict.get(organ, pd.DataFrame())
        if data.empty:
            continue
        for interaction in interactions_to_plot:
            group_data = data[data["interaction"] == interaction]
            control_vals = group_data[group_data["cond_test"] == "control"]["mean"].dropna()
            fibrosis_vals = group_data[group_data["cond_test"] == "fibrosis"]["mean"].dropna()
            if len(control_vals) > 0 and len(fibrosis_vals) > 0:
                _, pval = mannwhitneyu(
                    control_vals,
                    fibrosis_vals,
                    alternative="two-sided",
                )
                raw_results.append(pval)
                test_keys.append((organ, interaction))
            else:
                raw_results.append(None)
                test_keys.append((organ, interaction))

    valid_indices = [i for i, p in enumerate(raw_results) if p is not None]
    if valid_indices:
        valid_pvals = [raw_results[i] for i in valid_indices]
        reject, corrected_pvals, _, _ = multipletests(
            valid_pvals, alpha=0.05, method="fdr_bh"
        )
    else:
        reject = []
        corrected_pvals = []

    corrected_dict: Dict[str, Dict[str, Dict[str, Optional[float]]]] = {}
    corrected_idx = 0

    for idx, (organ, interaction) in enumerate(test_keys):
        corrected_dict.setdefault(organ, {})
        raw_p = raw_results[idx]
        if raw_p is not None and valid_indices:
            corrected_dict[organ][interaction] = {
                "raw_p": raw_p,
                "corrected_p": corrected_pvals[corrected_idx],
                "significant": bool(reject[corrected_idx]),
            }
            corrected_idx += 1
        else:
            corrected_dict[organ][interaction] = {
                "raw_p": raw_p,
                "corrected_p": None,
                "significant": False,
            }

    return corrected_dict


ccc = pd.read_pickle(ccc_path)
condition_map = {
    "PSC": "fibrosis",
    "CKD": "fibrosis",
    "FZ": "fibrosis",
    "IPF": "fibrosis",
    "CTRL": "control",
    "Healthy": "control",
    "normal": "control",
}

per_organ_tables: Dict[str, pd.DataFrame] = {}
per_organ_rows = []
cosine_sim_full_dict: Dict[str, pd.DataFrame] = {}

for organ in organs:
    cosine_sim = pd.read_csv(cosine_sim_paths[organ], index_col=0)
    cosine_sim["organ"] = organ
    cosine_sim["cond_test"] = cosine_sim["cond_test"].replace(condition_map)
    cosine_sim = cosine_sim[cosine_sim["cond_test"].isin(["fibrosis", "control"])]
    cosine_sim_full_dict[organ] = cosine_sim.copy()

    cosine_sim_fibrosis = cosine_sim[cosine_sim["cond_test"] == "fibrosis"].copy()

    spatial_summary = (
        cosine_sim_fibrosis.groupby("interaction", as_index=False)
        .agg({"ligand": "first", "receptor": "first", "morans": "mean"})
        .rename(columns={"morans": "morans_mean"})
    )
    if not spatial_summary.empty:
        idx = spatial_summary.groupby("interaction")["morans_mean"].idxmax()
        spatial_summary = spatial_summary.loc[idx].copy()
        spatial_summary["spatial"] = spatial_summary["morans_mean"].rank(pct=True)

    ccc_organ = ccc[organ]
    idx = ccc_organ.groupby("interaction")["ligand_eff"].idxmax()
    ccc_top = ccc_organ.loc[idx].copy()
    ccc_top["expression"] = ccc_top["interaction_eff"].rank(pct=True)

    ccc_expanded = expand_complex_interactions(ccc_top)
    idx = ccc_expanded.groupby("interaction")["ligand_eff"].idxmax()
    ccc_expanded = ccc_expanded.loc[idx].copy()

    merge_cols = ["ligand", "receptor"]
    merged = spatial_summary.merge(
        ccc_expanded, how="left", on=merge_cols, suffixes=("", "_expr")
    )
    merged["mean_rank_organ"] = merged[["spatial", "expression"]].mean(axis=1)
    merged["organ"] = organ
    merged["interaction"] = merged["interaction"].fillna(
        merged["ligand"] + "^" + merged["receptor"]
    )

    columns_for_output: List[str] = [
        "organ",
        "interaction",
        "ligand",
        "receptor",
        "morans_mean",
        "spatial",
        "expression",
        "interaction_eff",
        "ligand_eff",
        "mean_rank_organ",
    ]
    available_columns = [col for col in columns_for_output if col in merged.columns]
    per_organ_table = merged[available_columns].set_index("interaction")
    per_organ_tables[organ] = per_organ_table
    per_organ_rows.append(per_organ_table.reset_index())

per_organ_df = (
    pd.concat(per_organ_rows, ignore_index=True) if per_organ_rows else pd.DataFrame()
)
per_interaction_output_path.parent.mkdir(parents=True, exist_ok=True)
per_organ_df.to_csv(per_interaction_output_path, index=False)

matrix_parts = {
    organ: table[["spatial", "expression"]].add_suffix(f" {real_names[organ]}")
    for organ, table in per_organ_tables.items()
    if not table.empty
}
if matrix_parts:
    merged_matrix = pd.concat(matrix_parts, axis=1)
    merged_matrix.columns = merged_matrix.columns.droplevel(0)
    merged_matrix["mean"] = merged_matrix.mean(axis=1)
else:
    merged_matrix = pd.DataFrame()

merged_output_path.parent.mkdir(parents=True, exist_ok=True)
merged_matrix.to_csv(merged_output_path)

heatmap_matrix = merged_matrix
heatmap_columns = list(heatmap_matrix.columns)
org_color_real = {
    real_name: org_colors.get(organ_key, "lightgrey")
    for organ_key, real_name in real_names.items()
}
org_colors_plot = []
type_colors_plot = []
for column in heatmap_columns:
    if " " in column:
        _, real_name = column.split(" ", 1)
        org_colors_plot.append(org_color_real.get(real_name, "lightgrey"))
        type_colors_plot.append("grey" if column.startswith("spatial") else "black")
    else:
        org_colors_plot.append("white")
        type_colors_plot.append("white")
col_colors_df = pd.DataFrame(
    [org_colors_plot, type_colors_plot],
    index=["organ", "type"],
    columns=heatmap_columns,
)

covered = heatmap_matrix.dropna(thresh=7) if not heatmap_matrix.empty else pd.DataFrame()
if covered.empty:
    top_genes: List[str] = []
else:
    sorted_means = covered.mean(axis=1).sort_values()
    if len(covered) >= 40:
        top_genes = sorted_means.iloc[-40:].index.tolist()
    else:
        top_genes = sorted_means.index.tolist()



top_genes_output_path.parent.mkdir(parents=True, exist_ok=True)
with top_genes_output_path.open("w") as handle:
    handle.write("\n".join(top_genes))


_top_interactors: Dict[str, List[str]] = {}
for organ in organs:
    organ_real = real_names[organ]
    columns = [f"spatial {organ_real}", f"expression {organ_real}"]
    if heatmap_matrix.empty or not set(columns).issubset(heatmap_matrix.columns):
        _top_interactors[organ] = []
        continue
    sub = heatmap_matrix[columns].copy()
    sub["mean"] = sub.mean(axis=1)
    sub = sub.dropna()
    if sub.empty:
        _top_interactors[organ] = []
        continue
    top_rows = sub.sort_values(by="mean").tail(10)
    _top_interactors[organ] = list(top_rows.index)

top_interactors = _top_interactors

top_interactors_output_path.parent.mkdir(parents=True, exist_ok=True)
with top_interactors_output_path.open("w") as handle:
    json.dump(top_interactors, handle, indent=2)

significance_dict = (
    test_significance_with_fdr(cosine_sim_full_dict, organs, top_genes)
    if top_genes
    else {}
)

top_genes_heatmap_path.parent.mkdir(parents=True, exist_ok=True)
top_interactors_plot_path.parent.mkdir(parents=True, exist_ok=True)
top_genes_boxplot_path.parent.mkdir(parents=True, exist_ok=True)
heatmap_export_df = pd.DataFrame()
top_interactors_plot_frames: List[pd.DataFrame] = []
top_genes_boxplot_frames: List[pd.DataFrame] = []

if top_genes:
    genes_for_plot = [gene for gene in top_genes if gene in heatmap_matrix.index]
    if genes_for_plot:
        plot_order = list(reversed(genes_for_plot))
        plot_data = heatmap_matrix.loc[plot_order]
        if not plot_data.empty:
            heatmap_export_df = plot_data.copy()
            g = sns.clustermap(
                plot_data,
                figsize = (7,20),
                vmin=0,
                vmax=1,
                cmap="plasma",
                row_cluster=False,
                col_cluster=False,
                cbar_pos=(0, .45, .03, .1), 
                cbar_kws={"label": "percentile"},
                col_colors=[org_colors_plot,type_colors_plot],
                yticklabels=True
            )
            g.fig.savefig(top_genes_heatmap_path, format="pdf", bbox_inches="tight")
            plt.close(g.fig)
        else:
            with PdfPages(top_genes_heatmap_path) as pdf:
                fig, ax = plt.subplots(figsize=(6, 4))
                ax.axis("off")
                ax.text(0.5, 0.5, "No top genes to display", ha="center", va="center")
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
    else:
        with PdfPages(top_genes_heatmap_path) as pdf:
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.axis("off")
            ax.text(0.5, 0.5, "No top genes to display", ha="center", va="center")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
else:
    with PdfPages(top_genes_heatmap_path) as pdf:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.axis("off")
        ax.text(0.5, 0.5, "No top genes to display", ha="center", va="center")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

if heatmap_export_df.empty:
    pd.DataFrame().to_csv(top_genes_heatmap_csv_path, index=False)
else:
    if heatmap_export_df.index.name is None:
        heatmap_export_df.index.name = "interaction"
    heatmap_export_df.to_csv(top_genes_heatmap_csv_path)

with PdfPages(top_interactors_plot_path) as pdf:
    for organ in organs:
        organ_real = real_names[organ]
        genes = [
            gene for gene in top_interactors.get(organ, []) if gene in heatmap_matrix.index
        ]
        if not genes:
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.axis("off")
            ax.text(
                0.5,
                0.5,
                f"No top interactors for {organ_real}",
                ha="center",
                va="center",
            )
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            continue
        plot_order = list(reversed(genes))
        plot_data = heatmap_matrix.loc[plot_order]
        if plot_data.empty:
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.axis("off")
            ax.text(
                0.5,
                0.5,
                f"No data available for {organ_real}",
                ha="center",
                va="center",
            )
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            continue
        top_interactors_plot_frames.append(
            plot_data.reset_index()
            .rename(columns={"index": "interaction"})
            .assign(organ=organ, organ_real=organ_real)
        )
        g = sns.clustermap(
            plot_data,
            figsize = (8,8),
            vmin=0,
            vmax=1,
            cmap="plasma",
            row_cluster=False,
            col_cluster=False,
            cbar_pos=(0, .45, .03, .1), 
            cbar_kws={"label": "percentile"},
            col_colors=[org_colors_plot,type_colors_plot],
        )
        g.fig.suptitle(organ_real)
        pdf.savefig(g.fig, bbox_inches="tight")
        plt.close(g.fig)

if top_interactors_plot_frames:
    top_interactors_export = pd.concat(top_interactors_plot_frames, ignore_index=True)
else:
    top_interactors_export = pd.DataFrame()
top_interactors_export.to_csv(top_interactors_plot_csv_path, index=False)

if top_genes:
    n_organs = len(organs)
    if n_organs == 0:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.axis("off")
        ax.text(0.5, 0.5, "No organs configured", ha="center", va="center")
        fig.savefig(top_genes_boxplot_path, format="pdf", bbox_inches="tight")
        plt.close(fig)
    else:
        ncols = 2 if n_organs > 1 else 1
        nrows = ceil(n_organs / ncols)
        fig, axs = plt.subplots(
            nrows,
            ncols,
            figsize=(32, 10),
            tight_layout=True,
        )
        for idx, organ in enumerate(organs):
            ax = axs[idx // ncols][idx % ncols]
            data = cosine_sim_full_dict.get(organ, pd.DataFrame())
            data = data[data["interaction"].isin(top_genes)].sort_values(by="interaction")
            if data.empty:
                ax.axis("off")
                ax.text(0.5, 0.5, f"No data for {real_names[organ]}", ha="center", va="center")
                continue
            palette = {"control": "white", "fibrosis": org_colors.get(organ, "grey")}
            export_df = data.copy()
            export_df["organ"] = organ
            export_df["organ_real"] = real_names[organ]
            top_genes_boxplot_frames.append(export_df)
            sns.boxplot(
                data=data,
                x="interaction",
                y="mean",
                hue="cond_test",
                hue_order=["control", "fibrosis"],
                width=0.8,
                ax=ax,
                palette=palette,
            )
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
            ax.set_ylabel("mean cosine similarity")
            ax.set_xlabel("")
            ax.set_title(real_names[organ])
            ax.legend(loc="upper left", bbox_to_anchor=(1, 1))

            organ_stats = significance_dict.get(organ, {}) if significance_dict else {}
            interactions = data["interaction"].unique()
            for i, interaction in enumerate(interactions):
                stats = organ_stats.get(interaction, {})
                pval = stats.get("corrected_p") if stats else None
                if pval is not None and pval <= 0.05:
                    group_data = data[data["interaction"] == interaction]
                    group_max = group_data["mean"].max()
                    if pd.isna(group_max):
                        continue
                    y = group_max + 0.02
                    h = 0.01
                    x1 = i - 0.2
                    x2 = i + 0.2
                    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1, c="black")
                    ax.text(i, y + h + 0.005, "*", ha="center", va="bottom", fontsize=14)

        for idx in range(len(organs), nrows * ncols):
            axs[idx // ncols][idx % ncols].axis("off")

        fig.savefig(top_genes_boxplot_path, format="pdf", bbox_inches="tight")
        plt.close(fig)
else:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.axis("off")
    ax.text(0.5, 0.5, "No top genes selected", ha="center", va="center")
    fig.savefig(top_genes_boxplot_path, format="pdf", bbox_inches="tight")
    plt.close(fig)

if top_genes_boxplot_frames:
    top_genes_boxplot_export = pd.concat(top_genes_boxplot_frames, ignore_index=True)
else:
    top_genes_boxplot_export = pd.DataFrame()
top_genes_boxplot_export.to_csv(top_genes_boxplot_csv_path, index=False)
