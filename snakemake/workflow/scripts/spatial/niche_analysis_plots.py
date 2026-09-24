from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy import stats
from statsmodels.stats.multitest import multipletests


plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update({"font.size": 18})
sns.set_style("whitegrid")

# Snakemake configuration and output ordering shared with the spatial rule.
organ_names = snakemake.params["organ_names"]
organ_colors = snakemake.params["organ_colors"]
celltype_colors = dict(snakemake.params["celltype_colors"])
sample_columns = snakemake.params["sample_columns"]
condition_names = snakemake.params["condition_names"]
condition_order = list(snakemake.params["condition_order"])
core_matrisome = snakemake.params["core_matrisome"]
celltype_colors["other"] = "grey"

real_organs = list(organ_names.values())
config_for_organ = {real: config for config, real in organ_names.items()}
all_celltypes = [
    "endothelial",
    "epithelial",
    "lymphoid",
    "mesenchymal",
    "myeloid",
    "other",
]

deconvolution_paths = dict(
    zip(real_organs, snakemake.input["deconvolution"])
)
colocalization_paths = dict(
    zip(real_organs, snakemake.input["colocalization"])
)
signature_plot_paths = dict(
    zip(real_organs, snakemake.output["signature_plots"])
)
signature_data_paths = dict(
    zip(real_organs, snakemake.output["signature_data"])
)
morans_plot_paths = dict(zip(real_organs, snakemake.output["morans_plots"]))
morans_data_paths = dict(zip(real_organs, snakemake.output["morans_data"]))


def clr_transform(data):
    """Apply a row-wise centered log-ratio transform to niche fractions.

    Zeros receive half the smallest nonzero fraction in their sample before
    normalization and log transformation.
    """
    values = data.to_numpy(dtype=float)
    transformed = np.empty_like(values)
    for index, row in enumerate(values):
        nonzero = row[row > 0]
        pseudocount = nonzero.min() / 2 if nonzero.size else 1e-6
        adjusted = np.where(row > 0, row, pseudocount)
        adjusted = adjusted / adjusted.sum()
        transformed[index] = np.log(adjusted) - np.log(adjusted).mean()
    return pd.DataFrame(
        transformed,
        index=data.index,
        columns=data.columns,
    )


def significance_stars(pvalue):
    if pvalue < 0.01:
        return "**"
    if pvalue < 0.05:
        return "*"
    return ""


def draw_niche_composition(
    axis,
    organ,
    plot_data,
    pvalues,
    niche_order,
    value_column,
    ylabel,
):
    """Draw one organ panel in a combined niche-composition figure."""
    fibrosis_color = organ_colors[config_for_organ[organ]]
    palette = {
        "control": "lightgray",
        "fibrosis": fibrosis_color,
    }
    data = plot_data.copy()
    data["cond_test"] = pd.Categorical(
        data["cond_test"],
        categories=condition_order,
        ordered=True,
    )

    sns.boxplot(
        data=data,
        x="niche",
        y=value_column,
        order=niche_order,
        hue="cond_test",
        hue_order=condition_order,
        palette=palette,
        showfliers=False,
        linewidth=1,
        saturation=0.9,
        width=0.65,
        ax=axis,
    )
    sns.stripplot(
        data=data,
        x="niche",
        y=value_column,
        order=niche_order,
        hue="cond_test",
        hue_order=condition_order,
        dodge=True,
        color="0.15",
        size=4,
        jitter=0.12,
        linewidth=0.4,
        edgecolor="white",
        legend=False,
        ax=axis,
    )

    adjusted_pvalues = pvalues.set_index("niche")["p_adj"]
    maximum = data[value_column].max()
    minimum = data[value_column].min()
    padding = 0.03 * (maximum - minimum if maximum > minimum else 1)
    for index, niche in enumerate(niche_order):
        marker = significance_stars(adjusted_pvalues.get(str(niche), 1))
        if marker:
            top = data.loc[
                data["niche"].astype(str) == str(niche), value_column
            ].max()
            axis.text(
                index,
                top + padding,
                marker,
                ha="center",
                va="bottom",
                fontsize=14,
                fontweight="bold",
            )

    axis.set_xlabel("")
    axis.set_ylabel(ylabel)
    axis.set_title(organ, loc="left")
    axis.tick_params(axis="x", rotation=45)
    for label in axis.get_xticklabels():
        label.set_ha("right")
    if value_column == "clr":
        axis.axhline(0, color="0.7", linewidth=0.8, zorder=0)
    if axis.get_legend() is not None:
        axis.get_legend().remove()
    sns.despine(ax=axis)


def composition_grid(
    data,
    pvalues,
    niche_orders,
    value_column,
    ylabel,
    output_path,
):
    """Create a four-organ grid for raw fractions or CLR values."""
    figure, axes = plt.subplots(
        4,
        1,
        figsize=(4, 14),
        sharey=True,
    )
    axes = axes.ravel()

    for axis, organ in zip(axes, real_organs):
        draw_niche_composition(
            axis,
            organ,
            data[data["organ"] == organ],
            pvalues[pvalues["organ"] == organ],
            niche_orders[organ],
            value_column,
            ylabel,
        )

    maximum = data[value_column].max()
    minimum = data[value_column].min()
    padding = 0.03 * (maximum - minimum if maximum > minimum else 1)
    axes[0].set_ylim(top=maximum + 6 * padding)

    handles = [
        Patch(
            facecolor="lightgray",
            edgecolor="0.4",
            label="control",
        )
    ]
    handles.extend(
        Patch(
            facecolor=organ_colors[config_for_organ[organ]],
            edgecolor="0.4",
            label=f"{organ} fibrosis",
        )
        for organ in real_organs
    )
    figure.legend(
        handles=handles,
        frameon=False,
        fontsize=16,
        loc="upper left",
        bbox_to_anchor=(1.0, 1.0),
    )
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, bbox_inches="tight")
    plt.close(figure)


# Load spot-level deconvolution tables and spatial colocalization results.
# Niche labels stay numerically ordered so plots do not sort them as strings.
deconvolution = {}
colocalization = {}
for organ in real_organs:
    deconv = pd.read_csv(deconvolution_paths[organ])
    deconv = deconv.rename(
        columns={sample_columns[organ]: "sample"}
    )
    if "sample" not in deconv.columns:
        raise ValueError(
            f"No sample column for {organ}; expected "
            f"{sample_columns[organ]}"
        )
    if "niche" not in deconv.columns:
        raise ValueError(f"No niche column for {organ}")

    niche_strings = deconv["niche"].astype(str)
    niche_order = sorted(niche_strings.unique(), key=int)
    deconv["niche"] = pd.Categorical(
        niche_strings,
        categories=niche_order,
        ordered=True,
    )
    deconvolution[organ] = deconv
    colocalization[organ] = pd.read_csv(colocalization_paths[organ])

# Summarize each organ once and reuse these tables across plots and CSV exports.
results = {}
pvalue_tables = []
fraction_tables = []
clr_tables = []

for organ in real_organs:
    deconv = deconvolution[organ]
    cosine = colocalization[organ].copy()
    celltypes = [
        celltype for celltype in all_celltypes if celltype in deconv.columns
    ]

    sample_conditions = (
        cosine.drop_duplicates("sample")
        .set_index("sample")["cond_test"]
        .replace(condition_names)
    )
    composition = pd.crosstab(
        deconv["sample"],
        deconv["niche"],
        normalize="index",
    )
    composition.columns = [str(column) for column in composition.columns]
    composition["cond_test"] = composition.index.map(sample_conditions)
    composition = composition[
        composition["cond_test"].isin(condition_order)
    ]
    niche_order = composition.columns.drop("cond_test").tolist()

    composition_clr = clr_transform(composition[niche_order])
    composition_clr["cond_test"] = composition["cond_test"].values

    # Long tables retain every sample-level point shown in the box/strip plots.
    fraction_long = composition.melt(
        id_vars="cond_test",
        var_name="niche",
        value_name="fraction",
        ignore_index=False,
    ).reset_index()
    fraction_long.insert(0, "organ", organ)
    clr_long = composition_clr.melt(
        id_vars="cond_test",
        var_name="niche",
        value_name="clr",
        ignore_index=False,
    ).reset_index()
    clr_long.insert(0, "organ", organ)
    fraction_tables.append(fraction_long)
    clr_tables.append(clr_long)

    test_rows = []
    for niche in niche_order:
        control_clr = composition_clr.loc[
            composition_clr["cond_test"] == condition_order[0], niche
        ]
        fibrosis_clr = composition_clr.loc[
            composition_clr["cond_test"] == condition_order[1], niche
        ]
        control_fraction = composition.loc[
            composition["cond_test"] == condition_order[0], niche
        ]
        fibrosis_fraction = composition.loc[
            composition["cond_test"] == condition_order[1], niche
        ]
        statistic, pvalue = stats.mannwhitneyu(
            control_clr,
            fibrosis_clr,
        )
        test_rows.append(
            {
                "organ": organ,
                "niche": str(niche),
                "n_control": len(control_clr),
                "n_fibrosis": len(fibrosis_clr),
                "mean_fraction_control": control_fraction.mean(),
                "mean_fraction_fibrosis": fibrosis_fraction.mean(),
                "mean_clr_control": control_clr.mean(),
                "mean_clr_fibrosis": fibrosis_clr.mean(),
                "mannwhitney_u": statistic,
                "p_value": pvalue,
            }
        )

    pvalues = pd.DataFrame(test_rows)
    # Control false discoveries across the niches tested within this organ.
    pvalues["p_adj"] = multipletests(
        pvalues["p_value"],
        method="fdr_bh",
    )[1]
    pvalue_tables.append(pvalues)

    signature = deconv.groupby("niche", observed=False)[celltypes].mean()
    signature["n_spots"] = deconv["niche"].value_counts()
    signature_z = (
        signature[celltypes] - deconv[celltypes].mean()
    ) / deconv[celltypes].std()

    # The Moran panels use only niche colocalization with the core matrisome.
    cosine["cond_test"] = cosine["cond_test"].replace(condition_names)
    cosine = cosine[cosine["cond_test"].isin(condition_order)]
    cosine_core = cosine[cosine["receptor"] == core_matrisome]

    results[organ] = {
        "niche_order": niche_order,
        "fraction_long": fraction_long,
        "clr_long": clr_long,
        "pvalues": pvalues,
        "celltypes": celltypes,
        "signature": signature,
        "signature_z": signature_z,
        "cosine_core": cosine_core,
    }

fraction_data = pd.concat(fraction_tables, ignore_index=True)
clr_data = pd.concat(clr_tables, ignore_index=True)
pvalue_data = pd.concat(pvalue_tables, ignore_index=True)
niche_orders = {
    organ: results[organ]["niche_order"] for organ in real_organs
}

# Export all sample-level composition values and the complete CLR test table.
Path(snakemake.output["composition_fraction_data"]).parent.mkdir(
    parents=True,
    exist_ok=True,
)
fraction_data.to_csv(
    snakemake.output["composition_fraction_data"],
    index=False,
)
clr_data.to_csv(
    snakemake.output["composition_clr_data"],
    index=False,
)
pvalue_data.to_csv(
    snakemake.output["composition_pvalues"],
    index=False,
)

# Draw raw-fraction and CLR views of the same four-organ comparison.
composition_grid(
    fraction_data,
    pvalue_data,
    niche_orders,
    value_column="fraction",
    ylabel="% of spots",
    output_path=snakemake.output["composition_fraction_plot"],
)
composition_grid(
    clr_data,
    pvalue_data,
    niche_orders,
    value_column="clr",
    ylabel="CLR",
    output_path=snakemake.output["composition_clr_plot"],
)

# Create one signature panel and one niche–ECM Moran panel per organ. Each CSV
# shares its basename with the corresponding PDF.
for organ in real_organs:
    result = results[organ]
    celltypes = result["celltypes"]
    signature = result["signature"]
    signature_z = result["signature_z"]

    signature_data = signature.rename(
        columns={
            celltype: f"fraction_{celltype}" for celltype in celltypes
        }
    ).join(
        signature_z.rename(
            columns={
                celltype: f"zscore_{celltype}" for celltype in celltypes
            }
        )
    )
    signature_data.index = signature_data.index.astype(str)
    signature_data.index.name = "niche"
    signature_data.reset_index().to_csv(
        signature_data_paths[organ],
        index=False,
    )

    # Left: raw mean fractions; right: standardized fractions.
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(12, 6),
        gridspec_kw={"width_ratios": [1, 1.7]},
    )
    sns.heatmap(
        signature_z,
        annot=True,
        fmt=".1f",
        cmap="RdBu_r",
        center=0,
        ax=axes[1],
        vmin=-2,
        vmax=2,
        cbar_kws={"label": "z-scored fraction"},
    )

    bottom = np.zeros(len(signature))
    x_values = [str(niche) for niche in signature.index]
    for celltype in celltypes:
        axes[0].bar(
            x_values,
            signature[celltype].to_numpy(),
            bottom=bottom,
            width=0.8,
            color=celltype_colors.get(celltype),
            label=celltype,
            edgecolor="white",
            linewidth=0.6,
        )
        bottom += signature[celltype].to_numpy()
    axes[0].set_xlabel("niche")
    axes[0].set_ylabel("mean cell type prop. [%]")
    axes[0].legend(
        frameon=False,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=14,
    )
    axes[0].spines[["top", "right"]].set_visible(False)
    figure.suptitle(organ)
    figure.tight_layout()
    figure.savefig(signature_plot_paths[organ], bbox_inches="tight")
    plt.close(figure)

    niche_morans = result["cosine_core"][
        ~result["cosine_core"]["ligand"].isin(all_celltypes)
    ].copy()
    niche_order = sorted(
        niche_morans["ligand"].unique(),
        key=lambda value: int(str(value).rsplit("_", 1)[-1]),
    )
    niche_morans.to_csv(morans_data_paths[organ], index=False)

    # Every row in niche_morans contributes one point to this boxplot.
    figure, axis = plt.subplots(figsize=(4, 4), tight_layout=True)
    sns.boxplot(
        data=niche_morans,
        x="ligand",
        y="morans",
        order=niche_order,
        color=organ_colors[config_for_organ[organ]],
        ax=axis,
    )
    #figure.suptitle(f"{organ}\nMoran's with ECM score")
    # Show only the niche number instead of the full "niche_<n>" label.
    axis.set_xticks(range(len(niche_order)))
    axis.set_xticklabels([str(niche).rsplit("_", 1)[-1] for niche in niche_order])
    axis.set_xlabel("")
    axis.set_ylabel("morans R")
    figure.savefig(morans_plot_paths[organ], bbox_inches="tight")
    plt.close(figure)
