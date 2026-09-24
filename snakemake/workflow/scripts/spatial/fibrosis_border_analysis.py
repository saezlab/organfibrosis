from pathlib import Path
import gc

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from matplotlib.lines import Line2D
from sklearn.neighbors import KDTree, NearestNeighbors


plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update({"font.size": 15})
sns.set_style("whitegrid")

organs = snakemake.params["organs"]
organ_names = snakemake.params["organ_names"]
celltype_colors = dict(snakemake.params["celltype_colors"])
fibro_niches = snakemake.params["fibro_niches"]
refine_params = snakemake.params["refine"]
example_slides = snakemake.params["example_slides"]
dot_sizes = snakemake.params["dot_sizes"]

profile_microns_out = dict(
    zip(organ_names.values(), snakemake.output["profiles_microns"])
)
profile_normalized_out = dict(
    zip(organ_names.values(), snakemake.output["profiles_normalized"])
)
example_slides_out = dict(
    zip(organ_names.values(), snakemake.output["example_slides"])
)

deconv_h5ads = {}
deconv_paths = []
for directory in snakemake.input["deconv_dirs"]:
    sample_paths = sorted(Path(directory).glob("*.h5ad"))
    if not sample_paths:
        raise ValueError(f"No DOT H5AD files found in {directory}")
    deconv_paths.extend(sample_paths)
for path in deconv_paths:
    organ = Path(path).parents[1].name
    deconv_h5ads.setdefault(organ, []).append(path)

ulm_estimates = {
    Path(path).name.removesuffix("_ulm_estimate.csv"): path
    for path in snakemake.input["ulm_estimates"]
}

keep_obs = [
    "niche",
    "cond_test",
    "sample",
    "endothelial",
    "lymphoid",
    "mesenchymal",
    "myeloid",
    "epithelial",
    "other",
]
celltypes = [
    "endothelial",
    "lymphoid",
    "mesenchymal",
    "myeloid",
    "epithelial",
    "other",
]
celltype_colors["other"] = "grey"


def refine_region_mask(
    adata,
    region_key,
    fibrosis_label,
    sample_key="sample",
    spatial_key="spatial",
    key_added="fib_refined",
):
    is_fibrotic = np.asarray(adata.obs[region_key].values) == fibrosis_label
    refined = is_fibrotic.copy()

    for sample in adata.obs[sample_key].unique():
        sample_mask = (adata.obs[sample_key] == sample).values
        coordinates = adata.obsm[spatial_key][sample_mask]
        n_neighbors = min(refine_params["n_neighs"] + 1, len(coordinates))
        _, neighbors = NearestNeighbors(n_neighbors=n_neighbors).fit(
            coordinates
        ).kneighbors(coordinates)
        neighbors = neighbors[:, 1:]
        current = is_fibrotic[sample_mask].copy()

        for _ in range(int(refine_params["smooth_iters"])):
            current = (
                current[neighbors].mean(axis=1)
                >= refine_params["smooth_thresh"]
            )
        if refine_params["fill_holes"]:
            current = current | (
                (~current) & (current[neighbors].mean(axis=1) == 1.0)
            )
        for _ in range(int(refine_params["expand"])):
            current = current | ((~current) & current[neighbors].any(axis=1))

        refined[sample_mask] = current

    adata.obs[key_added] = refined


def signed_border_distance(
    adata,
    region_key,
    fibrosis_label,
    sample_key="sample",
    spatial_key="spatial",
):
    raw = np.full(adata.n_obs, np.nan)

    for sample in adata.obs[sample_key].unique():
        sample_mask = (adata.obs[sample_key] == sample).values
        coordinates = adata.obsm[spatial_key][sample_mask]
        is_fibrotic = (
            np.asarray(adata.obs[region_key].values[sample_mask])
            == fibrosis_label
        )
        if is_fibrotic.all() or (~is_fibrotic).all():
            continue

        fibrotic_coordinates = coordinates[is_fibrotic]
        other_coordinates = coordinates[~is_fibrotic]
        distance = np.empty(sample_mask.sum())
        distance[is_fibrotic] = -KDTree(other_coordinates).query(
            fibrotic_coordinates
        )[0].ravel()
        distance[~is_fibrotic] = KDTree(fibrotic_coordinates).query(
            other_coordinates
        )[0].ravel()
        raw[sample_mask] = distance

    adata.obs["border_signed_raw"] = raw
    normalized = np.full(adata.n_obs, np.nan)

    for sample in adata.obs[sample_key].unique():
        sample_mask = (adata.obs[sample_key] == sample).values
        values = raw[sample_mask]
        output = np.full_like(values, np.nan)
        finite = ~np.isnan(values)
        output[finite] = 0.0
        negative = finite & (values < 0)
        positive = finite & (values > 0)
        if negative.any():
            output[negative] = values[negative] / (-values[negative].min())
        if positive.any():
            output[positive] = values[positive] / values[positive].max()
        normalized[sample_mask] = output

    adata.obs["border_signed"] = normalized


def add_signed_microns(
    adata,
    sample_key="sample",
    spatial_key="spatial",
    spot_pitch_um=100.0,
):
    micron_distance = np.full(adata.n_obs, np.nan)

    for sample in adata.obs[sample_key].unique():
        sample_mask = (adata.obs[sample_key] == sample).values
        coordinates = adata.obsm[spatial_key][sample_mask]
        if len(coordinates) < 2:
            continue
        distances, _ = NearestNeighbors(n_neighbors=2).fit(
            coordinates
        ).kneighbors(coordinates)
        microns_per_unit = spot_pitch_um / float(np.median(distances[:, 1]))
        micron_distance[sample_mask] = (
            adata.obs.loc[sample_mask, "border_signed_raw"].to_numpy()
            * microns_per_unit
        )

    adata.obs["border_signed_um"] = micron_distance


def supported_distance_range(distance, trim_percent=5.0):
    distance = distance[~np.isnan(distance)]
    if len(distance) == 0 or np.ptp(distance) == 0:
        return np.nan, np.nan
    negative = distance[distance < 0]
    positive = distance[distance > 0]
    lower = (
        float(np.percentile(negative, trim_percent))
        if len(negative)
        else float(distance.min())
    )
    upper = (
        float(np.percentile(positive, 100 - trim_percent))
        if len(positive)
        else float(distance.max())
    )
    return lower, upper


def core_tissue_header(ax, x_offset):
    transform = ax.get_xaxis_transform()
    ax.text(
        -x_offset,
        1.02,
        "← core",
        transform=transform,
        fontsize=12,
        color="0.4",
        ha="right",
        va="bottom",
        clip_on=False,
    )
    ax.text(
        x_offset,
        1.02,
        "tissue →",
        transform=transform,
        fontsize=12,
        color="0.4",
        ha="left",
        va="bottom",
        clip_on=False,
    )


def axis_settings(normalized):
    if normalized:
        return 0.067, (-1.05, 1.05), (
            "normalized signed distance (−1 core → 0 border → +1 tissue)"
        ), 0.03, 0.05
    return 130, (-1000, 1500), (
        "signed distance to fibrosis border (µm)"
    ), 25, 25


def binned_profile(distance, lower, upper, bin_width):
    edges = np.arange(lower, upper + bin_width, bin_width)
    centers = 0.5 * (edges[:-1] + edges[1:])
    indices = np.clip(
        np.digitize(distance, edges[1:-1]), 0, len(centers) - 1
    )
    return centers, indices


def plot_celltypes_along_axis(
    adata,
    dist_key,
    normalized,
    ax,
    min_per_bin=10,
):
    distance = np.asarray(adata.obs[dist_key].values, dtype=float)
    bin_width, x_limits, x_label, header_offset, end_padding = axis_settings(
        normalized
    )
    lower, upper = supported_distance_range(distance)
    centers, indices = binned_profile(distance, lower, upper, bin_width)

    profiles = {celltype: np.full(len(centers), np.nan) for celltype in celltypes}
    for index in range(len(centers)):
        selected = (
            (indices == index) & (distance >= lower) & (distance <= upper)
        )
        if selected.sum() >= min_per_bin:
            for celltype in celltypes:
                profiles[celltype][index] = adata.obs[celltype].values[
                    selected
                ].mean()

    order = (
        adata.obs[celltypes].mean().sort_values(ascending=False).index.tolist()
    )
    for celltype in order:
        finite = np.isfinite(profiles[celltype])
        ax.plot(
            centers[finite],
            profiles[celltype][finite],
            color=celltype_colors[celltype],
            linewidth=2.0,
            marker="o",
            markersize=3,
        )
        if finite.any():
            ax.text(
                x_limits[1] + end_padding,
                profiles[celltype][finite][-1],
                celltype,
                color=celltype_colors[celltype],
                fontsize=12,
                va="center",
            )

    ax.axvline(0, color="0.35", linestyle="--", linewidth=1.0)
    ax.set_xlabel(x_label)
    ax.set_ylabel(f"mean cell type \n prop. [%]")
    ax.set_xlim(x_limits)
    core_tissue_header(ax, header_offset)


def plot_score_along_axis(
    adata,
    score_key,
    dist_key,
    normalized,
    color,
    ax,
    min_per_bin=10,
):
    distance = np.asarray(adata.obs[dist_key].values, dtype=float)
    score = np.asarray(adata.obs[score_key].values, dtype=float)
    bin_width, x_limits, x_label, header_offset, _ = axis_settings(normalized)
    lower, upper = supported_distance_range(distance)
    centers, indices = binned_profile(distance, lower, upper, bin_width)
    mean = np.full(len(centers), np.nan)

    for index in range(len(centers)):
        selected = (
            (indices == index)
            & (distance >= lower)
            & (distance <= upper)
            & ~np.isnan(score)
        )
        if selected.sum() >= min_per_bin:
            mean[index] = score[selected].mean()

    finite = np.isfinite(mean)
    ax.plot(
        centers[finite],
        mean[finite],
        color=color,
        linewidth=2.0,
        marker="o",
        markersize=3,
    )
    ax.axvline(0, color="0.35", linestyle="--", linewidth=1.0)
    ax.set_xlim(x_limits)
    ax.set_xlabel(x_label)
    ax.margins(x=0.03)
    core_tissue_header(ax, header_offset)


def plot_example_slide_panel(adata, organ, sample, dot_size):
    slide_mask = (adata.obs["organ"] == organ) & (
        adata.obs["sample"].astype(str) == sample
    )
    if slide_mask.sum() == 0:
        raise ValueError(f"Sample {sample} not found for organ {organ}")

    subset = adata[slide_mask].copy()
    subset.obs["niche_grouping"] = np.where(
        subset.obs["fibrosis"].astype(str) == "fibrosis",
        "fibrosis",
        "not_fibrosis",
    )
    refine_region_mask(
        subset,
        region_key="niche_grouping",
        fibrosis_label="fibrosis",
    )
    signed_border_distance(
        subset,
        region_key="fib_refined",
        fibrosis_label=True,
    )
    add_signed_microns(subset)

    coordinates = np.asarray(subset.obsm["spatial"])
    x, y = coordinates[:, 0], coordinates[:, 1]
    figure, axes = plt.subplots(5, 1, figsize=(4.2, 18))

    score = subset.obs["NABA_CORE_MATRISOME"].astype(float)
    valid_score = np.isfinite(score)
    vmin, vmax = (
        np.nanpercentile(score[valid_score], [5, 95])
        if valid_score.any()
        else (0, 1)
    )
    score_scatter = axes[0].scatter(
        x[valid_score],
        y[valid_score],
        c=score[valid_score],
        s=dot_size,
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
        linewidths=0,
        rasterized=True,
    )
    #axes[0].set_title("ECM score")
    score_colorbar = figure.colorbar(
        score_scatter, ax=axes[0], fraction=0.045, pad=0.03
    )
    score_colorbar.set_label("enrichment score", fontsize=14)

    niches = subset.obs["niche"].astype(str)
    unique_niches = pd.Index(niches.dropna().unique())
    colors = mpl.colormaps["Set2"].colors[:5]
    color_map = {
        value: colors[index % len(colors)]
        for index, value in enumerate(unique_niches)
    }
    axes[1].scatter(
        x,
        y,
        c=[color_map[value] for value in niches],
        s=dot_size,
        linewidths=0,
        rasterized=True,
    )
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=5,
            markerfacecolor=color_map[value],
            markeredgecolor="none",
            label=str(value),
        )
        for value in sorted(
            unique_niches,
            key=lambda value: (
                (0, int(value), "")
                if str(value).isdigit()
                else (1, 0, str(value))
            ),
        )
    ]
    axes[1].legend(
        handles=handles,
        title="Niche",
        loc="upper center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=min(len(handles), 5),
        frameon=False,
        fontsize=14,
        title_fontsize=14,
        handletextpad=0.2,
        columnspacing=0.7,
        borderpad=0.1,
    )
    #axes[1].set_title("Niches")

    binary = (
        subset.obs["niche_grouping"].astype(str) == "fibrosis"
    ).astype(int)
    binary_cmap = mpl.colors.ListedColormap(["#d8c7a0", "#1b5e20"])
    axes[2].scatter(
        x,
        y,
        c=binary,
        s=dot_size,
        cmap=binary_cmap,
        vmin=0,
        vmax=1,
        linewidths=0,
        rasterized=True,
    )
    #axes[2].set_title("Fibrosis niche")

    axes[3].scatter(
        x,
        y,
        c=subset.obs["fib_refined"].astype(bool).astype(int),
        s=dot_size,
        cmap=binary_cmap,
        vmin=0,
        vmax=1,
        linewidths=0,
        rasterized=True,
    )
    #axes[3].set_title("Cleaned fibrosis niche")

    distance = subset.obs["border_signed_um"].astype(float)
    valid_distance = np.isfinite(distance)
    max_distance = (
        np.nanmax(np.abs(distance[valid_distance]))
        if valid_distance.any()
        else 1.0
    )
    distance_scatter = axes[4].scatter(
        x[valid_distance],
        y[valid_distance],
        c=distance[valid_distance],
        s=dot_size,
        cmap="RdBu",
        vmin=-max_distance,
        vmax=max_distance,
        linewidths=0,
        rasterized=True,
    )
    #axes[4].set_title("Signed border distance")
    distance_colorbar = figure.colorbar(
        distance_scatter, ax=axes[4], fraction=0.045, pad=0.03
    )
    distance_colorbar.set_label("distance to border [µm]", fontsize=14)

    for axis in axes:
        axis.set_aspect("equal")
        axis.axis("off")
        axis.invert_yaxis()

    figure.tight_layout()
    return figure


adata_parts = []
for organ_key in organs:
    organ = organ_names[organ_key]
    for h5ad_path in deconv_h5ads[organ]:
        adata = sc.read_h5ad(h5ad_path, backed="r")
        sample = Path(h5ad_path).stem
        counts = adata.X[:]
        if not sp.issparse(counts):
            counts = sp.csr_matrix(counts)
        counts = counts.astype("float32")

        count_values = counts[:50].data
        if count_values.size and not np.allclose(
            count_values, np.round(count_values)
        ):
            raise ValueError(f"{sample}: selected matrix is not integer counts")

        obs = adata.obs[
            [column for column in keep_obs if column in adata.obs.columns]
        ].copy()
        obs["sample"] = sample
        obs["organ"] = organ
        obs.index = [f"{sample}_{barcode}" for barcode in obs.index]
        adata_parts.append(
            ad.AnnData(
                X=counts,
                obs=obs,
                var=pd.DataFrame(index=adata.var_names),
                obsm={"spatial": adata.obsm["spatial"]},
            )
        )
        del adata, counts
        gc.collect()

combined = ad.concat(adata_parts, join="outer", index_unique=None)
niches = combined.obs["niche"].astype(str)
combined.obs["fibrosis"] = pd.Categorical(
    np.where(
        [
            niche in fibro_niches[organ]
            for niche, organ in zip(niches, combined.obs["organ"])
        ],
        "fibrosis",
        "not_fibrosis",
    ),
    categories=["not_fibrosis", "fibrosis"],
    ordered=True,
)

id_columns = {"organ", "sample", "cell", "condition", "spot_id"}
ulm_parts = []
for organ in organ_names.values():
    ulm = pd.read_csv(ulm_estimates[organ])
    if "spot_id" not in ulm.columns:
        ulm["spot_id"] = (
            ulm["sample"].astype(str) + "_" + ulm["cell"].astype(str)
        )
    ulm["organ"] = organ
    ulm_parts.append(ulm.set_index("spot_id"))

ulm = pd.concat(ulm_parts)
if not ulm.index.is_unique:
    raise ValueError("Duplicate spot IDs in ULM estimate tables")
geneset_columns = [
    column for column in ulm.columns if column not in id_columns
]
combined.obs[geneset_columns] = ulm.reindex(combined.obs_names)[geneset_columns]

if not combined.obs_names.is_unique:
    raise ValueError("Duplicate spot IDs in combined spatial data")
sc.pp.normalize_total(combined, target_sum=1e4)
sc.pp.log1p(combined)

for organ in organ_names.values():
    subset = combined[
        (combined.obs["organ"] == organ)
        & (combined.obs["cond_test"] == "fibrosis")
    ].copy()
    subset.obs["niche_grouping"] = subset.obs["fibrosis"]
    refine_region_mask(
        subset,
        region_key="niche_grouping",
        fibrosis_label="fibrosis",
    )
    signed_border_distance(
        subset,
        region_key="fib_refined",
        fibrosis_label=True,
    )
    add_signed_microns(subset)

    for column in (
        "border_signed_raw",
        "border_signed_um",
        "border_signed",
    ):
        combined.obs.loc[subset.obs_names, column] = subset.obs[column].values

    for distance_key, normalized, output_path in (
        ("border_signed_um", False, profile_microns_out[organ]),
        ("border_signed", True, profile_normalized_out[organ]),
    ):
        figure, axes = plt.subplots(4, 1, figsize=(5.5, 12))
        plot_celltypes_along_axis(
            subset,
            dist_key=distance_key,
            normalized=normalized,
            ax=axes[0],
        )
        plot_score_along_axis(
            subset,
            "NABA_CORE_MATRISOME",
            dist_key=distance_key,
            normalized=normalized,
            color="#4EAB85",
            ax=axes[1],
        )
        plot_score_along_axis(
            subset,
            "HALLMARK_INFLAMMATORY_RESPONSE",
            dist_key=distance_key,
            normalized=normalized,
            color="#9A275A",
            ax=axes[2],
        )
        plot_score_along_axis(
            subset,
            "TGFb",
            dist_key=distance_key,
            normalized=normalized,
            color="#487AB8",
            ax=axes[3],
        )
        axes[1].set_ylabel("ECM score")
        axes[2].set_ylabel("inflammation score")
        axes[3].set_ylabel("TGFβ activity score")

        # Widen the liver inflammation range; other organs keep autoscaling.
        if organ == "liver":
            axes[2].set_ylim(-0.5, 1)

        figure.tight_layout()
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_path, bbox_inches="tight")
        plt.close(figure)

for organ in sorted(organ_names.values()):
    figure = plot_example_slide_panel(
        combined,
        organ,
        sample=example_slides[organ],
        dot_size=dot_sizes[organ],
    )
    output_path = example_slides_out[organ]
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close(figure)
