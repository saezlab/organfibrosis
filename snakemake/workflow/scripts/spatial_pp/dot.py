"""DOT deconvolution and organ-wide niche clustering from dot.ipynb."""

from pathlib import Path
import os

# Use OpenMP to avoid a Numba/TBB shutdown hang.
# Set before importing Scanpy/DOT, which can initialize Numba.
os.environ["NUMBA_THREADING_LAYER"] = snakemake.params.numba_threading_layer
print(
    f"DOT: requested Numba threading layer: {snakemake.params.numba_threading_layer}",
    flush=True,
)

import matplotlib

# File-only plotting for the cluster, before Scanpy/DOT import pyplot.
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from dotpy import DOT, setup_reference, setup_spatial
from dotpy.visualization import plot_spatial_weights
from sklearn.cluster import KMeans


device = snakemake.params.device
if device == "cuda" and not torch.cuda.is_available():
    raise RuntimeError("DOT requires CUDA; request a GPU or set params.device to 'cpu'.")

cond_rename = snakemake.params.condition_map
for output_dir in snakemake.output:
    Path(output_dir).mkdir(parents=True, exist_ok=True)


def run_dot(spatial_processed, ref_processed, sample, organ, spatial_adata):
    # Run DOT with batching
    dot = DOT(
        spatial_processed,
        ref_processed,
        batch_size=snakemake.params.batch_size,
        device=device
    )

    dot.fit(
        mode=snakemake.params.mode,
        max_spot_size=snakemake.params.max_spot_size,
        ratios_weight=snakemake.params.ratios_weight,
        iterations=snakemake.params.iterations,
        checkpoint_dir=str(
            Path(snakemake.output[f'{organ}_checkpoints']) / str(sample)
        ),
        checkpoint_freq=snakemake.params.checkpoint_freq,
        verbose=snakemake.params.verbose
    )

    # Get results
    weights = dot.get_weights(normalize=True)
    cell_types = list(dot.get_cell_types())

    # Visualize results
    plot_spatial_weights(
        spatial_adata.obsm[snakemake.params.spatial_key],
        weights,
        cell_types=cell_types,
        ncols=snakemake.params.plot_ncols,
        save_path=str(
            Path(snakemake.output[f'{organ}_maps']) / f'{sample}_cell_type_maps.png'
        ),
    )

    # Release figures after the map has been saved.
    plt.close("all")

    weights_df = pd.DataFrame(weights, columns=cell_types, index=spatial_adata.obs_names)
    weights_df.to_csv(Path(snakemake.output[f'{organ}_weights']) / f'{sample}_weights.csv')

    # attach the deconvolution weights to the spatial object right after computing them
    spatial_adata.obs[cell_types] = weights_df.values
    if 'counts' not in spatial_adata.layers:
        spatial_adata.layers['counts'] = spatial_adata.X.copy()
    spatial_adata.uns['deconv_celltypes'] = cell_types
    return spatial_adata


KEEP_OBS = snakemake.params.keep_obs


def compute_niches(weights, ct_cols, n_niches, seed):
    """Hellinger (sqrt) transform of the cell-type proportions, then k-means into
    `n_niches` niches. Returns the niche label per row (aligned to `weights`)."""
    X = np.sqrt(weights[ct_cols].to_numpy(float))
    return KMeans(
        n_niches, n_init=snakemake.params.niche_n_init, random_state=seed
    ).fit_predict(X).astype(str)


def _sanitize_keys(adata):
    """Replace '/' in all column/key names h5ad uses as HDF5 keys."""
    adata.obs.columns = [str(c).replace('/', '_') for c in adata.obs.columns]
    adata.var.columns = [str(c).replace('/', '_') for c in adata.var.columns]
    for attr in ('obsm', 'varm'):
        m = getattr(adata, attr)
        for key, val in m.items():
            if isinstance(val, pd.DataFrame):
                val.columns = [str(c).replace('/', '_') for c in val.columns]
    return adata


def _slim(adata, sample):
    """Raw counts in X, only the columns of interest in obs, nothing in var."""
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts'].copy()
    adata.layers.clear()

    if 'sample_uid' in adata.obs.columns:
        adata.obs = adata.obs.rename(columns={'sample_uid': 'sample'})
    if 'sample' not in adata.obs.columns:
        adata.obs['sample'] = str(sample)

    # sanitize before subsetting so a cell type like 'a/b' is kept as 'a_b'
    adata = _sanitize_keys(adata)
    keep = [c.replace('/', '_') for c in KEEP_OBS]
    adata.obs = adata.obs[[c for c in keep if c in adata.obs.columns]].copy()
    adata.var = pd.DataFrame(index=adata.var_names)
    return adata


def save_organ_niches(organ, adata_dict):
    """Pool the per-sample DOT weights of one organ, Hellinger-transform, k-means
    into organ-wide niches, attach the niche label to each sample's AnnData and
    write one slim spatial h5ad (raw counts in X) per sample. These files are read
    by heart/lung/kidney/liver.py."""
    samples = list(adata_dict)
    ct_cols = list(adata_dict[samples[0]].uns['deconv_celltypes'])

    # organ-wide niches on the pooled cell-type proportions
    pooled = pd.concat([adata_dict[s].obs[ct_cols + ['cond_test']] for s in samples], ignore_index=True)
    niche_labels = compute_niches(
        pooled, ct_cols, n_niches=snakemake.params.n_niches,
        seed=snakemake.params.niche_seed,
    )
    niche_names = sorted(pd.unique(niche_labels), key=lambda x: int(x))

    outdir = Path(snakemake.output[f'{organ}_h5ads'])

    start = 0
    for s in samples:
        adata = adata_dict[s]
        n = adata.n_obs
        niche_s = niche_labels[start:start + n]
        start += n

        adata.obs['niche'] = pd.Categorical([str(x) for x in niche_s], categories=niche_names)
        adata = _slim(adata, s)
        adata.uns['niche_names'] = list(niche_names)

        adata = adata.to_memory() if adata.isbacked else adata
        adata.write(outdir / f'{s}.h5ad')

    return pooled


# Lung

# Load data
ref_adata = sc.read_h5ad(snakemake.input.lung_reference)
ref_adata.X = ref_adata.layers['counts'].copy()
visium_path = snakemake.input.lung_spatial

# Process reference and spatial data
ref_processed = setup_reference(
    ref_adata,
    cell_type_key=snakemake.params.cell_type_key,
    subcluster_size=snakemake.params.subcluster_size,
    max_genes=snakemake.params.max_genes,
    verbose=snakemake.params.verbose,
    random_state=snakemake.params.reference_seed,
    remove_mt=snakemake.params.remove_mt,
)

spatial_adata = sc.read_h5ad(visium_path)
spatial_adata.X = spatial_adata.layers['counts'].copy()
dataset_names = list(spatial_adata.uns['spatial'].keys())

spatial_adata.obs["cond_test"] = spatial_adata.obs["treatment"].replace(cond_rename)


adata_dict = {}
for dataset_name in dataset_names:
    spatial_adata_sub = spatial_adata[spatial_adata.obs['sample'] == dataset_name].copy()

    spatial_processed = setup_spatial(
        spatial_adata_sub,
        spatial_key=snakemake.params.spatial_key,
        th_spatial=snakemake.params.th_spatial,
        remove_mt=snakemake.params.remove_mt,
        verbose=snakemake.params.verbose
    )

    adata_dict[dataset_name] = run_dot(spatial_processed, ref_processed, dataset_name, 'lung', spatial_adata_sub)

pooled = save_organ_niches('lung', adata_dict)


# Heart

# Load data
ref_adata = sc.read_h5ad(snakemake.input.heart_reference)
ref_adata.X = ref_adata.layers['counts'].copy()

# Process reference and spatial data
ref_processed = setup_reference(
    ref_adata,
    cell_type_key=snakemake.params.cell_type_key,
    subcluster_size=snakemake.params.subcluster_size,
    max_genes=snakemake.params.max_genes,
    verbose=snakemake.params.verbose,
    random_state=snakemake.params.reference_seed,
    remove_mt=snakemake.params.remove_mt,
)

metadata = pd.read_csv(snakemake.input.heart_metadata)

metadata_sub = metadata[metadata['major_labl'].isin(snakemake.params.heart_regions)]

adata_dict = {}
for count, slide in enumerate(metadata_sub['slide_name']):
    spatial_adata = sc.read(Path(snakemake.input.heart_slides) / f'{slide}.h5ad')
    region = cond_rename[metadata['major_labl'][metadata['slide_name'] == slide].item()]

    spatial_adata.obs.loc[:, 'cond_test'] = region

    # raw counts already in X

    spatial_processed = setup_spatial(
        spatial_adata,
        spatial_key=snakemake.params.spatial_key,
        th_spatial=snakemake.params.th_spatial,
        remove_mt=snakemake.params.remove_mt,
        verbose=snakemake.params.verbose
    )

    adata_dict[slide] = run_dot(spatial_processed, ref_processed, slide, 'heart', spatial_adata)

pooled = save_organ_niches('heart', adata_dict)


# Kidney

# Load data
ref_adata = sc.read_h5ad(snakemake.input.kidney_reference)
ref_adata.X = ref_adata.layers['counts'].copy()

# Process reference and spatial data
ref_processed = setup_reference(
    ref_adata,
    cell_type_key=snakemake.params.cell_type_key,
    subcluster_size=snakemake.params.subcluster_size,
    max_genes=snakemake.params.max_genes,
    verbose=snakemake.params.verbose,
    random_state=snakemake.params.reference_seed,
    remove_mt=snakemake.params.remove_mt,
)

spatial_adata = sc.read_h5ad(snakemake.input.kidney_spatial)

spatial_adata = spatial_adata[
    spatial_adata.obs['condition_full'].isin(snakemake.params.kidney_conditions)
]

spatial_adata.obs['cond_test'] = 'fibrosis'
spatial_adata.obs.loc[
    spatial_adata.obs['condition_full'].isin(snakemake.params.kidney_control_conditions),
    'cond_test',
] = 'control'

exclude_datasets = snakemake.params.kidney_exclude_samples

dataset_names = list(spatial_adata.obs['uid_slide'].unique())

dataset_names_filtered = list(set(dataset_names) - set(exclude_datasets))

adata_dict = {}
for dataset_name in dataset_names_filtered:
    spatial_adata_sub = spatial_adata[spatial_adata.obs['uid_slide'] == dataset_name].copy()

    # counts already in X

    spatial_processed = setup_spatial(
        spatial_adata_sub,
        spatial_key=snakemake.params.spatial_key,
        th_spatial=snakemake.params.th_spatial,
        remove_mt=snakemake.params.remove_mt,
        verbose=snakemake.params.verbose
    )

    adata_dict[dataset_name] = run_dot(spatial_processed, ref_processed, dataset_name, 'kidney', spatial_adata_sub)

pooled = save_organ_niches('kidney', adata_dict)


# Liver

# Load data
ref_adata = sc.read_h5ad(snakemake.input.liver_reference)
ref_adata.X = ref_adata.layers['counts'].copy()

# Process reference and spatial data
ref_processed = setup_reference(
    ref_adata,
    cell_type_key=snakemake.params.cell_type_key,
    subcluster_size=snakemake.params.subcluster_size,
    max_genes=snakemake.params.max_genes,
    verbose=snakemake.params.verbose,
    random_state=snakemake.params.reference_seed,
    remove_mt=snakemake.params.remove_mt,
)

adata_dict = {}
for sample, spatial_path in zip(snakemake.params.liver_samples, snakemake.input.liver_spatial):
    spatial_adata = sc.read_h5ad(spatial_path)
    spatial_adata.X = spatial_adata.layers['counts'].copy()
    spatial_adata.obs["cond_test"] = spatial_adata.obs["disease"].replace(cond_rename)

    spatial_processed = setup_spatial(
        spatial_adata,
        spatial_key=snakemake.params.spatial_key,
        th_spatial=snakemake.params.th_spatial,
        remove_mt=snakemake.params.remove_mt,
        verbose=snakemake.params.verbose
    )

    adata_dict[sample] = run_dot(spatial_processed, ref_processed, sample, 'liver', spatial_adata)

pooled = save_organ_niches('liver', adata_dict)

print("DOT: all outputs written.", flush=True)
