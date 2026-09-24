"""Prepare DOT reference and liver Visium inputs from prepare_liver.ipynb."""

import gzip
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.image import imread
from scipy import sparse
from scipy.io import mmread


for output_file in snakemake.output:
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)


# Liver single-cell reference


# read in anndata
def read_geo_10x(matrix_file):
    matrix_file = Path(matrix_file)

    prefix = str(matrix_file)
    prefix = prefix.replace("_matrix.mtx.gz", "")
    prefix = prefix.replace("matrix.mtx.gz", "")  # handles the odd C51_flush file

    barcode_file = Path(prefix + "_barcodes.tsv.gz")

    if Path(prefix + "_features.tsv.gz").exists():
        feature_file = Path(prefix + "_features.tsv.gz")
    elif Path(prefix + "_genes.tsv.gz").exists():
        feature_file = Path(prefix + "_genes.tsv.gz")
    else:
        raise FileNotFoundError(f"No features/genes file found for {matrix_file}")

    # read matrix
    X = mmread(matrix_file).T.tocsr()

    # read metadata
    barcodes = pd.read_csv(barcode_file, header=None, sep="\t")
    features = pd.read_csv(feature_file, header=None, sep="\t")

    adata = ad.AnnData(X=sparse.csr_matrix(X))

    adata.obs_names = barcodes.iloc[:, 0].astype(str).values

    # works for both genes.tsv and features.tsv
    if features.shape[1] >= 2:
        adata.var_names = features.iloc[:, 1].astype(str).values
        adata.var["gene_ids"] = features.iloc[:, 0].astype(str).values
    else:
        adata.var_names = features.iloc[:, 0].astype(str).values

    adata.var_names_make_unique()

    # sample metadata
    fname = Path(prefix).name
    parts = fname.split("_")

    read_orientation = None
    for part in parts:
        if part in {"3pr", "5pr"}:
            read_orientation = part
            break

    adata.obs["gsm_id"] = parts[0]
    adata.obs["sample_id"] = parts[1] if len(parts) > 1 else np.nan
    adata.obs["read_orientation"] = read_orientation if read_orientation is not None else np.nan
    adata.obs["sample"] = fname

    return adata


adatas = []
for mtx in Path(snakemake.input.liver_raw).glob("*matrix.mtx.gz"):
    if any(pattern in str(mtx) for pattern in snakemake.params.liver_exclude_patterns):
        continue
    adatas.append(read_geo_10x(mtx))


raw_snscadata = ad.concat(
    adatas,
    join="outer",
    label="sample",
    index_unique="-"
)


normedsc = sc.read(snakemake.input.liver_normalized)


# make some edits to make the nomenclature fit
raw_snscadata.obs_names = [
    f"{barcode.split('-')[0]}_{read_orientation}_{sample}"
    for barcode, read_orientation, sample in zip(
        raw_snscadata.obs_names,
        raw_snscadata.obs["read_orientation"],
        raw_snscadata.obs["sample_id"]
    )
]


def extract_orientation(barcode: str) -> str:
    match = re.search(r"(?i)([35]pr(?:v2)?)", str(barcode))
    if match:
        return match.group(1).lower().replace("v2", "")
    return ""


normedsc.obs_names = [
    f"{barcode.split('_')[-1].split('-')[0]}_{extract_orientation(barcode)}_{sample}"
    for barcode, sample in zip(
        normedsc.obs_names,
        normedsc.obs["donor_id"],
    )
]

normedsc.var_names = normedsc.var['feature_name']


# Keep only cells and genes shared between both objects
common_obs = normedsc.obs_names.intersection(raw_snscadata.obs_names)
common_var = normedsc.var_names.intersection(raw_snscadata.var_names)


# Subset both to the same order
adata_sub = normedsc[common_obs, common_var].copy()
raw_sub = raw_snscadata[common_obs, common_var].copy()


assert adata_sub.n_obs == raw_sub.n_obs
assert adata_sub.n_vars == raw_sub.n_vars
assert (adata_sub.obs_names == raw_sub.obs_names).all()
assert (adata_sub.var_names == raw_sub.var_names).all()


adata_sub.layers['counts'] = raw_sub.X.copy()
adata_sub.var = adata_sub.var.drop('feature_name', axis=1)


# Homogenize cell type labels
adata_sub.obs.loc[
    adata_sub.obs["cell_type"].isin(snakemake.params.liver_mesenchymal_types),
    "annotation_MOFA",
] = "mesenchymal"
adata_sub.obs.loc[
    adata_sub.obs["cell_type"].isin(snakemake.params.liver_myeloid_types), "annotation_MOFA"
] = "myeloid"
adata_sub.obs.loc[
    adata_sub.obs["cell_type"].isin(snakemake.params.liver_lymphoid_types), "annotation_MOFA"
] = "lymphoid"
adata_sub.obs.loc[
    adata_sub.obs["cell_type"].isin(snakemake.params.liver_endothelial_types), "annotation_MOFA"
] = "endothelial"
adata_sub.obs.loc[
    adata_sub.obs["cell_type"].isin(snakemake.params.liver_epithelial_types), "annotation_MOFA"
] = "epithelial"
adata_sub.obs["annotation_MOFA"] = adata_sub.obs["annotation_MOFA"].fillna("other")


adata_sub.obs["cond_test"] = "control"
adata_sub.obs.loc[adata_sub.obs["disease"] != snakemake.params.liver_control_condition, "cond_test"] = "fibrosis"


adata_sub.obs = adata_sub.obs.rename(
    columns={
        "cell_type": "cell_type1",
        "sample_uuid": "sample",
        "donor_age": "age",
    }
)


adata_sub.write(snakemake.output.liver_reference)


del raw_snscadata
del normedsc


# Liver Visium slides


def read_geo_visium(folder):
    """
    Read a GEO Visium dataset into an AnnData object.

    Expected structure:
    folder/
        matrix.mtx.gz
        barcodes.tsv.gz
        features.tsv.gz
        spatial/
            tissue_positions_list.csv(.gz)
            scalefactors_json.json(.gz)
            tissue_hires_image.png(.gz)
            tissue_lowres_image.png(.gz)
    """

    folder = Path(folder)
    spatial = folder / "spatial"

    # Read count matrix
    adata = sc.read_10x_mtx(
        folder,
        var_names="gene_symbols",
        make_unique=True,
    )

    # Read coordinates
    posfile = spatial / "tissue_positions_list.csv"
    if not posfile.exists():
        posfile = spatial / "tissue_positions_list.csv.gz"

    coords = pd.read_csv(
        posfile,
        header=None,
        names=[
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ],
    ).set_index("barcode")

    # reorder to match AnnData
    coords = coords.loc[adata.obs_names]

    adata.obs["in_tissue"] = coords["in_tissue"].astype(bool)
    adata.obs["array_row"] = coords["array_row"]
    adata.obs["array_col"] = coords["array_col"]

    adata.obsm["spatial"] = coords[
        ["pxl_col_in_fullres", "pxl_row_in_fullres"]
    ].to_numpy()

    # Read scalefactors
    sf = spatial / "scalefactors_json.json"
    if sf.exists():
        with open(sf) as f:
            scalefactors = json.load(f)
    else:
        with gzip.open(spatial / "scalefactors_json.json.gz", "rt") as f:
            scalefactors = json.load(f)

    # Read images
    def load_image(fname):
        p = spatial / fname
        if p.exists():
            return imread(p)

        pgz = spatial / (fname + ".gz")
        if pgz.exists():
            with gzip.open(pgz, "rb") as f:
                return imread(f)

        return None

    hires = load_image("tissue_hires_image.png")
    lowres = load_image("tissue_lowres_image.png")

    library_id = folder.name

    adata.uns["spatial"] = {
        library_id: {
            "images": {
                "hires": hires,
                "lowres": lowres,
            },
            "scalefactors": scalefactors,
            "metadata": {},
        }
    }

    return adata


for normed, raw_path, normed_path, output_path in zip(
    snakemake.params.liver_samples,
    snakemake.input.liver_visium_raw,
    snakemake.input.liver_visium_normalized,
    snakemake.output.liver_spatial,
):
    adata_raw = read_geo_visium(raw_path)

    adata = sc.read(normed_path)
    adata.var_names = adata.var['feature_name'].str.rsplit('_', n=1).str[0]
    adata.var.index.names = ['index']
    sc.pp.filter_cells(adata, min_genes=snakemake.params.visium_min_genes)
    sc.pp.filter_genes(adata, min_cells=snakemake.params.visium_min_cells)

    adata.var_names_make_unique()

    # Keep only cells and genes shared between both objects
    common_obs = adata.obs_names.intersection(adata_raw.obs_names)
    common_var = adata.var_names.intersection(adata_raw.var_names)

    # Subset both to the same order
    adata_sub = adata[common_obs, common_var].copy()
    raw_sub = adata_raw[common_obs, common_var].copy()

    assert adata_sub.n_obs == raw_sub.n_obs
    assert adata_sub.n_vars == raw_sub.n_vars
    assert (adata_sub.obs_names == raw_sub.obs_names).all()
    assert (adata_sub.var_names == raw_sub.var_names).all()

    assert adata.obs_names.is_unique
    assert adata.var_names.is_unique
    assert adata_raw.obs_names.is_unique
    assert adata_raw.var_names.is_unique

    # Add raw counts as a layer
    adata_sub.layers["counts"] = raw_sub.X.copy()
    adata_sub.layers["normed"] = adata_sub.X.copy()

    adata_sub.obs['sample'] = list(adata_sub.uns['spatial'].keys())[0]

    adata_sub.write(output_path)


# Lung reference


adata = sc.read(snakemake.input.lung_atlas)


# Subset to the selected diagnoses.
adata = adata[adata.obs['diagnosis'].isin(snakemake.params.lung_conditions)]


adata.obs['annotation_MOFA'] = None
# Homogenize cell type labels
adata.obs.loc[
    adata.obs["cell type"].isin(snakemake.params.lung_mesenchymal_types),
    "annotation_MOFA",
] = "mesenchymal"
adata.obs.loc[
    adata.obs["cell type"].isin(snakemake.params.lung_myeloid_types), "annotation_MOFA"
] = "myeloid"
adata.obs.loc[
    adata.obs["cell type"].isin(snakemake.params.lung_lymphoid_types), "annotation_MOFA"
] = "lymphoid"
adata.obs.loc[
    adata.obs["cell type"].isin(snakemake.params.lung_endothelial_types), "annotation_MOFA"
] = "endothelial"
adata.obs.loc[
    adata.obs["cell type"].isin(snakemake.params.lung_epithelial_types), "annotation_MOFA"
] = "epithelial"
adata.obs["annotation_MOFA"] = adata.obs["annotation_MOFA"].fillna("other")


adata.obs["cond_test"] = "control"
adata.obs.loc[adata.obs["diagnosis"] != snakemake.params.lung_control_condition, "cond_test"] = "fibrosis"


adata.obs = adata.obs.rename(
    columns={
        "diagnosis": "grouping",
        "cell type": "cell_type1",
        "cell subtype": "cell_type2",
        "sampleID": "sample",
        "donor_age": "age",
    }
)


adata.var_names = adata.var['gene_name']


adata.write(snakemake.output.lung_reference)


del adata


# Heart reference


input_file = snakemake.input.heart_atlas


# Read AnnData object
adata = sc.read_h5ad(input_file)

# Initialize and standardize AnnData obs columns
adata.obs['annotation_MOFA'] = None
adata.obs['cond_test'] = 'control'
adata.obs['tech'] = "10x 3' v3"
adata.obs = adata.obs.rename(
    columns={
        'sample_id': 'sample',
        'disease_code': 'grouping',
        'cell_type': 'cell_type1'
    }
)

# Set modality and calculate QC metrics
adata.obs['modality'] = 'sn'
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Rename columns for uniform nomenclature across atlases
adata.obs['nCount_RNA'] = adata.X.sum(axis=1)
adata.obs['nFeature_RNA'] = np.count_nonzero(adata.X.astype('int16').toarray(), axis=1)


# Filter cells and genes for quality
adata_sub = adata[adata.obs['pct_counts_mt'] < snakemake.params.heart_mt_cutoff, :]
sc.pp.filter_cells(adata_sub, min_genes=snakemake.params.reference_min_genes)
sc.pp.filter_cells(adata_sub, min_counts=snakemake.params.reference_min_counts)
sc.pp.filter_cells(adata_sub, max_counts=snakemake.params.reference_max_counts)
sc.pp.filter_genes(adata_sub, min_cells=snakemake.params.reference_min_cells)

# Ensure UMAP coordinates are present
if ('X_umap' not in adata_sub.obsm_keys()) and ('UMAP' in adata_sub.obsm_keys()):
    adata_sub.obsm['X_umap'] = adata_sub.obsm['UMAP'].to_numpy()
elif ('X_umap' not in adata_sub.obsm_keys()) and ('UMAP_HARMONY' in adata_sub.obsm_keys()):
    adata_sub.obsm['X_umap'] = adata_sub.obsm['UMAP_HARMONY']

# Unify annotations to common annotation across tissues

# Remove specific regions for Kuppe_2022
for region in snakemake.params.heart_exclude_regions:
    adata_sub = adata_sub[adata_sub.obs['major_labl'] != region]
adata_sub.obs.loc[
    adata_sub.obs['cell_type1'].isin(snakemake.params.heart_pericyte_types),
    'annotation_MOFA'
] = 'pericytesSMCs'
adata_sub.obs.loc[
    adata_sub.obs['cell_type1'].str.contains(snakemake.params.heart_fibroblast_pattern, case=False),
    'annotation_MOFA'
] = 'fibroblast'
adata_sub.obs.loc[
    adata_sub.obs['cell_type1'].isin(snakemake.params.heart_myeloid_types),
    'annotation_MOFA'
] = 'myeloid'
adata_sub.obs.loc[
    adata_sub.obs['cell_type1'].str.contains(snakemake.params.heart_lymphoid_pattern, case=False),
    'annotation_MOFA'
] = 'lymphoid'
adata_sub.obs.loc[
    adata_sub.obs['cell_type1'].str.contains(snakemake.params.heart_endothelial_pattern, case=False),
    'annotation_MOFA'
] = 'endothelial'
adata_sub.obs['cond_test'] = 'control'
adata_sub.obs.loc[
    adata_sub.obs['major_labl'] != snakemake.params.heart_control_region,
    'cond_test'
] = 'fibrosis'
adata_sub.obs['grouping'] = adata_sub.obs['major_labl'].replace(
    {'CTRL': 'control', 'FZ': 'MI-FZ'}
)
adata_sub.layers['counts'] = adata_sub.X.copy()
adata_sub.obs['study'] = 'Kuppe_2022'
adata_sub.obs = adata_sub.obs.rename(columns={'patient_region_id': 'sample'})
del adata_sub.obsm['PCA']
del adata_sub.obsm['HARMONY']

# Set ctype and mesenchymal annotation
adata_sub.obs['ctype'] = adata_sub.obs['annotation_MOFA']
adata_sub.obs.loc[
    adata_sub.obs['annotation_MOFA'].isin(snakemake.params.heart_mesenchymal_types),
    'annotation_MOFA'
] = 'mesenchymal'

adata_sub.obs['annotation_MOFA'] = adata_sub.obs['annotation_MOFA'].fillna('other')

adata_sub.write(snakemake.output.heart_reference)


# Kidney reference


def qc(adata_sub, studyname, mt_cutoff):
    """
    Calculate QC values from an AnnData object and perform basic quality filtering.
    Returns filtered AnnData object and QC values in a separate DataFrame.
    """
    # MT genes were already excluded by authors, but pct is written in obs
    if studyname in (["Lake_2023_sn", "Wilson_2022"]):
        adata_sub.obs = adata_sub.obs.rename(columns={"percent.mt": "pct_counts_mt"})
    else:
        # Calculate QC metrics
        adata_sub.var["mt"] = adata_sub.var_names.str.startswith(
            "MT-"
        )  # annotate mitochondrial genes
        sc.pp.calculate_qc_metrics(
            adata_sub, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )
    # Rename columns for uniform nomenclature across atlases
    adata_sub.obs["nCount_RNA"] = adata_sub.X.sum(axis=1)
    if str(type(adata_sub.X)) in [
        "<class 'scipy.sparse._csr.csr_matrix'>",
        "<class 'scipy.sparse._csc.csc_matrix'>",
    ]:
        adata_sub.obs["nFeature_RNA"] = np.count_nonzero(
            adata_sub.X.astype("int16").toarray(), axis=1
        )
    else:
        adata_sub.obs["nFeature_RNA"] = np.count_nonzero(
            adata_sub.X.astype("int16"), axis=1
        )
    qc_df = adata_sub.obs[["pct_counts_mt", "nFeature_RNA", "nCount_RNA"]]
    qc_df["study"] = studyname
    # Filter cells and genes
    adata_sub = adata_sub[adata_sub.obs["pct_counts_mt"] < mt_cutoff, :]
    sc.pp.filter_cells(adata_sub, min_genes=snakemake.params.reference_min_genes)
    sc.pp.filter_cells(adata_sub, min_counts=snakemake.params.reference_min_counts)
    sc.pp.filter_cells(adata_sub, max_counts=snakemake.params.reference_max_counts)
    sc.pp.filter_genes(adata_sub, min_cells=snakemake.params.reference_min_cells)
    return adata_sub, qc_df


study = "Lake_2023_sc"
adata = sc.read_h5ad(snakemake.input.kidney_atlas)
adata.obs["annotation_MOFA"] = None
adata.obs["cond_test"] = "control"


# subset to diseases of interest
adata_sub = adata[adata.obs["disease"].isin(snakemake.params.kidney_conditions)]
del adata
adata_sub, qc_df = qc(adata_sub, study, snakemake.params.kidney_mt_cutoff)
modality = "sc"

# identify cell types
adata_sub.obs.loc[
    adata_sub.obs["subclass.l2"].isin(snakemake.params.kidney_fibroblast_types),
    "annotation_MOFA",
] = "fibroblast"
adata_sub.obs.loc[adata_sub.obs["subclass.l1"] == snakemake.params.kidney_endothelial_type, "annotation_MOFA"] = (
    "endothelial"
)
adata_sub.obs.loc[
    adata_sub.obs["subclass.l2"].isin(snakemake.params.kidney_pericyte_types),
    "annotation_MOFA",
] = "pericytesSMCs"
adata_sub.obs.loc[adata_sub.obs["subclass.l1"] == snakemake.params.kidney_myeloid_type, "annotation_MOFA"] = (
    "myeloid"
)
adata_sub.obs.loc[
    adata_sub.obs["subclass.l1"].isin(snakemake.params.kidney_lymphoid_types), "annotation_MOFA"
] = "lymphoid"
adata_sub.obs.loc[
    adata_sub.obs["subclass.l1"].isin(
        snakemake.params.kidney_epithelial_types
    ),
    "annotation_MOFA",
] = "epithelial"
adata_sub.obs["sample"] = adata_sub.obs["SampleID"].astype("str")
adata_sub.obs["study"] = study
adata_sub.obs["modality"] = modality
adata_sub.obs["cond_test"] = "control"
adata_sub.obs.loc[adata_sub.obs["disease"] != snakemake.params.kidney_control_condition, "cond_test"] = "fibrosis"
adata_sub.obs["grouping"] = adata_sub.obs["EnrollementCategory"]
adata_sub.obs["grouping"] = adata_sub.obs["grouping"].replace(
    {"LD": "control", "HRT": "control"}
)

adata_sub.obs["tech"] = "10x 3' v3"
adata_sub.obs["ctype"] = adata_sub.obs["annotation_MOFA"]
adata_sub.obs.loc[
    adata_sub.obs["annotation_MOFA"].isin(snakemake.params.kidney_mesenchymal_types),
    "annotation_MOFA",
] = "mesenchymal"
adata_sub.layers["counts"] = adata_sub.X.copy()
adata_sub.obs['annotation_MOFA'] = adata_sub.obs['annotation_MOFA'].fillna('other')
adata_sub.write(snakemake.output.kidney_reference)
