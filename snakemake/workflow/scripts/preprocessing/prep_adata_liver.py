# Preprocessing script for liver single-cell/nucleus data from multiple studies
# Author: Leonie Küchenhoff


import os
import re
from glob import glob

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages
from scipy.sparse import csr_matrix

# Snakemake input/output/params
input_file = snakemake.input[0]
cmap_cells = snakemake.params[0]
cmap_conditions = snakemake.params[1]
study = snakemake.wildcards["study"]
output = snakemake.output.adata
output_pdf = snakemake.output.umap
qc_outputs = snakemake.output.qc_df
ctype_count = snakemake.output.ctype_count
ctype_count_patient = snakemake.output.ctype_count_patient
metadata_file = snakemake.params[2]
anno_file = snakemake.params[3]
meta_columns = snakemake.params[4]


# Helper functions
def convert_age_range_to_mean(value):
    """
    Converts an age range represented as a string in the format 'X-Y Years' or 'X-Y' 
    to the mean of X and Y.
    If the value is NaN or does not match the expected format, 
    it returns the value unchanged.
    """
    if pd.isna(value):
        return value
    if "Years" in value:
        match = re.match(r"(\d+)-(\d+) Years", value)
    else:
        match = re.match(r"(\d+)-(\d+)", value)
    if match:
        lower, upper = int(match.group(1)), int(match.group(2))
        return (lower + upper) / 2
    return value


def read_df(h5, key, skip_check=False):
    """
    Reads a DataFrame from an HDF5 group, handling categorical columns if present.
    """
    h5df = h5[key]
    df = pd.DataFrame(
        index=h5df["index"][()].astype(str), columns=h5df["colnames"][()].astype(str)
    )
    for key in h5df.keys():
        print(key)
        # Check if key is a factor and if so convert
        if key + "_levels" in h5df.keys():
            idx = h5df[key][()].astype("int16")
            df[key] = h5df[key + "_levels"][()][idx].astype(str)
        elif "_levels" not in key and "colnames" not in key and "index" not in key:
            df[key] = h5df[key][()]
    if not skip_check:
        # Check that not all rows are NA
        assert (df.isna().sum(axis=0) != df.shape[0]).any(), (
            "At least one column has all NA values"
        )
    return df


def readh5_to_anndata(file, main_assay="RNA"):
    """
    Reads from h5 file and returns an AnnData object.
    Currently only supports a single type of data (i.e. data OR counts OR scaled data from Seurat).
    Args:
        file (str): Path to h5 file
        main_assay (str): Name of assay in h5 file to use as X
    Returns:
        AnnData: AnnData with main_assay as X and all other assays as layers.
    """
    with h5py.File(file, "r") as f:
        # Check that obs is in f
        if "obs" not in f.keys():
            raise ValueError("Obs not in h5 file")
        # Create AnnData object
        print("Loading main assay %s" % main_assay)
        X = f["assays"][main_assay]
        if X.attrs["datatype"][0] == "SparseMatrix":
            mat = csr_matrix(
                (X["values"], X["indices"], X["indptr"]),
                shape=(X["dims"][0], X["dims"][1]),
            )
            adata = sc.AnnData(X=mat, dtype=np.float64)
            print("INFO: first 20 cells and genes in the X matrix")
            print(adata.X[0:20, 0:20])
            adata.obs_names = X["obs_names"][()].astype(str)
            adata.var_names = X["var_names"][()].astype(str)
        elif X.attrs["datatype"][0] == "Matrix":
            mat = X["data"][()]
            adata = sc.AnnData(
                X=mat,
                dtype=np.float64,
                obs=pd.DataFrame(index=X["obs_names"][()].astype(str)),
                var=pd.DataFrame(index=X["var_names"][()].astype(str)),
            )
        # Add obs
        print("Loading obs")
        obs = read_df(f, "obs")
        adata.obs = obs.loc[adata.obs.index]
        print("INFO: head of obs")
        print(adata.obs.head())
        # Load reductions
        if "reductions" in f.keys():
            print("Loading reductions")
            reductions = f["reductions"]
            for red in reductions.keys():
                print("\tLoading reduction: %s" % red)
                embedding = reductions[red]
                adata.obsm["X_" + red] = np.array(
                    [embedding[key][()] for key in embedding.keys()], dtype=np.float64
                ).T
    return adata


def qc(adata_sub, studyname, mt_cutoff=25):
    """
    Calculate QC values from an AnnData object and perform basic quality filtering.
    Returns filtered AnnData object and QC values in a separate DataFrame.
    """
    if studyname == "Gribben_2024":
        # MT genes were already excluded by authors, but pct is written in obs
        adata_sub.obs = adata_sub.obs.rename(
            columns={"percent.mt.RNA": "pct_counts_mt"}
        )
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
    sc.pp.filter_cells(adata_sub, min_genes=200)
    sc.pp.filter_cells(adata_sub, min_counts=500)
    sc.pp.filter_cells(adata_sub, max_counts=40000)
    sc.pp.filter_genes(adata_sub, min_cells=3)
    return adata_sub, qc_df


def ctype_counter(adata, output, output_patient):
    """
    Summarize cell type counts and save to output files.
    """
    counts = {
        "endothelial": [],
        "myeloid": [],
        "lymphoid": [],
        "mesenchymal": [],
        "epithelial": [],
    }
    ctypes = ["endothelial", "myeloid", "lymphoid", "epithelial", "mesenchymal"]
    for i in ctypes:
        count = adata.obs[adata.obs["annotation_MOFA"] == i].shape[0]
        counts[i] = count
    counts["other"] = adata.obs[~adata.obs["annotation_MOFA"].isin(ctypes)].shape[0]
    counts["pericytesSMCs"] = adata.obs[adata.obs["ctype"] == "pericytesSMCs"].shape[0]
    counts["fibroblast"] = adata.obs[adata.obs["ctype"] == "fibroblast"].shape[0]
    # Save cell type counts
    counts_df = pd.DataFrame({study: counts})
    counts_df.to_csv(output)
    # Also get counts per patient
    if pd.api.types.is_categorical_dtype(adata.obs["annotation_MOFA"]):
        adata.obs["annotation_MOFA"] = adata.obs["annotation_MOFA"].cat.add_categories(
            ["other"]
        )
    obs_data = adata.obs.copy()
    obs_data["annotation_MOFA"].fillna("other", inplace=True)
    counts = obs_data.groupby(["sample", "annotation_MOFA"]).size().reset_index()
    pivoted_df = counts.pivot(index="sample", columns="annotation_MOFA", values=0)
    counts_withcond = pivoted_df.merge(
        obs_data[["sample", "grouping", "cond_test"]], on="sample", how="left"
    ).drop_duplicates()
    counts_withcond.to_csv(output_patient)


def read_txt(path):
    """
    Reads a tab-separated text file, skipping the first row and using the 
    first column as index.
    """
    file = pd.read_csv(path, sep="\t", index_col=0)
    file = file.iloc[1:, :]
    return file


# Read in data

# Ramachandran
if study == "Ramachandran_2019":
    adata = sc.read_10x_mtx(input_file)
    anno = read_txt(anno_file)
    meta = read_txt(metadata_file)

    assert adata.obs.index.equals(meta.index)
    assert meta.index.equals(anno.index)

    adata.obsm["umap"] = np.array(anno[["X", "Y"]]).astype(np.float16)
    adata.obs = meta

    adata_sub = adata[adata.obs["paper"] == "ramachandran"]

    adata_sub = adata_sub[adata_sub.obs["indication"] != "Tumor"]

    adata_sub.obs["annotation_MOFA"] = None

    adata_sub, qc_df = qc(adata_sub, study)
    qc_df.to_csv(qc_outputs)

    # assign cell types
    adata_sub.obs.loc[
        adata_sub.obs["cell_type__ontology_label"] == "stromal cell", "annotation_MOFA"
    ] = "fibroblast"

    adata_sub.obs.loc[
        adata_sub.obs["cell_type__ontology_label"].isin(
            ["hepatocyte", "epithelial cell"]
        ),
        "annotation_MOFA",
    ] = "epithelial"

    adata_sub.obs.loc[
        adata_sub.obs["cell_type__ontology_label"].isin(["endothelial cell"]),
        "annotation_MOFA",
    ] = "endothelial"
    adata_sub.obs.loc[
        adata_sub.obs["cell_type__ontology_label"].isin(
            ["myeloid cell", "mast cell", "erythrocyte", "plasmacytoid dendritic cell"]
        ),
        "annotation_MOFA",
    ] = "myeloid"

    adata_sub.obs.loc[
        adata_sub.obs["cell_type__ontology_label"].isin(["lymphocyte", "B cell"]),
        "annotation_MOFA",
    ] = "lymphoid"

    adata_sub.obs["study"] = study
    adata_sub.obs["cond_test"] = "control"
    adata_sub.obs.loc[adata_sub.obs["indication"] != "healthy", "cond_test"] = (
        "fibrosis"
    )
    adata_sub.obs["grouping"] = adata_sub.obs["indication"]

    adata_sub.obs["grouping"] = adata_sub.obs["grouping"].replace(
        {
            "healthy": "control",
            "Cirrhotic": "cirrhosis",
            "alcohol": "ALD",
            "mild_steatosis": "mild-steatosis",
            "low_steatosis": "low-steatosis",
            "NAFLD": "MASLD",
        }
    )
    adata_sub.obs = adata_sub.obs.rename(
        columns={
            "cell_type__ontology_label": "cell_type1",
            "library_preparation_protocol__ontology_label": "tech",
        }
    )

    adata_sub.obs["ctype"] = adata_sub.obs["annotation_MOFA"]
    adata_sub.obs.loc[
        adata_sub.obs["annotation_MOFA"].isin(["pericytesSMCs", "fibroblast"]),
        "annotation_MOFA",
    ] = "mesenchymal"

    ctype_counter(adata_sub, ctype_count, ctype_count_patient)
    adata_sub.layers["counts"] = adata_sub.X.copy()
    adata_sub.obs["modality"] = "sc"

    adata_sub.obs = adata_sub.obs.drop(["organ"], axis=1)

    # Keep only relevant metadata columns
    cols_to_keep = list(set(adata_sub.obs.columns) & set(meta_columns))
    adata_sub.obs = adata_sub.obs[cols_to_keep]

    adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()

    adata_sub.write(output)


# Gribben
if study == "Gribben_2024":
    adata_sub = readh5_to_anndata(input_file)
    adata_sub.obs["annotation_MOFA"] = None

    adata_sub, qc_df = qc(adata_sub, study)
    qc_df.to_csv(qc_outputs)

    # assign cell types
    adata_sub.obs.loc[
        adata_sub.obs["cell.annotation"] == "Stellate", "annotation_MOFA"
    ] = "fibroblast"
    adata_sub.obs.loc[
        adata_sub.obs["cell.annotation"] == "Endothelial", "annotation_MOFA"
    ] = "endothelial"
    adata_sub.obs.loc[
        adata_sub.obs["cell.annotation"].isin(["Macrophages", "Neutrophils"]),
        "annotation_MOFA",
    ] = "myeloid"
    adata_sub.obs.loc[
        adata_sub.obs["cell.annotation"].isin(["B-cell 1", "B-cell 2", "Lymphocytes"]),
        "annotation_MOFA",
    ] = "lymphoid"
    adata_sub.obs.loc[
        adata_sub.obs["cell.annotation"].isin(["Cholangiocytes", "Hepatocytes"]),
        "annotation_MOFA",
    ] = "epithelial"

    adata_sub.obs["study"] = study
    adata_sub.obs["cond_test"] = "control"
    adata_sub.obs.loc[
        adata_sub.obs["Disease.status"] != "Healthy control", "cond_test"
    ] = "fibrosis"
    adata_sub.obs["grouping"] = adata_sub.obs["Disease.status"]

    adata_sub.obs["grouping"] = adata_sub.obs["grouping"].replace(
        {
            "Healthy control": "control",
            "NAFLD": "MASLD",
            "NASH w/o cirrhosis": "MASH w/o cirrhosis",
            "NASH with cirrhosis": "MASH",
            "end stage": "cirrhosis",
        }
    )
    adata_sub.obs = adata_sub.obs.rename(
        columns={
            "cell.annotation": "cell_type1",
            "Gender": "sex",
            "Lobe": "region",
            "orig.ident": "sample",
            "Age": "age",
            "Fibrosis.score..F0.4.": "Fibrosis",
        }
    )
    adata_sub.obs["sex"] = adata_sub.obs["sex"].replace({"M": "Male", "F": "Female"})

    adata_sub.obs["ctype"] = adata_sub.obs["annotation_MOFA"]
    adata_sub.obs.loc[
        adata_sub.obs["annotation_MOFA"].isin(["pericytesSMCs", "fibroblast"]),
        "annotation_MOFA",
    ] = "mesenchymal"

    ctype_counter(adata_sub, ctype_count, ctype_count_patient)
    adata_sub.layers["counts"] = adata_sub.X.copy()
    adata_sub.obs["modality"] = "sn"
    adata_sub.obs["tech"] = "10x 3' v3.1"

    # Keep only relevant metadata columns
    cols_to_keep = list(set(adata_sub.obs.columns) & set(meta_columns))
    adata_sub.obs = adata_sub.obs[cols_to_keep]

    adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()

    adata_sub.write(output)


# Watson
if study == "Watson_2025":
    # Only contains the normed data, but also the metadata
    # Therefore, the raw counts also need to be read
    watson = sc.read(metadata_file)
    samples = pd.read_csv(anno_file)

    # Get all .csv.gz files in the folder
    full_files = glob(os.path.join(input_file, "*.mtx.gz"))
    files = [os.path.basename(i) for i in full_files]
    rel_files = [
        file for file in files if file.split("_")[0] in samples["accession"].unique()
    ]

    adatas = []
    for file in rel_files:
        prefix = re.match(r"^(GSM\d+_\d+-\d+)_", file).group(1)
        sample_name = samples.loc[
            samples["accession"] == prefix.split("_")[0], "name"
        ].item()

        adata = sc.read_10x_mtx(input_file, prefix=f"{prefix}_")
        adata.obs["file"] = prefix
        adata.obs["sample"] = sample_name

        adata.obs_names = sample_name + "_" + adata.obs_names
        adatas.append(adata)

    adata_combined = ad.concat(adatas, join="outer")
    # Unique IDs are annotated differently
    watson.obs_names = [i[:-2].replace(".", "-") for i in watson.obs_names]
    new_obs = adata_combined.obs.merge(
        watson.obs, left_index=True, right_index=True, how="left"
    )
    assert np.all(new_obs.index == adata_combined.obs_names)
    adata_combined.obs = new_obs
    # Remove cells without metadata information
    adata_sub = adata_combined[watson.obs_names, :]
    assert np.all(adata_sub.obs_names == watson.obs_names)
    adata_sub.obsm = watson.obsm

    adata_sub.obs["cell_type1"] = (
        adata_sub.obs["cell_type_final_injured"]
        .astype("object")
        .fillna(adata_sub.obs["cell_type_final_healthy"].astype("object"))
    )
    adata_sub.obs["annotation_MOFA"] = None

    adata_sub, qc_df = qc(adata_sub, study)
    qc_df.to_csv(qc_outputs)

    # assign cell types
    adata_sub.obs.loc[
        adata_sub.obs["cell_type1"].isin(["HSC", "HSC_1", "HSC_2"]), "annotation_MOFA"
    ] = "fibroblast"

    adata_sub.obs.loc[adata_sub.obs["cell_type1"].isin(["LSEC"]), "annotation_MOFA"] = (
        "endothelial"
    )

    adata_sub.obs.loc[adata_sub.obs["cell_type1"].isin(["VSMC"]), "annotation_MOFA"] = (
        "pericytesSMCs"
    )

    adata_sub.obs.loc[
        adata_sub.obs["cell_type1"].isin(["Mac_2", "Mac_1", "Mac"]), "annotation_MOFA"
    ] = "myeloid"

    adata_sub.obs.loc[
        adata_sub.obs["cell_type1"].isin(["Immune Cells"]), "annotation_MOFA"
    ] = "lymphoid"

    adata_sub.obs.loc[
        adata_sub.obs["cell_type1"].isin(
            [
                "Central_Hep",
                "Portal_Hep",
                "IJ_1_Hep",
                "IJ_2_Hep",
                "Cholangiocyte",
                "Hep_1",
                "Hep_2",
                "Hep_3",
            ]
        ),
        "annotation_MOFA",
    ] = "epithelial"

    adata_sub.obs["study"] = study
    adata_sub.obs["cond_test"] = "control"
    adata_sub.obs.loc[adata_sub.obs["Condition"] != "Normal", "cond_test"] = "fibrosis"
    adata_sub.obs["grouping"] = adata_sub.obs["Condition"]

    adata_sub.obs["grouping"] = adata_sub.obs["grouping"].replace(
        {"Normal": "control", "Disease": "fibrosis"}
    )
    adata_sub.obs["batch"] = ["batch_" + str(i) for i in adata_sub.obs["batch"]]
    adata_sub.obs = adata_sub.obs.rename(columns={"Sex": "sex"})
    adata_sub.obs["sex"] = adata_sub.obs["sex"].replace({"F": "Female", "M": "Male"})
    adata_sub.obs["age"] = adata_sub.obs["Age"].apply(convert_age_range_to_mean)
    adata_sub.obs["ctype"] = adata_sub.obs["annotation_MOFA"]
    adata_sub.obs.loc[
        adata_sub.obs["annotation_MOFA"].isin(["pericytesSMCs", "fibroblast"]),
        "annotation_MOFA",
    ] = "mesenchymal"

    ctype_counter(adata_sub, ctype_count, ctype_count_patient)
    adata_sub.layers["counts"] = adata_sub.X.copy()
    adata_sub.obs["modality"] = "sn"
    adata_sub.obs["tech"] = "10x 3' v3.1"

    # Keep only relevant metadata columns
    cols_to_keep = list(set(adata_sub.obs.columns) & set(meta_columns))
    adata_sub.obs = adata_sub.obs[cols_to_keep]

    adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()

    adata_sub.write(output)


# Buonomo
if study == "Buonomo_2022":
    files = glob(os.path.join(input_file, "*.csv.gz"))
    # List to store individual AnnData objects
    adatas = []

    for file in files:
        # Read the count matrix
        df = pd.read_csv(file, index_col=0).astype("float32")

        # Convert to AnnData object
        adata = ad.AnnData(df.T)

        # Optionally add a batch column based on the filename
        adata.obs["batch"] = os.path.basename(file)

        adatas.append(adata)

    # Concatenate all AnnData objects
    adata_combined = ad.concat(adatas, join="outer", merge="first")

    mdata = pd.read_csv(metadata_file, index_col=0)
    new_obs = adata_combined.obs.merge(
        mdata, left_index=True, right_index=True, how="left"
    )
    assert np.all(new_obs.index == adata_combined.obs.index)
    adata_combined.obs = new_obs
    # Remove cells without metadata information
    adata_sub = adata_combined[~adata_combined.obs["UMAP_1"].isna()]
    adata_sub.obs["annotation_MOFA"] = None

    adata_sub, qc_df = qc(adata_sub, study)
    qc_df.to_csv(qc_outputs)

    # assign cell types
    adata_sub.obs.loc[
        (adata_sub.obs["ano_l1"] == "Mesenchyme")
        & (adata_sub.obs["ano_l2"].isin(["Mes-1", "Mes-2", "Mesothelia", "HSC"])),
        "annotation_MOFA",
    ] = "fibroblast"
    adata_sub.obs.loc[
        adata_sub.obs["ano_l1"].isin(["Endothelial"]), "annotation_MOFA"
    ] = "endothelial"
    adata_sub.obs.loc[
        adata_sub.obs["ano_l2"].isin(["Arterial EC", "Vascular EC-1", "Vascular EC-2"]),
        "annotation_MOFA",
    ] = "endothelial"
    adata_sub.obs.loc[adata_sub.obs["ano_l2"].isin(["VSMC"]), "annotation_MOFA"] = (
        "pericytesSMCs"
    )
    adata_sub.obs.loc[
        adata_sub.obs["ano_l1"].isin(
            ["Macrophage", "Monocytes", "mDC", "Erythroid", "pDCs"]
        ),
        "annotation_MOFA",
    ] = "myeloid"
    adata_sub.obs.loc[
        adata_sub.obs["ano_l1"].isin(["NK", "Plasma Bcells", "Bcells", "CD4", "CD8"]),
        "annotation_MOFA",
    ] = "lymphoid"
    adata_sub.obs.loc[
        adata_sub.obs["ano_l1"].isin(["Epithelial", "Hepatocyte"]), "annotation_MOFA"
    ] = "epithelial"
    adata_sub.obs.loc[
        adata_sub.obs["ano_l2"].isin(["Cholangiocytes"]), "annotation_MOFA"
    ] = "epithelial"

    adata_sub.obs["study"] = study
    adata_sub.obs["cond_test"] = "control"
    adata_sub.obs.loc[adata_sub.obs["fraction"] != "NonCirrhotic", "cond_test"] = (
        "fibrosis"
    )
    adata_sub.obs["grouping"] = adata_sub.obs["fraction"]

    adata_sub.obs["grouping"] = adata_sub.obs["grouping"].replace(
        {"NonCirrhotic": "control", "Cirrhotic": "cirrhosis"}
    )
    adata_sub.obs = adata_sub.obs.rename(
        columns={"ano_l1": "cell_type1", "ano_l2": "cell_type2"}
    )
    adata_sub.obs["ctype"] = adata_sub.obs["annotation_MOFA"]
    adata_sub.obs.loc[
        adata_sub.obs["annotation_MOFA"].isin(["pericytesSMCs", "fibroblast"]),
        "annotation_MOFA",
    ] = "mesenchymal"

    ctype_counter(adata_sub, ctype_count, ctype_count_patient)
    adata_sub.layers["counts"] = adata_sub.X.copy()
    adata_sub.obs["modality"] = "sc"
    adata_sub.obs["tech"] = "10x 3' v2"

    # Keep only relevant metadata columns
    cols_to_keep = list(set(adata_sub.obs.columns) & set(meta_columns))
    adata_sub.obs = adata_sub.obs[cols_to_keep]

    adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()

    adata_sub.write(output)


if "X_umap" in adata_sub.obsm:
    print(f"embedding used from authors for {study}")
else:
    # UMAP embedding still needs to be generated
    # Normalize the data
    sc.pp.normalize_total(adata_sub, target_sum=1e4)
    sc.pp.log1p(adata_sub)
    # Identify the 2000 most highly variable genes
    sc.pp.highly_variable_genes(adata_sub, min_mean=0.0125, max_mean=3, min_disp=0.5)
    # Filter higly variable genes
    adata_sub.raw = adata_sub
    adata_sub = adata_sub[:, adata_sub.var.highly_variable]

    # Regress and scale the data
    sc.pp.regress_out(adata_sub, ["nCount_RNA", "pct_counts_mt"])
    sc.pp.scale(adata_sub, max_value=10)

    # Generate PCA features
    sc.tl.pca(adata_sub, svd_solver="arpack")

    # Find cell neighbors
    sc.pp.neighbors(adata_sub, n_neighbors=25, n_pcs=30)

    # Generate UMAP features
    sc.tl.umap(adata_sub)

    adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()


with PdfPages(output_pdf) as output_pdf:
    for i in [adata_sub, adata_new]:
        fig, axs = plt.subplots(1, 3, tight_layout=True, figsize=(15, 4))
        sc.pl.umap(
            i, color=["annotation_MOFA"], ax=axs[0], show=False, palette=cmap_cells
        )
        sc.pl.umap(
            i, color=["cond_test"], ax=axs[1], show=False, palette=cmap_conditions
        )
        sc.pl.umap(i, color=["sample"], ax=axs[2], show=False, legend_loc=None)
        output_pdf.savefig(fig)
