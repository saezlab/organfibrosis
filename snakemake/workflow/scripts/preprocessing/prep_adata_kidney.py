# Preprocessing script for kidney single-cell/nucleus data from multiple studies
# Author: Leonie Küchenhoff

import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages

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
meta_columns = snakemake.params[3]

# Load data
adata = sc.read_h5ad(input_file)
adata.obs["annotation_MOFA"] = None
adata.obs["cond_test"] = "control"


# Utility functions
def convert_age_range_to_mean(value):
    """
    Converts an age range represented as a string in the format 'X-Y Years' or 'X-Y' to the mean of X and Y.
    If the value is NaN or does not match the expected format, it returns the value unchanged.
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


def convert_egfr_range_to_mean(gfr):
    """
    Converts an eGFR range represented as a string in the format 'XX-YY ml/min/1.73m2', to the mean of X and Y.
    If the value is NaN or does not match the expected format, it returns np.nan.
    """
    if pd.isna(gfr):
        return np.nan
    elif gfr == ">60":
        return 80
    try:
        range_values = gfr.split(" ")[0]  # Get the part before the first space
        low, high = map(float, range_values.split("-"))
        return (low + high) / 2
    except:
        return np.nan


def qc(adata_sub, studyname, mt_cutoff=25):
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


# Preprocessing pipeline
if study == "Lake_2023_sc":
    # subset to diseases of interest
    adata_sub = adata[adata.obs["disease"].isin(["chronic kidney disease", "normal"])]
    del adata
    adata_sub, qc_df = qc(adata_sub, study, 35)
    modality = "sc"
    qc_df.to_csv(qc_outputs)
    # identify cell types
    adata_sub.obs.loc[
        adata_sub.obs["subclass.l2"].isin(["MYOF", "aFIB", "FIB", "dFIB"]),
        "annotation_MOFA",
    ] = "fibroblast"
    adata_sub.obs.loc[adata_sub.obs["subclass.l1"] == "EC", "annotation_MOFA"] = (
        "endothelial"
    )
    adata_sub.obs.loc[
        adata_sub.obs["subclass.l2"].isin(["VSMC", "MC", "REN", "VSMC/P", "dVSMC"]),
        "annotation_MOFA",
    ] = "pericytesSMCs"
    adata_sub.obs.loc[adata_sub.obs["subclass.l1"] == "MYL", "annotation_MOFA"] = (
        "myeloid"
    )
    adata_sub.obs.loc[
        adata_sub.obs["subclass.l1"].isin(["T", "B"]), "annotation_MOFA"
    ] = "lymphoid"
    adata_sub.obs.loc[
        adata_sub.obs["subclass.l1"].isin(
            ["PC", "DTL/ATL", "TAL", "IC", "PT", "DCT/CNT", "POD/PEC"]
        ),
        "annotation_MOFA",
    ] = "epithelial"
    adata_sub.obs["sample"] = adata_sub.obs["SampleID"].astype("str")
    adata_sub.obs["study"] = study
    adata_sub.obs["modality"] = modality
    adata_sub.obs["cond_test"] = "control"
    adata_sub.obs.loc[adata_sub.obs["disease"] != "normal", "cond_test"] = "fibrosis"
    adata_sub.obs["grouping"] = adata_sub.obs["EnrollementCategory"]
    adata_sub.obs["grouping"] = adata_sub.obs["grouping"].replace(
        {"LD": "control", "HRT": "control"}
    )
    # convert binned age into mean age of bin
    adata_sub.obs["age"] = adata_sub.obs["Age"].apply(convert_age_range_to_mean)
    # convert binned eGFR into mean eGFR of bin
    adata_sub.obs["eGFR"] = adata_sub.obs["eGFR"].apply(convert_egfr_range_to_mean)
    adata_sub.obs = adata_sub.obs.rename(
        columns={
            "subclass.l1": "cell_type1",
            "subclass.l2": "cell_type2",
            "Gender": "sex",
        }
    )
    adata_sub.obs["tech"] = "10x 3' v3"
    adata_sub.obs["ctype"] = adata_sub.obs["annotation_MOFA"]
    adata_sub.obs.loc[
        adata_sub.obs["annotation_MOFA"].isin(["pericytesSMCs", "fibroblast"]),
        "annotation_MOFA",
    ] = "mesenchymal"
    ctype_counter(adata_sub, ctype_count, ctype_count_patient)
    adata_sub.layers["counts"] = adata_sub.X.copy()
    # keep only relevant metadata columns
    cols_to_keep = list(set(adata_sub.obs.columns) & set(meta_columns))
    adata_sub.obs = adata_sub.obs[cols_to_keep]
    adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()
    adata_sub.write(output)

if study == "Lake_2023_sn":
    # subset to diseases of interest
    adata_sub = adata[adata.obs["condition_1"].isin(["HR", "CKD"])]
    del adata
    adata_sub, qc_df = qc(adata_sub, study)
    modality = "sn"
    qc_df.to_csv(qc_outputs)
    # identify cell types
    adata_sub.obs.loc[adata_sub.obs["subclass.l1"] == "FIB", "annotation_MOFA"] = (
        "fibroblast"
    )
    adata_sub.obs.loc[adata_sub.obs["subclass.l1"] == "EC", "annotation_MOFA"] = (
        "endothelial"
    )
    adata_sub.obs.loc[adata_sub.obs["subclass.l1"] == "VSM/P", "annotation_MOFA"] = (
        "pericytesSMCs"
    )
    adata_sub.obs.loc[
        adata_sub.obs["subclass.l2"].isin(
            ["MAC-M2", "N", "MDC", "cycMNP", "ncMON", "cDC", "MAST", "pDC"]
        ),
        "annotation_MOFA",
    ] = "myeloid"
    adata_sub.obs.loc[
        adata_sub.obs["subclass.l2"].isin(["T", "B", "NKT", "PL"]), "annotation_MOFA"
    ] = "lymphoid"
    adata_sub.obs.loc[
        adata_sub.obs["class"] == "epithelial cells", "annotation_MOFA"
    ] = "epithelial"
    adata_sub.obs["sample"] = (
        adata_sub.obs["patient"].astype("str")
        + "_"
        + adata_sub.obs["region"].astype("str")
    )
    adata_sub.obs["study"] = study
    adata_sub.obs["modality"] = modality
    adata_sub.obs["cond_test"] = "control"
    adata_sub.obs.loc[adata_sub.obs["condition_1"] != "HR", "cond_test"] = "fibrosis"
    adata_sub.obs["grouping"] = adata_sub.obs["condition_2"]
    adata_sub.obs["grouping"] = adata_sub.obs["grouping"].replace(
        {
            "HR-B": "control",
            "HR": "control",
            "HR-D": "control",
            "HR-N": "control",
            "HR-T": "control",
        }
    )
    # convert binned age into mean age of bin
    adata_sub.obs["age"] = adata_sub.obs["clin_age_bin"].apply(
        convert_age_range_to_mean
    )
    # convert binned eGFR into mean eGFR of bin
    adata_sub.obs["eGFR"] = adata_sub.obs["clin_base_egfr_bin"].apply(
        convert_egfr_range_to_mean
    )
    adata_sub.obs = adata_sub.obs.rename(
        columns={
            "subclass.l1": "cell_type1",
            "subclass.l2": "cell_type2",
            "clin_sex": "sex",
            "tech_protocol": "tech",
        }
    )
    adata_sub.obs["ctype"] = adata_sub.obs["annotation_MOFA"]
    adata_sub.obs.loc[
        adata_sub.obs["annotation_MOFA"].isin(["pericytesSMCs", "fibroblast"]),
        "annotation_MOFA",
    ] = "mesenchymal"
    ctype_counter(adata_sub, ctype_count, ctype_count_patient)
    adata_sub.layers["counts"] = adata_sub.X.copy()
    # keep only relevant metadata columns
    cols_to_keep = list(set(adata_sub.obs.columns) & set(meta_columns))
    adata_sub.obs = adata_sub.obs[cols_to_keep]
    adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()
    adata_sub.write(output)

if study in ["Abedini_2024_sn", "Abedini_2024_sc"]:
    meta = pd.read_csv(metadata_file, skiprows=1, usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8])
    adata.obs = (
        adata.obs.reset_index()
        .merge(
            meta,
            how="left",
            left_on=["clin_sex", "clin_age", "uid_donor_id"],
            right_on=["Sex", "Age", "ID"],
        )
        .set_index("index")
    )
    adata_sub = adata.copy()
    del adata
    if study == "Abedini_2024_sc":
        adata_sub, qc_df = qc(adata_sub, "Abedini_2024_sc", 35)
        modality = "sc"
    elif study == "Abedini_2024_sn":
        adata_sub, qc_df = qc(adata_sub, "Abedini_2024_sn")
        modality = "sn"
    qc_df.to_csv(qc_outputs)
    # identify cell types
    adata_sub.obs.loc[adata_sub.obs["celltype_1"] == "FIB", "annotation_MOFA"] = (
        "fibroblast"
    )
    adata_sub.obs.loc[adata_sub.obs["celltype_1"] == "EC", "annotation_MOFA"] = (
        "endothelial"
    )
    adata_sub.obs.loc[adata_sub.obs["celltype_1"] == "VSM_P", "annotation_MOFA"] = (
        "pericytesSMCs"
    )
    adata_sub.obs.loc[
        adata_sub.obs["celltype_2"].isin(
            ["MAC-M2", "N", "MDC", "cycMNP", "ncMON", "cDC", "MAST", "pDC"]
        ),
        "annotation_MOFA",
    ] = "myeloid"
    adata_sub.obs.loc[
        adata_sub.obs["celltype_2"].isin(["T", "B", "PL", "NKC_T", "NKT", "cycNKC_T"]),
        "annotation_MOFA",
    ] = "lymphoid"
    adata_sub.obs.loc[
        adata_sub.obs["celltype_1"].isin(
            ["PT", "IC", "TAL", "PC", "PapE", "CNT", "DCT", "POD", "DTL", "ATL", "PEC"]
        ),
        "annotation_MOFA",
    ] = "epithelial"
    adata_sub.obs["study"] = study
    adata_sub.obs["modality"] = modality
    adata_sub.obs["cond_test"] = "control"
    adata_sub.obs.loc[adata_sub.obs["condition_1"] != "HR", "cond_test"] = "fibrosis"
    adata_sub.obs["grouping"] = adata_sub.obs["Group"]
    adata_sub.obs["grouping"] = adata_sub.obs["grouping"].replace(
        {"HKD": "H-CKD", "Control": "control"}
    )
    adata_sub.obs = adata_sub.obs.rename(
        columns={
            "celltype_1": "cell_type1",
            "celltype_2": "cell_type2",
            "clin_sex": "sex",
            "clin_age": "age",
            "uid_donor_id": "sample",
        }
    )
    adata_sub.obs["sample"] = (
        adata_sub.obs["sample"].astype("str")
        + "_"
        + adata_sub.obs["modality"].astype("str")
    )
    adata_sub.obs["tech"] = "10x 3' v3"
    adata_sub.obs["ctype"] = adata_sub.obs["annotation_MOFA"]
    adata_sub.obs.loc[
        adata_sub.obs["annotation_MOFA"].isin(["pericytesSMCs", "fibroblast"]),
        "annotation_MOFA",
    ] = "mesenchymal"
    ctype_counter(adata_sub, ctype_count, ctype_count_patient)
    adata_sub.layers["counts"] = adata_sub.X.copy()
    # keep only relevant metadata columns
    cols_to_keep = list(set(adata_sub.obs.columns) & set(meta_columns))
    adata_sub.obs = adata_sub.obs[cols_to_keep]
    adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()
    adata_sub.write(output)

if study == "McCown_2025_sn":
    # subset to diseases of interest
    adata_sub = adata[adata.obs["Prep"] == "SN"].copy()
    meta = pd.read_csv(metadata_file)
    # Merge with metadata, but only with non-existing columns
    new_columns = list(set(meta.columns) - set(adata_sub.obs.columns))
    adata_sub.obs = (
        adata_sub.obs.reset_index()
        .merge(meta[new_columns + ["EdgarID"]], how="left", on=["EdgarID"])
        .set_index("index")
    )
    adata_sub = adata_sub[adata_sub.obs["Cohort"] != "N/A"]
    adata_sub.var_names = list(adata_sub.var.index)
    del adata
    adata_sub, qc_df = qc(adata_sub, study)
    modality = "sn"
    qc_df.to_csv(qc_outputs)
    # identify cell types
    adata_sub.obs.loc[adata_sub.obs["LR_Cluster"] == "EC", "annotation_MOFA"] = (
        "endothelial"
    )
    adata_sub.obs.loc[
        (adata_sub.obs["LR_Cluster"] == "vSMC/MC")
        | (adata_sub.obs["LR_Cluster"] == "PC"),
        "annotation_MOFA",
    ] = "pericytesSMCs"
    adata_sub.obs.loc[adata_sub.obs["LR_Cluster"] == "FIB", "annotation_MOFA"] = (
        "fibroblast"
    )
    adata_sub.obs.loc[adata_sub.obs["celltype"] == "MAC/MON", "annotation_MOFA"] = (
        "myeloid"
    )
    adata_sub.obs.loc[(adata_sub.obs["celltype"] == "NKC/NKT"), "annotation_MOFA"] = (
        "lymphoid"
    )
    adata_sub.obs.loc[
        adata_sub.obs["LR_Cluster"].isin(
            [
                "PT",
                "PC",
                "DTL/aPT",
                "IC",
                "PEC/POD",
                "DCT",
                "TAL",
                "POD",
                "ATL",
                "PEC",
                "CNT",
                "tPC-IC",
                "cycDTL",
            ]
        ),
        "annotation_MOFA",
    ] = "epithelial"
    adata_sub.obs["sample"] = adata_sub.obs["ID"].astype("str") + adata_sub.obs[
        "Method"
    ].astype("str")
    adata_sub.obs["study"] = study
    adata_sub.obs["modality"] = modality
    adata_sub.obs["grouping"] = list(adata_sub.obs["Cohort"])
    adata_sub.obs.loc[adata_sub.obs["Cohort"].str.contains("IgAN"), "grouping"] = "IgA"
    adata_sub.obs.loc[adata_sub.obs["Cohort"].str.contains("LD"), "grouping"] = (
        "control"
    )
    adata_sub.obs.loc[adata_sub.obs["Cohort"].str.contains("Trans"), "grouping"] = (
        "control"
    )
    adata_sub.obs.loc[adata_sub.obs["Cohort"] == "T2D", "grouping"] = "DKD"
    adata_sub.obs.loc[adata_sub.obs["Cohort"].str.contains("FSGS"), "grouping"] = "FSGS"
    adata_sub.obs.loc[
        adata_sub.obs["Cohort"].isin(
            [
                "AIN/T2D",
                "C1QN",
                "C3_Glomerulopathy",
                "COVAN",
                "Immune_Cmpx_GN",
                "MN_Vax",
                "MPGN",
            ]
        ),
        "grouping",
    ] = "Other"
    adata_sub.obs["cond_test"] = "control"
    adata_sub.obs.loc[adata_sub.obs["grouping"] != "control", "cond_test"] = "fibrosis"
    adata_sub.obs = adata_sub.obs.rename(
        columns={
            "LR_Cluster": "cell_type1",
            "celltype": "cell_type2",
            "Sex": "sex",
            "Age": "age",
            "eGFRatBx_NEPTUNE": "eGFR",
            "Project1": "batch",
            "Method": "tech",
            "InterstitialFibrosis":"Fibrosis"
        }
    )
    # replace string 'N/A' with actual nan values for easier data handeling down the line
    adata_sub.obs["age"] = adata_sub.obs["age"].replace("N/A", np.nan)
    adata_sub.obs["sex"] = adata_sub.obs["sex"].replace("N/A", np.nan)
    adata_sub.obs["eGFR"] = adata_sub.obs["eGFR"].replace("N/A", np.nan)
    adata_sub.obs["ctype"] = adata_sub.obs["annotation_MOFA"]
    adata_sub.obs.loc[
        adata_sub.obs["annotation_MOFA"].isin(["pericytesSMCs", "fibroblast"]),
        "annotation_MOFA",
    ] = "mesenchymal"
    ctype_counter(adata_sub, ctype_count, ctype_count_patient)
    adata_sub.layers["counts"] = adata_sub.X.copy()
    # keep only relevant metadata columns
    cols_to_keep = list(set(adata_sub.obs.columns) & set(meta_columns))
    adata_sub.obs = adata_sub.obs[cols_to_keep]
    adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()
    adata_sub.write(output)

if study == "Muto_2022":
    adata_sub, qc_df = qc(adata, "Muto_2022")
    del adata
    qc_df.to_csv(qc_outputs)
    # identify cell types
    adata_sub.obs.loc[
        adata_sub.obs["cell_type2_pred"].isin(
            ["MYOF", "M-FIB", "aFIB", "dFIB", "FIB", "dM-FIB"]
        ),
        "annotation_MOFA",
    ] = "fibroblast"
    adata_sub.obs.loc[
        adata_sub.obs["cell_type2_pred"].isin(
            ["EC-PTC", "EC-GC", "EC-DVR", "EC-AEA", "EC-LYM", "EC-AVR", "cycEC"]
        ),
        "annotation_MOFA",
    ] = "endothelial"
    adata_sub.obs.loc[
        adata_sub.obs["cell_type2_pred"].isin(["REN", "MC", "VSMC/P", "dVSMC", "VSMC"]),
        "annotation_MOFA",
    ] = "pericytesSMCs"
    adata_sub.obs.loc[
        adata_sub.obs["cell_type2_pred"].isin(
            ["MAC-M2", "N", "MDC", "cycMNP", "ncMON", "cDC", "MAST"]
        ),
        "annotation_MOFA",
    ] = "myeloid"
    adata_sub.obs.loc[
        adata_sub.obs["cell_type2_pred"].isin(["T", "B", "PL", "NKT"]),
        "annotation_MOFA",
    ] = "lymphoid"
    adata_sub.obs.loc[
        adata_sub.obs["cell_type2_pred"].isin(
            [
                "PT-S1",
                "C-TAL",
                "DCT2",
                "IC-B",
                "aPT",
                "CCD-IC-A",
                "CNT-PC",
                "DCT1",
                "PT-S3",
                "aTAL1",
                "dC-IC-A",
                "POD",
                "M-TAL",
                "DTL1",
                "dM-TAL",
                "CNT-IC-A",
                "DTL2",
                "dPOD",
                "dOMCD-PC",
                "dPT/DTL",
                "aTAL2",
                "OMCD-PC",
                "tPC-IC",
                "cycPT",
                "OMCD-IC-A",
                "IMCD",
                "dIMCD",
                "MD",
                "CNT",
                "CCD-PC",
                "PEC",
                "PT-S2",
                "dCNT",
                "dPT",
            ]
        ),
        "annotation_MOFA",
    ] = "epithelial"
    # adapt nomenclature to make unify with other studies
    adata_sub.obs["study"] = "Muto_2022"
    adata_sub.obs["modality"] = "sn"
    adata_sub.obs["cond_test"] = "control"
    adata_sub.obs.loc[adata_sub.obs["disease"] != "control", "cond_test"] = "fibrosis"
    adata_sub.obs["grouping"] = adata_sub.obs["disease"]
    adata_sub.obs = adata_sub.obs.drop(columns=["cell_type2"])
    adata_sub.obs = adata_sub.obs.rename(
        columns={
            "celltype": "cell_type1",
            "cell_type2_pred": "cell_type2",
            "patient": "sample",
            "gender": "sex",
        }
    )
    adata_sub.obs["tech"] = "10x 5' v2"
    adata_sub.obs["ctype"] = adata_sub.obs["annotation_MOFA"]
    adata_sub.obs.loc[
        adata_sub.obs["annotation_MOFA"].isin(["pericytesSMCs", "fibroblast"]),
        "annotation_MOFA",
    ] = "mesenchymal"
    ctype_counter(adata_sub, ctype_count, ctype_count_patient)
    adata_sub.layers["counts"] = adata_sub.X.copy()
    # keep only relevant metadata columns
    cols_to_keep = list(set(adata_sub.obs.columns) & set(meta_columns))
    adata_sub.obs = adata_sub.obs[cols_to_keep]
    adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()
    adata_sub.write(output)

if study == "Li_2024":
    # control samples are almost all Cortex, therefore other regions will be excluded
    adata_sub = adata[
        (adata.obs["patient_category"].isin(["Health", "CKD"]))
        & (adata.obs["renal_region"] == "C")
    ]
    adata_sub, qc_df = qc(adata_sub, "Li_2024")
    qc_df.to_csv(qc_outputs)
    # identify cell types
    adata_sub.obs.loc[
        adata_sub.obs["celltype"].isin(["Fib1", "Fib2", "Fib3"]), "annotation_MOFA"
    ] = "fibroblast"
    adata_sub.obs.loc[adata_sub.obs["celltype"] == "ENDO", "annotation_MOFA"] = (
        "endothelial"
    )
    adata_sub.obs.loc[
        adata_sub.obs["celltype"].isin(["SMC1", "SMC2", "SMC3", "PC1", "PC2"]),
        "annotation_MOFA",
    ] = "pericytesSMCs"
    adata_sub.obs.loc[adata_sub.obs["celltype"] == "Ma", "annotation_MOFA"] = "myeloid"
    adata_sub.obs.loc[adata_sub.obs["celltype"] == "BT", "annotation_MOFA"] = "lymphoid"
    adata_sub.obs.loc[
        adata_sub.obs["celltype"].isin(
            [
                "CNT",
                "tL1",
                "DCT",
                "TAL2",
                "PT",
                "PT_dediff",
                "PT_VCAM1",
                "tL2",
                "TAL3",
                "ICA",
                "ICB",
                "POD",
                "JGA",
                "tL-TAL",
                "TAL1",
                "Uro1",
                "PEC",
                "Uro2",
            ]
        ),
        "annotation_MOFA",
    ] = "epithelial"
    # adapt nomenclature to make unify with other studies
    adata_sub.obs["study"] = "Li_2024"
    adata_sub.obs["modality"] = "shareseq"
    adata_sub.obs["cond_test"] = "control"
    adata_sub.obs.loc[adata_sub.obs["patient_category"] != "Health", "cond_test"] = (
        "fibrosis"
    )
    adata_sub.obs["grouping"] = adata_sub.obs["patient_category"]
    adata_sub.obs["grouping"] = adata_sub.obs["grouping"].replace(
        {"HKD": "H-CKD", "Health": "control"}
    )
    adata_sub.obs = adata_sub.obs.rename(
        columns={
            "celltype": "cell_type1",
            "patient_id": "sample",
            "renal_region": "region",
        }
    )
    adata_sub.obs["tech"] = "shareseq"
    adata_sub.obs["ctype"] = adata_sub.obs["annotation_MOFA"]
    adata_sub.obs.loc[
        adata_sub.obs["annotation_MOFA"].isin(["pericytesSMCs", "fibroblast"]),
        "annotation_MOFA",
    ] = "mesenchymal"
    ctype_counter(adata_sub, ctype_count, ctype_count_patient)
    adata_sub.layers["counts"] = adata_sub.X.copy()
    # keep only relevant metadata columns
    cols_to_keep = list(set(adata_sub.obs.columns) & set(meta_columns))
    adata_sub.obs = adata_sub.obs[cols_to_keep]
    adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()
    adata_sub.write(output)

if study == "Wilson_2022":
    adata_sub, qc_df = qc(adata, "Wilson_2022")
    del adata
    qc_df.to_csv(qc_outputs)
    add_metadata = pd.read_csv(metadata_file, sep="\t")
    keep = [
        "donor_id",
        "diabetes",
        "gfr_ckdepi",
        "serum_creatinine",
        "ACR_mg_per_g",
        "preop_dx",
        "hypertension",
        "global_glomerulosclerosis",
        "ifta",
    ]
    add_metadata = (
        add_metadata.filter(keep, axis=1)
        .set_index("donor_id")
        .rename({col: "clin_" + col for col in add_metadata.columns}, axis=1)
        .reset_index()
        .rename({"donor_id": "uid_donor_id"}, axis=1)
    )
    adata_sub.obs.loc[:, "donor_id"] = adata_sub.obs.loc[:, "donor_id"].str.replace(
        "_", ""
    )
    obs = pd.merge(
        adata_sub.obs.reset_index(),
        add_metadata,
        left_on="donor_id",
        right_on="uid_donor_id",
        how="left",
    ).set_index("index")
    obs[["condition_1"]] = None
    for ii, row in obs.iterrows():
        if "control" in row["uid_donor_id"] or "healthy" in row["uid_donor_id"]:
            obs.loc[ii, "condition_1"] = "control"
        elif "diabetic" in row["uid_donor_id"]:
            obs.loc[ii, "condition_1"] = "DKD"
    assert np.all(adata_sub.obs.index == obs.index)
    adata_sub.obs = obs
    # identify cell types
    adata_sub.obs.loc[
        adata_sub.obs["cell_type2_pred"].isin(
            ["MYOF", "M-FIB", "aFIB", "dFIB", "FIB", "dM-FIB"]
        ),
        "annotation_MOFA",
    ] = "fibroblast"
    adata_sub.obs.loc[
        adata_sub.obs["cell_type2_pred"].isin(
            ["EC-PTC", "EC-GC", "EC-DVR", "EC-AEA", "EC-LYM", "EC-AVR", "cycEC"]
        ),
        "annotation_MOFA",
    ] = "endothelial"
    adata_sub.obs.loc[
        adata_sub.obs["cell_type2_pred"].isin(["REN", "MC", "VSMC/P", "dVSMC", "VSMC"]),
        "annotation_MOFA",
    ] = "pericytesSMCs"
    adata_sub.obs.loc[
        adata_sub.obs["cell_type2_pred"].isin(
            ["MAC-M2", "N", "MDC", "cycMNP", "ncMON", "cDC", "MAST"]
        ),
        "annotation_MOFA",
    ] = "myeloid"
    adata_sub.obs.loc[
        adata_sub.obs["cell_type2_pred"].isin(["T", "B", "PL", "NKT"]),
        "annotation_MOFA",
    ] = "lymphoid"
    adata_sub.obs.loc[
        adata_sub.obs["cell_type2_pred"].isin(
            [
                "PT-S1",
                "C-TAL",
                "DCT2",
                "IC-B",
                "aPT",
                "CCD-IC-A",
                "CNT-PC",
                "DCT1",
                "PT-S3",
                "aTAL1",
                "dC-IC-A",
                "POD",
                "M-TAL",
                "DTL1",
                "dM-TAL",
                "CNT-IC-A",
                "DTL2",
                "dPOD",
                "dOMCD-PC",
                "dPT/DTL",
                "aTAL2",
                "OMCD-PC",
                "tPC-IC",
                "cycPT",
                "OMCD-IC-A",
                "IMCD",
                "dIMCD",
                "MD",
                "CNT",
                "CCD-PC",
                "PEC",
                "PT-S2",
                "dCNT",
                "dPT",
                "cycCNT",
            ]
        ),
        "annotation_MOFA",
    ] = "epithelial"
    adata_sub.obs["study"] = "Wilson_2022"
    adata_sub.obs["modality"] = "sn"
    adata_sub.obs["cond_test"] = "control"
    adata_sub.obs.loc[adata_sub.obs["condition_1"] != "control", "cond_test"] = (
        "fibrosis"
    )
    adata_sub.obs["grouping"] = adata_sub.obs["condition_1"]
    adata_sub.obs = adata_sub.obs.drop(columns=["cell_type2"])
    adata_sub.obs = adata_sub.obs.rename(
        columns={
            "label_before": "cell_type1",
            "cell_type2_pred": "cell_type2",
            "clin_sex": "sex",
            "sample_uuid": "sample",
            "clin_ifta": "IFTA",
            "clin_gfr_ckdepi": "eGFR",
        }
    )
    adata_sub.obs["cell_type2"] = adata_sub.obs["cell_type2"].astype("str")
    adata_sub.obs["tech"] = "10x 5' v2"
    adata_sub.obs["ctype"] = adata_sub.obs["annotation_MOFA"]
    adata_sub.obs.loc[
        adata_sub.obs["annotation_MOFA"].isin(["pericytesSMCs", "fibroblast"]),
        "annotation_MOFA",
    ] = "mesenchymal"
    ctype_counter(adata_sub, ctype_count, ctype_count_patient)
    adata_sub.layers["counts"] = adata_sub.X.copy()
    # keep only relevant metadata columns
    cols_to_keep = list(set(adata_sub.obs.columns) & set(meta_columns))
    adata_sub.obs = adata_sub.obs[cols_to_keep]
    adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()
    adata_sub.write(output)

# UMAP embedding
if study in [
    "Abedini_2024_sn",
    "Abedini_2024_sc",
    "Muto_2022",
    "Li_2024",
    "Wilson_2022",
    "Lake_2023_sc",
]:
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
    # Compute distances in the PCA space, and find cell neighbors
    sc.pp.neighbors(adata_sub, n_neighbors=25, n_pcs=30)
    # Generate UMAP features
    sc.tl.umap(adata_sub)
    adata_new = adata_sub[adata_sub.obs["annotation_MOFA"].notna()].copy()

# --- Save UMAP plots ---
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
