# Script for pseudobulk differential expression analysis of disease fibroblasts across organs
#
# This script reads a single-cell AnnData object, annotates disease-associated fibroblasts,
# aggregates pseudobulk counts per sample and substate, and performs DE analysis using DESeq2.
# Results are saved as CSV files for each study. Filtering options and harmonization results
# are supported. The script is designed for use in a Snakemake workflow.

import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# snakmake inputs
fibroblast_path = snakemake.input["adata_path"]
scvi_path = snakemake.input["scvi_path"]

# snakmake outputs
out_paths = snakemake.output["deg_paths"]
anno_path = snakemake.output["anno_path"]


# snakmake parameters
filter_by_expr = snakemake.params["filters"].get("filter_by_expr")
filter_by_prop = snakemake.params["filters"].get("filter_by_prop")
min_cells = snakemake.params["filters"].get("min_cells")
min_counts = snakemake.params["filters"].get("min_counts")
resolutions = snakemake.params["resolutions"]
diseasefibs = snakemake.params["diseasefibs"]
organs = snakemake.params["organs"]


def read_and_add_harmon_results(adata_to_add, path, suffix, id_study=False, pcs=50):
    """
    function to read in PCA, UMAP, and leiden clustering results and
    add them to adata object with the suffix so that several results can be added
    input:
    adata_to_add - adata object where results will be added
    path - path at wich PCA, UMAP and leiden results are saved
    suffix - suffix is added to names, i.e. UMPA results is saved at X_umap_suffix
    id_study - if id column is just the cell id or a merged cell id + study
    pcs - no. of pcs

    output:
    adata object with adde results
    """
    harm_results = pd.read_csv(
        path,
        index_col=0,
        dtype = {'id':'str','leiden_0.2':str, 
                        'leiden_0.3':str,'leiden_0.4':str,
                        'leiden_0.5':str,'leiden_0.7':str},
    )

    if id_study == False:
        assert np.all(adata_to_add.obs_names == harm_results["id"])
    elif id_study == True:
        harm_results = harm_results.set_index("Unnamed: 1")
        harm_results = harm_results.loc[
            (np.array(adata_to_add.obs.index)
              + '_' 
              + np.array(adata_to_add.obs['study']))
              ]
        print("success")
    else:
        assert np.all(
            (
                np.array(adata_to_add.obs.index)
                + "_"
                + np.array(adata_to_add.obs["study"])
            )
            == harm_results["id"]
        )
        harm_results = harm_results.set_index("id", drop=False)
        print("success")

    adata_to_add.obsm[f"X_harmony_{suffix}"] = np.array(harm_results.iloc[:, :pcs])
    adata_to_add.obsm[f"X_umap_{suffix}"] = np.array(
        harm_results.iloc[:, pcs + 1 : pcs + 3]
    )

    adata_to_add.obs[[f'leiden_0.2_{suffix}',
                    f'leiden_0.3_{suffix}',
                    f'leiden_0.4_{suffix}',
                    f'leiden_0.5_{suffix}',
                    f'leiden_0.7_{suffix}']] = harm_results[['leiden_0.2',
                                                            'leiden_0.3',
                                                            'leiden_0.4',
                                                            'leiden_0.5',
                                                            'leiden_0.7']].to_numpy()

    return adata_to_add


# read anndata and add harmonization results
adata = sc.read(fibroblast_path)
adata = read_and_add_harmon_results(adata, scvi_path, "scvi", id_study=True, pcs=30)

# annotate fibroblasts
adata.obs["substate"] = "rest"
for organ in organs:
    resolution = resolutions[organ]
    clusters = diseasefibs[organ]
    adata.obs.loc[
        (adata.obs["organ"] == organ)
        & (adata.obs[f"leiden_{resolution}_scvi"].isin(clusters)),
        "substate",
    ] = "diseasefib"


adata.obs["sample_study"] = (
    adata.obs["sample"].astype("str") + "_" + adata.obs["study"].astype("str")
)
# save annotation as csv
adata.obs.loc[:, ["sample", "study", "substate"]].to_csv(anno_path)

# OCEAN data was SoupX correted, following steps cannot all handle floats
adata.layers["counts"] = np.round(adata.layers["counts"])
# get pseudobulks
pdata = dc.get_pseudobulk(
    adata,
    sample_col="sample_study",
    groups_col="substate",
    layer="counts",
    mode="sum",
    min_cells=min_cells,
    min_counts=min_counts,
)


# for each study, calculate DEGs and save as csv

for path in out_paths:
    results_df = pd.DataFrame()
    study = path.split("/")[-1].split(".")[0]
    # check which cell types are still represented in data
    ctype_pdata = pdata[pdata.obs["study"] == study].copy()
    # Ensure grouping column has proper categorical type
    ctype_pdata.obs["substate"] = ctype_pdata.obs["substate"].astype(
        pd.CategoricalDtype(categories=["rest", "diseasefib"], ordered=True)
    )

    # edgeR filtering
    if filter_by_expr is not None:
        genes = dc.filter_by_expr(ctype_pdata, group="substate", **filter_by_expr)
        ctype_pdata = ctype_pdata[:, genes].copy()

    # proportion filtering
    if filter_by_prop is not None:
        genes = dc.filter_by_prop(ctype_pdata, group="substate", **filter_by_prop)
        ctype_pdata = ctype_pdata[:, genes].copy()

    # only do dge analysis when there are at least 3 diseasefib and 3 rest samples:
    if (ctype_pdata[ctype_pdata.obs["substate"] == "diseasefib"].shape[0] >= 3) & (
        ctype_pdata[ctype_pdata.obs["substate"] == "rest"].shape[0] >= 3
    ):
        print(ctype_pdata.obs.head())
        ctype_pdata
        dds = DeseqDataSet(
            adata=ctype_pdata,
            design_factors="substate",
            ref_level=["substate", "rest"],
            refit_cooks=True,
        )
        # Compute LFCs
        dds.deseq2()
        # Extract contrast between condition vs normal
        stat_res = DeseqStats(dds, contrast=["substate", "diseasefib", "rest"])
        # Shrink LFCs
        stat_res.summary()
        stat_res.lfc_shrink(coeff=f"substate[T.diseasefib]")
        # Extract results
        results_df = stat_res.results_df

    else:
        results_df = pd.DataFrame(
            columns=["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
        )
        print(f"{study} did not have enough samples")

    results_df.to_csv(path)
    print(f"{study} results saved at {path}.")


print("done.")
