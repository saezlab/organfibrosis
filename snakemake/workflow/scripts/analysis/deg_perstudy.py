# Author: Leonie Küchenhoff
# This script performs differential gene expression analysis for each cell type across multiple studies,
# and saves the results in a CSV file.
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import scanpy as sc
import decoupler as dc

# Get input and output paths from snakemake
pb_path = snakemake.input[0]
meta_path = snakemake.input[1]
out_path = snakemake.output[0]
filter_by_expr = snakemake.params[0].get("filter_by_expr")
filter_by_prop = snakemake.params[0].get("filter_by_prop")


metadata = pd.read_csv(meta_path)




pdata = sc.read_csv(pb_path)

index_names = pdata.obs.index.str.split("_")[:].tolist()

pdata.obs["ctype"] = [item[-1] for item in index_names]
pdata.obs["sample"] = [
    "_".join(item[:-1]) if len(item) > 1 else item[1] for item in index_names
]
pdata.obs = pdata.obs.merge(metadata, on="sample", how="left")

pdata.obs.index = pdata.obs["sample"]

results_df = pd.DataFrame()
for count, ctype in enumerate(pdata.obs["ctype"].unique()):
    # check which cell types are still represented in data
    ctype_pdata = pdata[pdata.obs["ctype"] == ctype].copy()

    # edgeR filtering
    if filter_by_expr is not None:
        genes = dc.filter_by_expr(ctype_pdata, group="cond_test", **filter_by_expr)
        ctype_pdata = ctype_pdata[:, genes].copy()

    # proportion filtering
    if filter_by_prop is not None:
        genes = dc.filter_by_prop(ctype_pdata, group="cond_test", **filter_by_prop)
        ctype_pdata = ctype_pdata[:, genes].copy()

    # only do dge analysis when there are at least 3 disease and 3 reference samples:
    if (ctype_pdata[ctype_pdata.obs["cond_test"] == "fibrosis"].shape[0] >= 3) & (
        ctype_pdata[ctype_pdata.obs["cond_test"] == "control"].shape[0] >= 3
    ):
        print(ctype_pdata.obs.head())
        ctype_pdata
        dds = DeseqDataSet(
            adata=ctype_pdata,
            design_factors="cond_test",
            ref_level=["cond_test", "control"],
            refit_cooks=True,
        )
        # Compute LFCs
        dds.deseq2()
        # Extract contrast between condition vs normal
        stat_res = DeseqStats(dds, contrast=["cond_test", "fibrosis", "control"])
        # Shrink LFCs
        stat_res.summary()
        stat_res.lfc_shrink(coeff=f"cond_test[T.fibrosis]")
        # Extract results
        results_df_new = stat_res.results_df

        results_df_new.columns = results_df_new.columns + f"_{ctype}"
        results_df = results_df.merge(
            results_df_new, left_index=True, right_index=True, how="outer"
        )


results_df.to_csv(out_path)
