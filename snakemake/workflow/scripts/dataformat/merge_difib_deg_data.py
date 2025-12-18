import pandas as pd
import os
import numpy as np

# snakmake inputs
deg_file_list = snakemake.input["deg_file_list"]

# snakemake params
real_names = snakemake.params["real_names"]

# snakmake outputs
out_path = snakemake.output["out_path"]

rows = []


for file_path in deg_file_list:
    name = file_path.split("/")[3]
    study = file_path.split("/")[4].replace(".csv", "")
    df = pd.read_csv(file_path, index_col=0)
    lfc_col = "log2FoldChange"
    padj_col = "padj"
    if lfc_col in df.columns and padj_col in df.columns:
        temp = df[[lfc_col, padj_col]].copy()
        temp = temp.rename(columns={lfc_col: "log2FC", padj_col: "padj"})
        temp["gene"] = temp.index
        temp["organ"] = real_names[name]
        temp["study"] = study
        temp["celltype"] = "disease fibroblasts"
        rows.append(temp.reset_index(drop=True))

# Combine and save
combined_df = pd.concat(rows, ignore_index=True)
combined_df["padj"] = combined_df["padj"].replace(0, np.nan)
combined_df["neglog10_padj"] = -np.log10(combined_df["padj"])
mask = combined_df["log2FC"].notna() & combined_df["neglog10_padj"].isna()
combined_df.loc[mask, "neglog10_padj"] = 0.0

combined_df.replace([np.inf, -np.inf], np.nan, inplace=True)

combined_df.to_parquet(out_path)
