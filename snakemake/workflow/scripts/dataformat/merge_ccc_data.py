import pandas as pd
import numpy as np

# snakmake inputs
ccc_results_path = snakemake.input["data_path"]

# snakemake params
organs = snakemake.params["organs"]
celltypes = snakemake.params["celltypes"]
organ_mapping = snakemake.params["organ_mapping"]

# snakmake outputs
out_path = snakemake.output["out_path"]

ccc_results = pd.read_pickle(ccc_results_path)

eff_cutoff_pos = 0.5
eff_cutoff_neg = -0.5

counts_dict_s = {}
counts_dict_t = {}
for organ in organs:
    ccc_result_organ = ccc_results[organ]
    sig = ccc_result_organ[
        (ccc_result_organ["interaction_eff"] > eff_cutoff_pos)
        & (
            ccc_result_organ["interaction_eff"] - ccc_result_organ["interaction_sd_eff"]
            > 0
        )
    ]
    heatmap_target = (
        sig.groupby(["interaction", "target"])
        .count()["ligand"]
        .reset_index()
        .pivot(columns="target", values="ligand", index="interaction")
    )
    heatmap_source = (
        sig.groupby(["interaction", "source"])
        .count()["ligand"]
        .reset_index()
        .pivot(columns="source", values="ligand", index="interaction")
    )

    counts_dict_s[organ] = heatmap_source
    counts_dict_t[organ] = heatmap_target

# Ensure all DataFrames have the same structure
dfs_s = [counts_dict_s[organ] for organ in organs]
dfs_t = [counts_dict_t[organ] for organ in organs]

# Outer join to align indices across all DataFrames
merged_df = pd.concat(dfs_s, axis=0, join="outer").groupby(level=0).sum(min_count=1)

# Fill NaN values with 0 to ensure proper summation
sum_df = (
    pd.concat([df.fillna(0) for df in dfs_s], axis=0, join="outer")
    .groupby(level=0)
    .sum()
)
sum_df = sum_df[celltypes]
sum_df = sum_df.add_prefix("source_")

sum_df_t = (
    pd.concat([df.fillna(0) for df in dfs_t], axis=0, join="outer")
    .groupby(level=0)
    .sum()
)
sum_df_t = sum_df_t[celltypes]
sum_df_t = sum_df_t.add_prefix("target_")

final_df = pd.DataFrame()
final_df[["ligand", "receptor"]] = sum_df.index.to_series().str.split("^", expand=True)
final_df = pd.concat([final_df, sum_df, sum_df_t], axis=1)

for org in organs:
    name = "n_" + organ_mapping[org]
    final_df[name] = counts_dict_s[org].sum(axis=1)

final_df = final_df.fillna(0)

final_df.to_pickle(out_path)
