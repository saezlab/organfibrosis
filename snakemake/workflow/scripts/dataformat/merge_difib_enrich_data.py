import pickle
import numpy as np
import pandas as pd

# snakemake inputs
pval_path = snakemake.input["pval_path"]
coef_path = snakemake.input["coef_path"]

# snakemake params
real_names = snakemake.params["real_names"]

# snakemake outputs
out_path = snakemake.output["out_path"]


with open(pval_path, "rb") as fp:
    pval_dict = pickle.load(fp)


with open(coef_path, "rb") as fp:
    ulm_dict = pickle.load(fp)

df = {
    "FC": ulm_dict["collecttri"],
    "pval": pval_dict["collecttri"],
    "neglogpval": -np.log10(pval_dict["collecttri"]),
}

fc_vals = df["FC"]
np_vals = df["neglogpval"]

# Reset index and clean organ names
df = fc_vals.reset_index().rename(columns={"index": "organ"})
df["organ"] = df["organ"].replace(real_names)
df_long = df.melt(id_vars="organ", var_name="TF", value_name="FC")

np_df = np_vals.reset_index().rename(columns={"index": "organ"})
np_df["organ"] = np_df["organ"].replace(real_names)
np_long = np_df.melt(id_vars="organ", var_name="TF", value_name="neglog10_padj")

# Merge the two long dataframes on organ and TF
df_long = pd.merge(df_long, np_long, on=["TF", "organ"])
df_long["celltype"] = "disease fibroblasts"

df_long.to_pickle(out_path)
