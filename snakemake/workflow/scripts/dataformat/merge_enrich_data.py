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


enrich_dict = {
    ctype: {
        "FC": ulm_dict[ctype]["collecttri"],
        "pval": pval_dict[ctype]["collecttri"],
        "neglogpval": -np.log10(pval_dict[ctype]["collecttri"]),
    }
    for ctype in ulm_dict
}

all_df = []

for ctype in enrich_dict:
    df = enrich_dict[ctype]
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
    df_long["celltype"] = ctype
    all_df.append(df_long)

all_df = pd.concat(all_df)

all_df.to_pickle(out_path)
