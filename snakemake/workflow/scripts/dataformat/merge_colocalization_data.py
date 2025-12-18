import pandas as pd
import os
import numpy as np

# snakemake inputs
data_path = snakemake.input["data_path"]

# snakemake params
real_names = snakemake.params["real_names"]

# snakemake outputs
out_path = snakemake.output["out_path"]

condition_map = {
    "PSC": "fibrosis",
    "CKD": "fibrosis",
    "FZ": "fibrosis",
    "IPF": "fibrosis",
    "CTRL": "control",
    "Healthy": "control",
    "normal": "control",
}

all_df = []

for file_path in data_path:
    cosine_sim = pd.read_csv(file_path, index_col=0)
    name = file_path.split("/")[-1].replace(".csv", "")
    cosine_sim["organ"] = real_names[name]
    cosine_sim["cond_test"] = cosine_sim["cond_test"].replace(condition_map)
    cosine_sim = cosine_sim[cosine_sim["cond_test"].isin(["fibrosis", "control"])]

    cosine_sim = cosine_sim[
        [
            "interaction",
            "ligand",
            "receptor",
            "morans",
            "mean",
            "std",
            "cond_test",
            "organ",
        ]
    ]
    cosine_sim = cosine_sim.rename(
        columns={"mean": "mean_cosine_similarity", "std": "std_cosine_similarity"}
    )

    all_df.append(cosine_sim)

merged = pd.concat(all_df)

merged.to_pickle(out_path)
