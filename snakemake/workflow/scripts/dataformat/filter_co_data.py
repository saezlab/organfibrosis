import pandas as pd
import pickle

# snakmake inputs
co_path = snakemake.input["data_path"]

# snakmake outputs
out_path = snakemake.output["out_path"]

with open(co_path, "rb") as fp:
    co = pickle.load(fp)

celltypes = co.keys()
final_fibs = {}

for ctype in celltypes:
    df = co[ctype].reset_index().set_index("gene")
    all_genes = df.groupby(["gene", "summary_row"]).count().reset_index()
    heart_genes = all_genes[all_genes["summary_row"] == "heart"]["gene"]
    lung_genes = all_genes[all_genes["summary_row"] == "lung"]["gene"]
    kidney_genes = all_genes[all_genes["summary_row"] == "kidney"]["gene"]
    liver_genes = all_genes[all_genes["summary_row"] == "liver"]["gene"]
    df["organ_nr"] = df.groupby("gene").size() - 4
    df.loc[heart_genes, "heart"] = df[df["summary_row"] == "heart"]["eff"][heart_genes]
    df.loc[lung_genes, "lung"] = df[df["summary_row"] == "lung"]["eff"][lung_genes]
    df.loc[kidney_genes, "kidney"] = df[df["summary_row"] == "kidney"]["eff"][
        kidney_genes
    ]
    df.loc[liver_genes, "liver"] = df[df["summary_row"] == "liver"]["eff"][liver_genes]
    # filter genes where model did not assign weights properly
    df = df.loc[df["summary_row"] == "random effect",].sort_values(by="gene")
    df["celltype"] = ctype
    final_fibs[ctype] = df[
        [
            "celltype",
            "eff",
            "sd_eff",
            "ci_low",
            "ci_upp",
            "heart",
            "lung",
            "kidney",
            "liver",
        ]
    ]

with open(out_path, "wb") as file:
    pickle.dump(final_fibs, file)
