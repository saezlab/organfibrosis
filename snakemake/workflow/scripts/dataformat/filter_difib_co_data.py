import pandas as pd

# snakmake inputs
fibs_path = snakemake.input["data_path"]

# snakmake outputs
out_path = snakemake.output["out_path"]

fibs = pd.read_pickle(fibs_path)

fibs = fibs.reset_index().set_index("gene")
all_genes = fibs.groupby(["gene", "summary_row"]).count().reset_index()
heart_genes = all_genes[all_genes["summary_row"] == "heart"]["gene"]
lung_genes = all_genes[all_genes["summary_row"] == "lung"]["gene"]
kidney_genes = all_genes[all_genes["summary_row"] == "kidney"]["gene"]
liver_genes = all_genes[all_genes["summary_row"] == "liver"]["gene"]
fibs["organ_nr"] = fibs.groupby("gene").size() - 4
fibs.loc[heart_genes, "heart"] = fibs[fibs["summary_row"] == "heart"]["eff"][
    heart_genes
]
fibs.loc[lung_genes, "lung"] = fibs[fibs["summary_row"] == "lung"]["eff"][lung_genes]
fibs.loc[kidney_genes, "kidney"] = fibs[fibs["summary_row"] == "kidney"]["eff"][
    kidney_genes
]
fibs.loc[liver_genes, "liver"] = fibs[fibs["summary_row"] == "liver"]["eff"][
    liver_genes
]
# filter genes where model did not assign weights properly
genes_to_keep = fibs.groupby(fibs.index)["w_re"].apply(lambda x: (x.dropna() > 0).all())

working_model = fibs.loc[fibs["summary_row"] == "random effect",].loc[genes_to_keep, :]
working_model["celltype"] = "disease fibroblasts"
working_model = working_model[
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

working_model.to_pickle(out_path)
