import pickle

import pandas as pd

# Snakemake inputs/outputs/params
organ_spec_dl_path = snakemake.input["organ_spec_dl"]
organs = snakemake.params["organs"]
ctypes = snakemake.params["views"]
organ_names = snakemake.params["organ_names"]
top_genes_csv = snakemake.output["top_genes_csv"]

with open(organ_spec_dl_path, "rb") as fp:
    co = pickle.load(fp)

celltypes = list(co.keys())
final_fibs = {}
top_upregulated_rows = []

for celltype in ctypes:
    print(celltype)
    for organ in organs:
        df = co[celltype][organ].reset_index().set_index("gene")
        all_genes = df.groupby(["gene", "summary_row"]).count().reset_index()
        heart_genes = all_genes[all_genes["summary_row"] == "heart"]["gene"]
        lung_genes = all_genes[all_genes["summary_row"] == "lung"]["gene"]
        kidney_genes = all_genes[all_genes["summary_row"] == "kidney"]["gene"]
        liver_genes = all_genes[all_genes["summary_row"] == "liver"]["gene"]
        df["organ_nr"] = df.groupby("gene").size() - 4
        # filter genes where model did not assign weights properly
        df = df.loc[df["summary_row"] == "random effect", ].sort_values(by="gene")
        print(organ)

        filtered = df[(df["organ_nr"] >= 3) & (df["ci_upp"] >= 0)]
        top_genes = filtered.sort_values(by="eff", ascending=False)
        print(list(top_genes.index[:3]))

        top_selection = top_genes.head(10).reset_index()[["gene"]]
        top_selection["celltype"] = celltype
        top_selection["organ"] = organ_names.get(organ, organ)
        top_selection["rank"] = list(range(1, len(top_selection) + 1))
        top_upregulated_rows.append(top_selection)

if top_upregulated_rows:
    top_10_upregulated = pd.concat(top_upregulated_rows, ignore_index=True)
else:
    top_10_upregulated = pd.DataFrame()

top_10_upregulated.to_csv(top_genes_csv, index=False)
