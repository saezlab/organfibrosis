"""Build full and spatial patient-scoring gene sets from organ effects and GMT files.

Converted from local_test_scripts/genesets_patient_scores.ipynb. Selection,
weighting, ordering, and duplicate handling follow the notebook.
"""

import pickle

import pandas as pd


def load_gmt_to_dataframe(gmt_file_path):
    with open(gmt_file_path, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            gene_set_name = parts[0]
            description = parts[1]
            genes = list(parts[2:])
    df = pd.DataFrame(genes, columns=["genesymbol"])
    df["geneset"] = gene_set_name
    df["description"] = description
    return df

def concatenate_gmt_files(gmt_files_list):
    dfs = []
    for file_path in gmt_files_list:
        dfs.append(load_gmt_to_dataframe(file_path))
    return pd.concat(dfs, ignore_index=True)

genesets = concatenate_gmt_files([snakemake.input.inflammatory, snakemake.input.core_matrisome])

genesets.loc[:, 'collection'] = 'lit_genesets'

genesets.loc[:, 'weight'] = '1'

with open(snakemake.input.organ_spec_dl, "rb") as fp:
    co = pickle.load(fp)

top_upregulated_rows = []
top_downregulated_rows = []

for celltype in snakemake.params.celltypes:
    for organ in snakemake.params.organs:
        df = co[celltype][organ].reset_index().set_index("gene")

        df["organ_nr"] = df.groupby("gene").size() - snakemake.params.summary_rows
        # filter genes where model did not assign weights properly
        df = df.loc[df["summary_row"] == "random effect", ].sort_values(by="gene")
        # get top genes
        filtered = df[(df["organ_nr"] >= snakemake.params.min_studies) & (df["ci_upp"] >= snakemake.params.ci_threshold)]
        top_genes = filtered.sort_values(by="eff", ascending=False)

        top_selection = top_genes.head(snakemake.params.top_up_per_celltype).reset_index()[["gene", "eff", 'sd_eff']]
        top_selection["celltype"] = celltype
        top_selection["organ"] = snakemake.params.organ_names.get(organ, organ)
        top_selection["rank"] = list(range(1, len(top_selection) + 1))
        top_selection['direction'] = 'upregulated'
        top_upregulated_rows.append(top_selection)

        # also get bottom genes
        filtered = df[(df["organ_nr"] >= snakemake.params.min_studies) & (df["ci_upp"] <= snakemake.params.ci_threshold)]
        bottom_genes = filtered.sort_values(by="eff", ascending=True)
        bottom_selection = bottom_genes.head(snakemake.params.top_down_per_celltype).reset_index()[["gene", "eff", 'sd_eff']]
        bottom_selection["celltype"] = celltype
        bottom_selection["organ"] = snakemake.params.organ_names.get(organ, organ)
        bottom_selection["rank"] = list(range(1, len(bottom_selection) + 1))
        bottom_selection['direction'] = 'downregulated'
        top_downregulated_rows.append(bottom_selection)

if top_upregulated_rows:
    top_10_upregulated = pd.concat(top_upregulated_rows, ignore_index=True)
else:
    top_10_upregulated = pd.DataFrame()

if top_downregulated_rows:
    top_10_downregulated = pd.concat(top_downregulated_rows, ignore_index=True)
else:
    top_10_downregulated = pd.DataFrame()

full_list = pd.concat([top_10_upregulated, top_10_downregulated], ignore_index=True)

out = (
    full_list.groupby("organ", group_keys=False)
      .apply(lambda x: pd.concat([
          x.nlargest(snakemake.params.top_up_per_organ, "eff"),
          x.nsmallest(snakemake.params.top_down_per_organ, "eff")
      ]))
      .drop_duplicates(subset=["organ", "gene"])
      .rename(columns={
          "organ": "geneset",
          "gene": "genesymbol",
          "eff": "weight"
      })
      [["geneset", "genesymbol", "weight"]]
)

out["collection"] = "organ_score"

out = out[[ "genesymbol","geneset", "collection", "weight"]]

full_set = pd.concat([genesets, out])

full_set.to_csv(snakemake.output.patient_scoring_genesets)

# for spatial enrichments
pro_col = concatenate_gmt_files([snakemake.input.proteoglycans, snakemake.input.collagens])

full_ecm = full_set[full_set['geneset']==snakemake.params.core_matrisome_name]
slim_ecm = full_ecm[full_ecm['genesymbol'].isin(pro_col['genesymbol'].unique())]

new_slim_set = pd.concat([full_set[full_set['geneset']!=snakemake.params.core_matrisome_name], slim_ecm])

new_slim_set.to_csv(snakemake.output.spatial_genesets)
