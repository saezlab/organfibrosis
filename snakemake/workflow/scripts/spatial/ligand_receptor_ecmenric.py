# Author: Leonie Küchenhoff
# Script compared diff. expressed genes of disease-assoc.
# fibroblasts to the spatial coexpression in scar areas.

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import numpy as np
import liana as li
import upsetplot
from collections.abc import Mapping
from pathlib import Path
from typing import Union

config = snakemake.config
geneset_inputs = snakemake.input["genesets"]

cosine_sim_input = snakemake.input.cosine_sim
if isinstance(cosine_sim_input, Mapping):
    cosine_sim_paths = dict(cosine_sim_input)
else:
    cosine_sim_paths = dict(zip(config["meta_organs"], cosine_sim_input))

cosine_results_output_path = Path(snakemake.output.cosine_summary)
scatter_org_plot_path = Path(snakemake.output.scatter_org_plot)
scatter_org_data_path = Path(snakemake.output.scatter_org_plot_csv)
upset_org_plot_path = Path(snakemake.output.upset_org_plot)
upset_org_data_path = Path(snakemake.output.upset_org_plot_csv)
scatter_crossorg_plot_path = Path(snakemake.output.scatter_crossorg_plot)
scatter_crossorg_data_path = Path(snakemake.output.scatter_crossorg_plot_csv)
upset_crossorg_plot_path = Path(snakemake.output.upset_crossorg_plot)
upset_crossorg_data_path = Path(snakemake.output.upset_crossorg_plot_csv)
top_genes_heatmap_plot_path = Path(snakemake.output.top_genes_heatmap)
top_genes_heatmap_data_path = Path(snakemake.output.top_genes_heatmap_csv)


co_path = snakemake.input.cross_organ_dl
organ_spec_up = snakemake.input.organ_spec_dl


plt.rcParams.update({"font.size": 20})

condition_colors = config["general_plotting"]["condition_colors"]
organs = config["meta_organs"]
all_studies = [config["datasets"][organ]["studies_ownmodel"] for organ in organs]
all_studies_for_cols = list(config["general_plotting"]["study_colors"].keys())
org_colors = config["general_plotting"]["organ_colors"]
org_colors_light = config["general_plotting"]["organ_colors_light"]
real_names = config["general_plotting"]["organ_names"]
study_colors = config["general_plotting"]["study_colors"]

org_color_real = {
    real_names_key: org_colors[real_names_dict_key]
    for real_names_dict_key, real_names_key in real_names.items()
}
resource = li.resource.select_resource(resource_name="consensus")
all_ligands = resource["ligand"].unique()
all_receptors = resource["receptor"].unique()


def load_gmt_to_dataframe(gmt_file_path):
    data = []  # To store each row of data
    with open(gmt_file_path, "r") as file:
        for line in file:
            parts = line.strip().split("\t")  # Splitting each line by tab
            gene_set_name = parts[0]
            description = parts[1]
            genes = list(parts[2:])  # The rest of the parts are genes
        # Convert the list of dictionaries to a DataFrame
        df = pd.DataFrame(genes, columns=["genesymbol"])
        df["geneset"] = gene_set_name
        df["description"] = description
    return df


def concatenate_gmt_files(gmt_files_list):
    dfs = []  # To store DataFrames loaded from each file
    for file_path in gmt_files_list:
        df = load_gmt_to_dataframe(file_path)
        dfs.append(df)
    # Concatenate all DataFrames into one
    concatenated_df = pd.concat(dfs, ignore_index=True)
    return concatenated_df


resource = li.resource.select_resource(resource_name="consensus")
ligands = resource["ligand"].unique()
receptors = resource["receptor"].unique()

ecm = concatenate_gmt_files(geneset_inputs)


ligands_filtered = set(ligands) - set(ecm["genesymbol"])
receptors_filtered = set(receptors) - set(ecm["genesymbol"])

all_filtered = ligands_filtered | receptors_filtered



def upset_to_dataframe(upset_data: Union[pd.Series, pd.DataFrame], value_col: str = "count") -> pd.DataFrame:
    """
    Normalize an upsetplot from_contents result to a flat dataframe.
    Handles both Series (older versions) and DataFrame (newer versions).
    """
    if isinstance(upset_data, pd.Series):
        return upset_data.to_frame(name=value_col).reset_index()
    df = upset_data.reset_index()
    non_index_cols = [col for col in df.columns if col not in upset_data.index.names]
    if len(non_index_cols) == 1 and non_index_cols[0] != value_col:
        df = df.rename(columns={non_index_cols[0]: value_col})
    return df


cosine_sim_dict = {}
for organ in organs:
    cosine_sim = pd.read_csv(cosine_sim_paths[organ], index_col=0)
    cosine_sim = cosine_sim[cosine_sim["receptor"] == 'NABA_CORE_MATRISOME']
    cosine_sim["organ"] = organ
    cosine_sim["cond_test"] = cosine_sim["cond_test"].replace(
        {
            "PSC": "fibrosis",
            "CKD": "fibrosis",
            "FZ": "fibrosis",
            "IPF": "fibrosis",
            "CTRL": "control",
            "Healthy": "control",
            "normal": "control",
        }
    )
    cosine_sim = cosine_sim[cosine_sim["cond_test"].isin(["fibrosis", "control"])]
    cosine_sim_dict[organ] = cosine_sim

all_cosine_results = pd.concat([j for i, j in cosine_sim_dict.items()])
all_cosine_results.to_csv(cosine_results_output_path)

with open(co_path, "rb") as fp:
    co = pickle.load(fp)

with open(organ_spec_up, "rb") as fp:
    up = pickle.load(fp)


high_both_orgspec = {}
data_allorgan_dict = {}


fig, axs = plt.subplots(2, 2, figsize=(9, 8), tight_layout=True)
ax = axs.ravel()
scatter_org_rows = []
for count, organ in enumerate(organs):
    organ_real = real_names[organ]
    top = 1000
    fibs = co.reset_index().set_index("gene")
    lig_rec_fibs = fibs[fibs.index.isin(list(all_filtered))]
    top_hits_up = (
        lig_rec_fibs.loc[lig_rec_fibs["summary_row"] == organ_real,]
        .sort_values(by="ci_low")
        .iloc[-top:]
        .index
    )

    print(top_hits_up[:5])
    print(organ)

    only_disease =  cosine_sim_dict[organ][cosine_sim_dict[organ]['cond_test'] == 'fibrosis']
    mean_morans = only_disease.groupby('ligand')['morans'].mean()
    print(lig_rec_fibs.head())
    print(mean_morans.head())

    to_plot = (
        lig_rec_fibs.loc[lig_rec_fibs["summary_row"] == organ_real,]
        .loc[top_hits_up]
        .merge(mean_morans, left_index=True, right_index=True, how = 'outer')
        #.merge(mean_morans, left_index=True, right_index=True)
    )
    print(to_plot.head())

    x_thresh =  np.percentile(to_plot[~to_plot['eff'].isna()]['eff'], 80)
    y_thresh =  np.percentile(to_plot[~to_plot['morans'].isna()]['morans'], 80)


    to_plot["expression"] = to_plot["eff"].rank(pct=True)
    to_plot["spatial"] = to_plot["morans"].rank(pct=True)

    # Make a flag column
    to_plot["highlight"] = (to_plot["eff"] > x_thresh) & (to_plot["morans"] > y_thresh)
    data_allorgan_dict[organ] = to_plot

    high_both_orgspec[organ] = list(to_plot[to_plot["highlight"]].index)

    clean_for_corr = to_plot[["eff", "morans"]].dropna()
    corr_val = clean_for_corr["eff"].corr(clean_for_corr["morans"]) if not clean_for_corr.empty else float("nan")
    if pd.notna(corr_val):
        print(f"[ligand_receptor] Organ-specific scatter {organ_real}: R = {corr_val:.3f}")
    else:
        print(f"[ligand_receptor] Organ-specific scatter {organ_real}: R could not be computed (insufficient data).")

    scatter_org_rows.append(
        to_plot.reset_index()
        .rename(columns={"index": "gene"})
        .assign(organ=organ, organ_real=organ_real, plot_type="organ_specific")
    )

    sns.scatterplot(
        data=to_plot,
        x="eff",
        y="morans",
        hue="highlight",
        palette={False: "darkgrey", True: org_color_real[organ_real]},
        ax=ax[count],
    )
    ax[count].set_ylabel(f"avg. morans I in {organ_real}")
    ax[count].set_xlabel(f"effect size in {organ_real}")
    ax[count].get_legend().remove()
    ax[count].set_xlim(-5, 5)
scatter_org_plot_path.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(scatter_org_plot_path, bbox_inches="tight")
if scatter_org_rows:
    scatter_org_df = pd.concat(scatter_org_rows, ignore_index=True)
else:
    scatter_org_df = pd.DataFrame()
scatter_org_df.to_csv(scatter_org_data_path, index=False)


upset_df = upsetplot.from_contents(high_both_orgspec)
upset_df.index.names = [real_names[i] for i in upset_df.index.names]

upsetplot.plot(upset_df, show_counts=True, sort_categories_by="input", element_size = 28)
plt.savefig(upset_org_plot_path)
upset_org_data = upset_to_dataframe(upset_df)
upset_org_data.to_csv(upset_org_data_path, index=False)


high_both = {}
data_all_dict = {}
scatter_crossorg_rows = []
fig, axs = plt.subplots(2, 2, figsize=(9, 8), sharex=True, tight_layout=True)
ax = axs.ravel()
for count, organ in enumerate(organs):
    print(organ)
    organ_real = real_names[organ]
    top = 1000
    fibs = co.reset_index().set_index("gene")
    lig_rec_fibs = fibs[fibs.index.isin(list(all_filtered))]
    only_disease = cosine_sim_dict[organ][
        cosine_sim_dict[organ]["cond_test"] == "fibrosis"
    ]
    mean_morans = only_disease.groupby("ligand")["morans"].mean()
    to_plot = lig_rec_fibs.loc[lig_rec_fibs["summary_row"] == "random effect",].merge(
        mean_morans, left_index=True, right_index=True, how="outer"
    )

    x_thresh =  np.percentile(to_plot[~to_plot['eff'].isna()]['eff'], 80)
    to_plot["x_percentile"] = to_plot["eff"].rank(pct=True)

    to_plot = to_plot.drop(columns="w_fe").dropna(axis=0)

    y_thresh =  np.percentile(to_plot[~to_plot['morans'].isna()]['morans'], 80)
    to_plot["y_percentile"] = to_plot["morans"].rank(pct=True)

    # Make a flag column
    to_plot["highlight"] = (to_plot["eff"] > x_thresh) & (to_plot["morans"] > y_thresh)

    high_both[organ] = list(to_plot[to_plot["highlight"]].index)

    data_all_dict[organ] = to_plot

    clean_for_corr = to_plot[["eff", "morans"]].dropna()
    corr_val = clean_for_corr["eff"].corr(clean_for_corr["morans"]) if not clean_for_corr.empty else float("nan")
    if pd.notna(corr_val):
        print(f"[ligand_receptor] Cross-organ scatter {organ_real}: R = {corr_val:.3f}")
    else:
        print(f"[ligand_receptor] Cross-organ scatter {organ_real}: R could not be computed (insufficient data).")

    scatter_crossorg_rows.append(
        to_plot.reset_index()
        .rename(columns={"index": "gene"})
        .assign(organ=organ, organ_real=organ_real, plot_type="cross_organ")
    )

    sns.scatterplot(
        data=to_plot,
        x="eff",
        y="morans",
        hue="highlight",
        palette={False: "darkgrey", True: org_color_real[organ_real]},
        ax=ax[count],
    )  # {True:org_color_real[organ_real], False:org_colors_light[organ]}, edgecolor= 'grey')
    ax[count].set_ylabel(f"avg. morans I in {organ_real}")
    ax[count].set_xlabel("combined organ effect size")
    ax[count].get_legend().remove()
scatter_crossorg_plot_path.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(scatter_crossorg_plot_path, bbox_inches="tight")
if scatter_crossorg_rows:
    scatter_crossorg_df = pd.concat(scatter_crossorg_rows, ignore_index=True)
else:
    scatter_crossorg_df = pd.DataFrame()
scatter_crossorg_df.to_csv(scatter_crossorg_data_path, index=False)


upset_df = upsetplot.from_contents(high_both)

upsetplot.plot(upset_df, show_counts=True, sort_categories_by="input", element_size = 28)
plt.savefig(upset_crossorg_plot_path)
upset_crossorg_data = upset_to_dataframe(upset_df)
upset_crossorg_data.to_csv(upset_crossorg_data_path, index=False)

merged = pd.concat(
    {
        k: v[["spatial", "expression"]].add_suffix(f" {real_names[k]}")
        for k, v in data_allorgan_dict.items()
    },
    axis=1,
)

merged.columns = merged.columns.droplevel(0)
covered = merged.dropna(thresh=6)
top_genes = covered.mean(axis=1).sort_values().iloc[-40:].index
merged["mean"] = merged.mean(axis=1)
org_colors_plot = [org_colors[i] for i in organs for _ in range(2)] + ["white"]
type_colors_plot = ["grey", "black"] * 4 + ["white"]


heatmap_gene_order = list(reversed(top_genes))
heatmap_plot_data = merged.loc[heatmap_gene_order]
fig = sns.clustermap(
    heatmap_plot_data,
    vmin=0,
    vmax=1,
    figsize=(7, 20),
    cmap="plasma",
    cbar_kws={"label": "percentile"},
    row_cluster=False,
    col_cluster=False,
    col_colors=[org_colors_plot, type_colors_plot],
    cbar_pos=(0, 0.45, 0.03, 0.1),
    label="percentile",
    yticklabels=True
)
plt.savefig(top_genes_heatmap_plot_path)
heatmap_plot_data.index.name = "interaction"
heatmap_plot_data.to_csv(top_genes_heatmap_data_path)
