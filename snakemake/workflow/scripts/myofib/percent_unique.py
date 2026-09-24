# Calculate and plot percentage of unique marker genes for "all" and "myofib" scenarios
import pickle

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import upsetplot

# Get file paths and parameters from Snakemake

organ_spec_dl_all = snakemake.input["organ_spec_dl_all"]
organ_spec_dl_myofib = snakemake.input["organ_spec_dl_myofib"]
output_pdf_comparison = snakemake.output["output_pdf_comparison"]
output_pdf_all_mesench = snakemake.output["plot_pdf_all_mesench"]
organs = snakemake.params["organs"]
all_studies = snakemake.params["all_studies"]
org_colors = snakemake.params["org_colors"]
real_names = snakemake.params["real_names"]
ctypes = snakemake.params["ctypes"]

org_color_real = {
    real_names_key: org_colors[real_names_dict_key]
    for real_names_dict_key, real_names_key in real_names.items()
}

# Plot parameters
plt.rcParams.update({'font.size': 18})

def get_genes(df, ctype, studies_organ, organ, gene_dict_up, gene_dict_down):
    fibs = df.reset_index().set_index("gene")
    all_genes = fibs.groupby(["gene", "summary_row"]).count().reset_index()
    for study in studies_organ:
        genes = all_genes[all_genes["summary_row"] == study]["gene"]
        fibs.loc[genes, f"{study}"] = fibs[fibs["summary_row"] == f"{study}"]["eff"][
            genes
        ]
    fibs["study_nr"] = fibs.groupby("gene").size() - 4
    genes_to_keep = fibs.groupby(fibs.index)["w_re"].apply(
        lambda x: (x.dropna() > 0).all()
    )
    working_model = fibs.loc[genes_to_keep]
    sig_results = working_model[
        (working_model["summary_row"] == "random effect")
        & (working_model["eff"] > 0.5)
        & (working_model["sd_eff"].notna())
        & (working_model["ci_low"] > 0)
    ]
    df_sub = sig_results.sort_values(by="gene")
    if ctype in gene_dict_up.keys():
        gene_dict_up[ctype][organ] = list(df_sub.index)
    else:
        gene_dict_up[ctype] = {}
        gene_dict_up[ctype][organ] = list(df_sub.index)
    sig_results = working_model[
        (working_model["summary_row"] == "random effect")
        & (working_model["eff"] < -0.5)
        & (working_model["sd_eff"].notna())
        & (working_model["ci_upp"] < 0)
    ]
    df_sub = sig_results.sort_values(by="gene")
    if ctype in gene_dict_down.keys():
        gene_dict_down[ctype][organ] = list(df_sub.index)
    else:
        gene_dict_down[ctype] = {}
        gene_dict_down[ctype][organ] = list(df_sub.index)


def get_percentage_unique(upset_df, organs, real_names):
    df = upset_df.reset_index().set_index("id")
    percentages = {}
    for organ in organs:
        is_true = df[organ]
        others = df.drop(columns=organ)
        only_this_true = is_true & others.eq(False).all(axis=1)
        total_true = is_true.sum()
        unique_true = only_this_true.sum()
        percentage = (unique_true / total_true * 100) if total_true > 0 else 0
        percentages[real_names[organ]] = round(percentage, 2)
    return percentages

def gene_dict_to_long_df(gene_dict, direction, real_names):
    rows = []

    for ctype, organ_dict in gene_dict.items():
        for organ, genes in organ_dict.items():
            for gene in genes:
                rows.append(
                    {
                        "ctype": ctype,
                        "organ": real_names[organ],
                        "direction": direction,
                        "gene": gene,
                    }
                )

    return pd.DataFrame(rows)


# Run both scenarios: 'all' and 'myofib'
results_scenarios = ["all", "myofib"]
percentage_unique_results = {}
gene_list_results = []
for results in results_scenarios:
    gene_dict_up = {}
    gene_dict_down = {}
    if results == "all":
        with open(organ_spec_dl_all, "rb") as fp:
            organ_spec = pickle.load(fp)
        for studies_organ, organ in zip(all_studies, organs):
            for ctype in ctypes:
                get_genes(
                    organ_spec[ctype][organ],
                    ctype,
                    studies_organ,
                    organ,
                    gene_dict_up,
                    gene_dict_down,
                )
    elif results == "myofib":
        with open(organ_spec_dl_myofib, "rb") as fp:
            organ_spec = pickle.load(fp)
        for studies_organ, organ in zip(all_studies, organs):
            ctype = "disease fibroblast"
            get_genes(
                organ_spec[organ],
                ctype,
                studies_organ,
                organ,
                gene_dict_up,
                gene_dict_down,
            )

    # Save gene lists for this scenario
    gene_lists_up = gene_dict_to_long_df(
        gene_dict_up,
        "up",
        real_names,
    )

    gene_lists_down = gene_dict_to_long_df(
        gene_dict_down,
        "down",
        real_names,
    )

    gene_lists = pd.concat(
        [gene_lists_up, gene_lists_down],
        ignore_index=True,
    )

    if results == "all":
        gene_lists = gene_lists[gene_lists['ctype'] == 'mesenchymal']
        gene_lists.to_csv(
            snakemake.output["csv_gene_lists_all"],
            index=False,
        )

    elif results == "myofib":
        gene_lists.to_csv(
            snakemake.output["csv_gene_lists_myofib"],
            index=False,
        )



    percentage_unique = {}
    for ctype in gene_dict_up.keys():
        percentage_unique[ctype] = {}
        for name, reg_genes in zip(["up", "down"], [gene_dict_up, gene_dict_down]):
            upset_df = upsetplot.from_contents(reg_genes[ctype])
            percentages = get_percentage_unique(upset_df, organs, real_names)
            percentage_unique[ctype][name] = percentages
    percentage_unique_results[results] = percentage_unique

# Convert results to DataFrames for plotting
unique_genes_summary_all = pd.concat(
    [
        pd.DataFrame.from_dict(percentage_unique_results["all"][i])
        for i in percentage_unique_results["all"]
    ],
    keys=percentage_unique_results["all"].keys(),
).reset_index()
myofib_summ = pd.concat(
    [
        pd.DataFrame.from_dict(percentage_unique_results["myofib"][i])
        for i in percentage_unique_results["myofib"]
    ],
    keys=percentage_unique_results["myofib"].keys(),
).reset_index()



# Plot and save as PDF
fig, ax = plt.subplots(2, tight_layout=True, sharex=True, figsize=(7, 7))
sns.barplot(
    data=unique_genes_summary_all,
    x="level_0",
    hue="level_1",
    y="up",
    ax=ax[0],
    palette=org_color_real,
)
sns.barplot(
    data=unique_genes_summary_all,
    x="level_0",
    hue="level_1",
    y="down",
    ax=ax[1],
    palette=org_color_real,
)
ax[0].set_xlabel("cell type")
ax[1].set_xlabel("cell type")
ax[0].get_legend().remove()
ax[1].legend(loc="upper left", bbox_to_anchor=(1, 1))

ax[0].set_title("percentage unique")
ax[0].set_ylim(0, 100)
ax[1].set_ylim(0, 100)
ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation=45, ha="right")
plt.savefig(output_pdf_all_mesench, bbox_inches="tight")


# Combine mesenchymal and disease fibroblast data for direct comparison in one plot
mesenchymal_up = unique_genes_summary_all[
    unique_genes_summary_all["level_0"] == "mesenchymal"
].copy()
mesenchymal_up["scenario"] = "all mesench. cells"  # Change legend label
myofib_up = myofib_summ.copy()
myofib_up["scenario"] = "disease fibroblasts"
# Concatenate for 'up' genes
combined_up = pd.concat(
    [
        mesenchymal_up[["level_1", "up", "scenario"]],
        myofib_up[["level_1", "up", "scenario"]],
    ]
)
# Concatenate for 'down' genes
combined_down = pd.concat(
    [
        mesenchymal_up[["level_1", "down", "scenario"]],
        myofib_up[["level_1", "down", "scenario"]],
    ]
)
# Plot both in one figure, side by side
fig, ax = plt.subplots(1, 2, tight_layout=True, sharey=True, figsize=(10, 4))
bar_up = sns.barplot(
    data=combined_up,
    x="level_1",
    y="up",
    hue="scenario",
    ax=ax[0],
    palette={'#801523', '#77B3A5'},
)
bar_down = sns.barplot(
    data=combined_down,
    x="level_1",
    y="down",
    hue="scenario",
    ax=ax[1],
    palette={'#801523', '#77B3A5'},
)
ax[0].set_title("up-regulated genes")
ax[1].set_title("down-regulated genes")
ax[0].set_ylabel("percentage unique (%)")
for a in ax:
    a.set_xlabel("organ")
    a.set_ylim(0, 100)
    a.set_xticklabels(a.get_xticklabels(), rotation=45, ha="right")
# Only show the legend on the first subplot
handles, labels = ax[0].get_legend_handles_labels()
ax[1].legend(handles, labels, loc="upper left", bbox_to_anchor=(1, 1))
ax[0].get_legend().remove()
plt.savefig(output_pdf_comparison, bbox_inches="tight")

# Save plotted data as CSV files

unique_genes_summary_all.to_csv(snakemake.output["csv_all"], index=False)
myofib_summ.to_csv(snakemake.output["csv_myofib"], index=False)
