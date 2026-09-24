# Script for meta-analysis and visualization of mixed effects model of differential expression results
#
# This script loads organ- and cross-organ-level model results of disease fibroblast vs. rest,
# processes and filters significant genes,
# generates forest plots and UpSet plots for up- and downregulated genes, and exports gene lists for downstream use.

import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import upsetplot
from matplotlib.backends.backend_pdf import PdfPages

# snakemake inputs
organ_spec_dl = snakemake.input["organ_spec_dl"]
cross_organ_dl = snakemake.input["cross_organ_dl"]
deg_files = snakemake.input["deg_file"]


# snakemake parameters
org_colors = snakemake.params["organ_colors"]
organs = org_colors.keys()
real_names = snakemake.params["organ_names"]

# snakemake outputs
cross_organ_per_dl = snakemake.output["cross_organ_dl"]
upsetplot_dl_organspec_genes = snakemake.output["upsetplot_dl_organspec_genes"]
gene_dict_up_output = snakemake.output["gene_dict_up_output"]
gene_dict_down_output = snakemake.output["gene_dict_down_output"]
count_organspec_genes = snakemake.output["count_organspec_genes"]
counts_up = snakemake.output["counts_up"]
counts_down = snakemake.output["counts_down"]
topx = snakemake.output["topx"]



input_df = pd.DataFrame(columns=["organ", "study", "deg_file"])
for count, deg_file in enumerate(deg_files):
    organ = deg_file.split("/")[-2]
    study = deg_file.split("/")[-1][0:-4]
    input_df.loc[count, "organ"] = organ
    input_df.loc[count, "study"] = study
    input_df.loc[count, "deg_file"] = deg_file

all_studies = input_df["study"].unique()

org_color_real = {
    real_names_key: org_colors[real_names_dict_key]
    for real_names_dict_key, real_names_key in real_names.items()
}


# plot paraneters
plt.rcParams.update({"font.size": 20})


def get_organ_for_study(study):
    """
    function to receive organ name when study name is given
    """
    organ = input_df[input_df["study"] == study]["organ"].iloc[0]
    return organ


def plot_forest_plot(
    df,
    low_lim=-0.1,
    high_lim=4,
    ctype="disease-associated fibroblast",
    heat_columns=["heart", "lung", "kidney", "liver"],
):
    if df.shape[0] == 0:
        return None
    # Create a figure with two subplots side by side.
    fig, (ax_forest, ax_heat) = plt.subplots(
        1,
        2,
        figsize=(8, len(df) * 0.7),
        gridspec_kw={"width_ratios": [2, 2]},
        sharey=True,
    )

    # Create a numeric index for the genes (to use on the y-axis)
    y_positions = np.arange(len(df))

    ### Forest Plot on ax_forest ###
    # Plot the effect sizes with error bars.
    # Calculate error bar lengths from the confidence intervals.
    lower_errors = df["eff"] - df["ci_low"]
    upper_errors = df["ci_upp"] - df["eff"]

    ax_forest.errorbar(
        df["eff"],
        y_positions + 0.5,
        xerr=[lower_errors, upper_errors],
        fmt="o",
        color="black",
        capsize=5,
    )

    # Add a vertical line at 0 for reference.
    ax_forest.axvline(x=0, color="grey", linestyle="--")

    # Set the y-ticks to show gene names.
    ax_forest.set_yticks(y_positions)
    ax_forest.set_yticklabels(df.index)
    ax_forest.set_xlabel("effect size")
    #ax_forest.set_title(f"{ctype} cells")

    # Invert y-axis so that the first gene is at the top.
    ax_forest.invert_yaxis()
    ax_forest.set_xlim(low_lim, high_lim)

    heat_data = df[heat_columns].fillna(-100)

    # Create a heatmap. We use a diverging colormap (or a binary one) to show True/False.
    sns.heatmap(
        heat_data,
        ax=ax_heat,
        cbar=True,
        cmap="seismic",
        linewidths=0.5,
        linecolor="gray",
        annot=True,
        fmt = '.1f',
        annot_kws={"size": 18},
        vmax=2.5,
        vmin=-2.5,
        mask=(heat_data == -100),
    )

    # Adjust the heatmap's appearance.
    ax_heat.set_title("organ effect sizes")
    ax_heat.set_xlabel("organ")
    # Remove the y-axis label from the heatmap (since it's redundant).
    ax_heat.set_ylabel("")

    # Tight layout for better spacing.
    plt.tight_layout()
    return fig


with open(organ_spec_dl, "rb") as fp:
    organ_spec = pickle.load(fp)


with open(cross_organ_dl, "rb") as fp:
    co = pickle.load(fp)


with PdfPages(cross_organ_per_dl) as output_pdf:
    fibs = co.reset_index().set_index("gene")
    all_genes = fibs.groupby(["gene", "summary_row"]).count().reset_index()
    heart_genes = all_genes[all_genes["summary_row"] == "heart"]["gene"]
    lung_genes = all_genes[all_genes["summary_row"] == "lung"]["gene"]
    kidney_genes = all_genes[all_genes["summary_row"] == "kidney"]["gene"]
    liver_genes = all_genes[all_genes["summary_row"] == "liver"]["gene"]
    fibs["organ_nr"] = fibs.groupby("gene").size() - 4
    fibs.loc[heart_genes, "heart"] = fibs[fibs["summary_row"] == "heart"]["eff"][
        heart_genes
    ]
    fibs.loc[lung_genes, "lung"] = fibs[fibs["summary_row"] == "lung"]["eff"][
        lung_genes
    ]
    fibs.loc[kidney_genes, "kidney"] = fibs[fibs["summary_row"] == "kidney"]["eff"][
        kidney_genes
    ]
    fibs.loc[liver_genes, "liver"] = fibs[fibs["summary_row"] == "liver"]["eff"][
        liver_genes
    ]
    # filter genes where model did not assign weights properly
    genes_to_keep = fibs.groupby(fibs.index)["w_re"].apply(
        lambda x: (x.dropna() > 0).all()
    )
    working_model = fibs.loc[genes_to_keep]

    # filter genes with sig. effect size
    sig_results = working_model[
        (working_model["summary_row"] == "random effect")
        & (working_model["eff"] > 0.5)
        & (working_model["sd_eff"].notna())
        & (working_model["ci_low"] > 0)
    ]

    df = sig_results.sort_values(by="gene")
    df = sig_results.sort_values(by="gene")
    df = df.fillna(0)
    df_sub = df[
        (df["organ_nr"] > 3)
        & (df["heart"] >= 0)
        & (df["lung"] >= 0)
        & (df["kidney"] >= 0)
        & (df["liver"] >= 0)
    ]
    plot = plot_forest_plot(
        df_sub, low_lim=-0.1, high_lim=7, ctype="disease-associated fibroblast"
    )
    output_pdf.savefig(bbox_inches="tight")

    # filter genes with sig. effect size
    sig_results_neg = working_model[
        (working_model["summary_row"] == "random effect")
        & (working_model["eff"] < -0.5)
        & (working_model["sd_eff"].notna())
        & (working_model["ci_upp"] < 0)
    ]

    df = sig_results_neg.sort_values(by="gene")
    df = sig_results_neg.sort_values(by="gene")
    df = df.fillna(0)
    df_sub = df[
        (df["organ_nr"] > 3)
        & (df["heart"] <= 0)
        & (df["lung"] <= 0)
        & (df["kidney"] <= 0)
        & (df["liver"] <= 0)
    ]
    plot = plot_forest_plot(
        df_sub, low_lim=-7, high_lim=0.1, ctype="disease-associated fibroblast"
    )
    output_pdf.savefig(bbox_inches="tight")



with PdfPages(topx) as output_pdf:
    top = 10

    fibs = co.reset_index().set_index('gene')
    all_genes = fibs.groupby(['gene', 'summary_row']).count().reset_index()
    heart_genes = all_genes[all_genes['summary_row'] == 'heart']['gene']
    lung_genes = all_genes[all_genes['summary_row'] == 'lung']['gene']
    kidney_genes = all_genes[all_genes['summary_row'] == 'kidney']['gene']
    liver_genes = all_genes[all_genes['summary_row'] == 'liver']['gene']
    fibs['organ_nr'] = fibs.groupby('gene').size() - 4
    fibs.loc[heart_genes, 'heart'] = fibs[fibs['summary_row'] == 'heart']['eff'][heart_genes]
    fibs.loc[lung_genes, 'lung'] = fibs[fibs['summary_row'] == 'lung']['eff'][lung_genes] 
    fibs.loc[kidney_genes, 'kidney'] = fibs[fibs['summary_row'] == 'kidney']['eff'][kidney_genes]
    fibs.loc[liver_genes, 'liver'] = fibs[fibs['summary_row'] == 'liver']['eff'][liver_genes]
    # filter genes where model did not assign weights properly
    test = fibs.loc[fibs['summary_row'] == 'random effect', ]
    sub = test[test['organ_nr']>3]

    top_hits_up = list(sub.sort_values(by = 'ci_low').iloc[-top:].index)
    top_hits_down = list(sub.sort_values(by = 'ci_upp').iloc[:top].index)

    working_model = fibs.loc[fibs['summary_row'] == 'random effect', ].loc[top_hits_up, :]
    df = working_model.sort_values(by = 'gene')
    plot = plot_forest_plot(df, low_lim= -0.5, high_lim=3.5, ctype = 'disease-assoc. fibroblast')
    output_pdf.savefig(bbox_inches="tight")


    working_model = fibs.loc[fibs['summary_row'] == 'random effect', ].loc[top_hits_down, :]
    df = working_model.sort_values(by = 'gene')
    plot = plot_forest_plot(df, low_lim= -3.5, high_lim=0.5, ctype = 'disease-assoc. fibroblast')
    output_pdf.savefig(bbox_inches="tight")





gene_dict_up = {}
gene_dict_down = {}
for organ in organs:
    studies_organ = input_df[input_df["organ"] == organ]["study"].unique()

    fibs = organ_spec[organ].reset_index().set_index("gene")
    all_genes = fibs.groupby(["gene", "summary_row"]).count().reset_index()
    for study in studies_organ:
        genes = all_genes[all_genes["summary_row"] == study]["gene"]
        fibs.loc[genes, f"{study}"] = fibs[fibs["summary_row"] == f"{study}"]["eff"][
            genes
        ]
    fibs["study_nr"] = fibs.groupby("gene").size() - 4

    # filter genes where model did not assign weights properly
    genes_to_keep = fibs.groupby(fibs.index)["w_re"].apply(
        lambda x: (x.dropna() > 0).all()
    )
    working_model = fibs.loc[genes_to_keep]

    # filter genes with sig. effect size
    sig_results = working_model[
        (working_model["summary_row"] == "random effect")
        & (working_model["eff"] > 0.5)
        & (working_model["sd_eff"].notna())
        & (working_model["ci_low"] > 0)
    ]
    df_sub = sig_results.sort_values(by="gene")
    gene_dict_up[organ] = list(df_sub.index)

    # now only dowwnregulated genes
    # filter genes with sig. effect size
    sig_results = working_model[
        (working_model["summary_row"] == "random effect")
        & (working_model["eff"] < -0.5)
        & (working_model["sd_eff"].notna())
        & (working_model["ci_upp"] < 0)
    ]
    df_sub = sig_results.sort_values(by="gene")
    gene_dict_down[organ] = list(df_sub.index)


with open(gene_dict_up_output, "wb") as fp:
    pickle.dump(gene_dict_up, fp)
with open(gene_dict_down_output, "wb") as fp:
    pickle.dump(gene_dict_down, fp)


with PdfPages(upsetplot_dl_organspec_genes) as output_pdf:
    for name, reg_genes in zip(["up", "down"], [gene_dict_up, gene_dict_down]):
        upset_df = upsetplot.from_contents(reg_genes)
        upset_df.index.rename(
            names=[real_names.get(name, name) for name in upset_df.index.names],
            inplace=True,
        )
        upsetplot.plot(upset_df, show_counts=True, sort_categories_by="input", element_size = 30)
        plt.suptitle(f"{name}regulated genes", size=18)

        output_pdf.savefig()


        # Find genes present in at least two organs

        atleast2 = np.where(upset_df.reset_index()[[real_names[i] for i in organs]].sum(axis=1) > 1)
        # Print summary of these genes for inspection

        print(f"{name}:")
        print(upset_df.iloc[atleast2])



with PdfPages(count_organspec_genes) as output_pdf:
    gene_number_dict = {}
    for name, dict in zip(['upregulated','downregulated'], [gene_dict_up, gene_dict_down]):
        # Prepare subplots
        fig, ax = plt.subplots(1, figsize=(4, 4), sharey=True, tight_layout=True)
        data = []
    
        for organ, genes in dict.items():

            organ_names = []
            gene_counts = []
            colors = []

            organ_names.append(real_names[organ])
            gene_counts.append(len(genes))
            colors.append(org_color_real.get(real_names[organ], 'gray'))  

            # Save data for the DataFrame
            data.append({
                'Organ': real_names[organ],
                'Gene Count': len(genes)
            })

        ax.bar(organ_names, gene_counts, color=colors)
        ax.set_xticklabels(organ_names, rotation=45, ha='right')
        ax.set_ylabel(f"Number of genes \n ({name})")

        fig.tight_layout()
        plt.show()

    
        count_summary = pd.DataFrame(data)
        if name == 'upregulated':
            count_summary.to_csv(counts_up)
        elif name == 'downregulated':
            count_summary.to_csv(counts_down)

        output_pdf.savefig()