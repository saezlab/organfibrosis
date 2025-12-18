# Import required libraries
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import upsetplot
from matplotlib.backends.backend_pdf import PdfPages

# Snakemake inputs
organ_spec_dl = snakemake.input["organ_spec_dl"]
cross_organ_dl = snakemake.input["cross_organ_dl"]
deg_files = snakemake.input["deg_file"]

# Snakemake parameters
org_colors = snakemake.params["organ_colors"]
organs = org_colors.keys()
views = snakemake.params["views"]
real_names = snakemake.params["organ_names"]

# Snakemake outputs
cross_organ_per_ctype_dl = snakemake.output["cross_organ_per_ctype_dl"]
upsetplot_dl_organspec_genes = snakemake.output["upsetplot_dl_organspec_genes"]
gene_dict_up_output = snakemake.output["gene_dict_up_output"]
gene_dict_down_output = snakemake.output["gene_dict_down_output"]
count_organspec_genes = snakemake.output["count_organspec_genes"]
counts_up = snakemake.output["counts_up"]
counts_down = snakemake.output["counts_down"]
topx = snakemake.output["topx"]

# Placeholder for cell type name
ctype = "not known yet"

# Build input dataframe describing all studies and organs
input_df = pd.DataFrame(columns=["organ", "study", "deg_file"])
for count, deg_file in enumerate(deg_files):
    organ = deg_file.split("/")[-2]
    study = "_".join(deg_file.split("/")[-1].split("_")[0:-3])
    input_df.loc[count, "organ"] = organ
    input_df.loc[count, "study"] = study
    input_df.loc[count, "deg_file"] = deg_file

all_studies = input_df["study"].unique()

# Map real organ names to colors for plotting
org_color_real = {
    real_names_key: org_colors[real_names_dict_key]
    for real_names_dict_key, real_names_key in real_names.items()
}

# Set plot parameters
plt.rcParams.update({"font.size": 20})


def get_organ_for_study(study):
    """
    Return organ name for a given study name.
    """
    organ = input_df[input_df["study"] == study]["organ"].iloc[0]
    return organ


def plot_forest_plot(
    df,
    low_lim=-0.1,
    high_lim=4,
    ctype="",
    heat_columns=["heart", "lung", "kidney", "liver"],
):
    """
    Create a forest plot with optional heatmap for organ effect sizes.

    Parameters:
    - df: DataFrame containing the data to plot.
    - low_lim, high_lim: Limits for the x-axis.
    - ctype: Cell type, used for the plot title.
    - heat_columns: List of columns to display in the heatmap.

    Returns:
    - fig: The created figure.
    """
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
    ax_forest.set_title(f"{ctype}")

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
        annot_kws={"size": 19},
        vmax=1,
        vmin=-1,
        mask=(heat_data == -100),
    )

    # Adjust the heatmap's appearance.
    ax_heat.set_title("organ effect sizes")
    ax_heat.set_xlabel("organ")
    # Remove the y-axis label from the heatmap (since it's redundant).
    ax_heat.set_ylabel("")
    plt.tight_layout()
    return fig


with open(organ_spec_dl, "rb") as fp:
    organ_spec = pickle.load(fp)


with open(cross_organ_dl, "rb") as fp:
    co = pickle.load(fp)


with PdfPages(cross_organ_per_ctype_dl) as output_pdf:
    for ctype in views:
        fibs = co[ctype].reset_index().set_index("gene")
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
        # epithelial cells are only present in three organs
        if ctype != "epithelial":
            df_sub = df[
                (df["organ_nr"] == 4)
                & (df["heart"] > 0)
                & (df["lung"] > 0)
                & (df["kidney"] > 0)
                & (df["liver"] > 0)
            ]
        else:
            df_sub = df[
                (df["organ_nr"] == 3)
                & (df["heart"] > 0)
                & (df["lung"] > 0)
                & (df["kidney"] > 0)
                & (df["liver"] > 0)
            ]
        plot = plot_forest_plot(df_sub, low_lim=-0.1, high_lim=5, ctype=ctype)
        output_pdf.savefig(bbox_inches="tight")

        # filter genes with sig. effect size
        sig_results_neg = working_model[
            (working_model["summary_row"] == "random effect")
            & (working_model["eff"] < -0.5)
            & (working_model["sd_eff"].notna())
            & (working_model["ci_upp"] < 0)
        ]

        df = sig_results_neg.sort_values(by="gene")
        # plot_forest_plot(df, low_lim= -4, high_lim=0.1, ctype = ctype)

        df = sig_results_neg.sort_values(by="gene")
        if ctype != "epithelial":
            df_sub = df[
                (df["organ_nr"] == 4)
                & (df["heart"] < 0)
                & (df["lung"] < 0)
                & (df["kidney"] < 0)
                & (df["liver"] < 0)
            ]
        else:
            df_sub = df[
                (df["organ_nr"] == 3)
                & (df["heart"] < 0)
                & (df["lung"] < 0)
                & (df["kidney"] < 0)
                & (df["liver"] < 0)
            ]
        plot = plot_forest_plot(df_sub, low_lim=-5, high_lim=0.1, ctype=ctype)
        output_pdf.savefig(bbox_inches="tight")


with PdfPages(topx) as output_pdf:
    for ctype in views:
        top = 10
        fibs = co[ctype].reset_index().set_index('gene')
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
        if ctype != 'epithelial':
            sub = test[test['organ_nr']>3]
        else:
            sub = test[test['organ_nr']>2]
        print(ctype)

        top_hits_up = list(sub.sort_values(by = 'ci_low').iloc[-top:].index)
        top_hits_down = list(sub.sort_values(by = 'ci_upp').iloc[:top].index)


        working_model = fibs.loc[fibs['summary_row'] == 'random effect', ].loc[top_hits_up, :]
        df = working_model.sort_values(by = 'gene')
        plot = plot_forest_plot(df, low_lim= -0.5, high_lim=3.0, ctype = ctype)
        output_pdf.savefig(bbox_inches="tight")


        working_model = fibs.loc[fibs['summary_row'] == 'random effect', ].loc[top_hits_down, :]
        df = working_model.sort_values(by = 'gene')
        plot = plot_forest_plot(df, low_lim= -3.0, high_lim=0.5, ctype = ctype)
        output_pdf.savefig(bbox_inches="tight")




gene_dict_up = {}
gene_dict_down = {}
for organ in organs:
    studies_organ = input_df[input_df["organ"] == organ]["study"].unique()
    for ctype in views:
        fibs = organ_spec[ctype][organ].reset_index().set_index("gene")
        all_genes = fibs.groupby(["gene", "summary_row"]).count().reset_index()
        for study in studies_organ:
            genes = all_genes[all_genes["summary_row"] == study]["gene"]
            fibs.loc[genes, f"{study}"] = fibs[fibs["summary_row"] == f"{study}"][
                "eff"
            ][genes]
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
        if ctype in gene_dict_up.keys():
            gene_dict_up[ctype][organ] = list(df_sub.index)
        else:
            gene_dict_up[ctype] = {}
            gene_dict_up[ctype][organ] = list(df_sub.index)

        # Now process only downregulated genes
        # Filter genes with significant negative effect size
        sig_results = working_model[
            (working_model["summary_row"] == "random effect")
            & (working_model["eff"] < -0.5)
            & (working_model["sd_eff"].notna())
            & (working_model["ci_upp"] < 0)
        ]
        df_sub = sig_results.sort_values(by="gene")
        # Save downregulated gene list for this cell type and organ
        if ctype in gene_dict_down.keys():
            gene_dict_down[ctype][organ] = list(df_sub.index)
        else:
            gene_dict_down[ctype] = {}
            gene_dict_down[ctype][organ] = list(df_sub.index)

# Save gene dictionaries as pickle files for downstream use
with open(gene_dict_up_output, "wb") as fp:
    pickle.dump(gene_dict_up, fp)
with open(gene_dict_down_output, "wb") as fp:
    pickle.dump(gene_dict_down, fp)

# Create UpSet plots for organ-specific up- and downregulated genes
with PdfPages(upsetplot_dl_organspec_genes) as output_pdf:
    for ctype in views:
        for name, reg_genes in zip(["up", "down"], [gene_dict_up, gene_dict_down]):
            upset_df = upsetplot.from_contents(reg_genes[ctype])
            upset_df.index.rename(
                names=[real_names.get(name, name) for name in upset_df.index.names],
                inplace=True,
            )
            upsetplot.plot(upset_df, show_counts=True, sort_categories_by="input", element_size = 30)
            plt.suptitle(f"{name}regulated genes in {ctype} cells", size=16)
            # Save the current plot to the PDF
            output_pdf.savefig()
            # Find genes present in at least two organs
            atleast2 = np.where(upset_df.reset_index()[[real_names[i] for i in organs]].sum(axis=1) > 1)
            # Print summary of these genes for inspection
            print(f"{name} in {ctype}:")
            print(upset_df.iloc[atleast2])


with PdfPages(count_organspec_genes) as output_pdf:
    for name, dict in zip(['upregulated','downregulated'], [gene_dict_up, gene_dict_down]):
        # Prepare subplots
        n_celltypes = len(dict)
        fig, axes = plt.subplots(1, n_celltypes, figsize=(3 * n_celltypes, 4), sharey=True, tight_layout=True)
        data = []
        for ax, (cell_type, organs) in zip(axes, dict.items()):
            organ_names = []
            gene_counts = []
            colors = []

            for organ, genes in organs.items():
                organ_names.append(real_names[organ])
                gene_counts.append(len(genes))
                colors.append(org_color_real.get(real_names[organ], 'gray'))  

                # Save data for the DataFrame
                data.append({
                    'Cell Type': cell_type,
                    'Organ': real_names[organ],
                    'Gene Count': len(genes)
                })

            ax.bar(organ_names, gene_counts, color=colors)
            ax.set_title(cell_type)
            
            ax.set_xticklabels(organ_names, rotation=45, ha='right')
        axes[0].set_ylabel(f"Number of genes \n ({name})")

        fig.tight_layout()
        plt.show()

    
        count_summary = pd.DataFrame(data)
        count_summary_filt = count_summary[
            ~((count_summary['Cell Type'] == 'epithelial') & (count_summary['Organ'] == 'heart'))
            ]
        if name == 'upregulated':
            count_summary_filt.to_csv(counts_up)
        elif name == 'downregulated':
            count_summary_filt.to_csv(counts_down)
            
        output_pdf.savefig()