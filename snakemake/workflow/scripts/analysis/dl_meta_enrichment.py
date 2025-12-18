# Import required libraries
import pickle

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# Snakemake inputs
organ_spec_dl = snakemake.input["organ_spec_dl"]

# Snakemake parameters
org_colors = snakemake.params["organ_colors"]
organs = org_colors.keys()
views = snakemake.params["views"]
real_names = snakemake.params["organ_names"]
enrichment = snakemake.params["enrichment"]
collectri_path = enrichment["collectri"]
progeny_path = enrichment["progeny"]

# Snakemake outputs
enrichment_on_dl_path = snakemake.output["enrichment_on_dl_path"]
enrichment_results_pval = snakemake.output["enrichment_results_pval"]
enrichment_results_coef = snakemake.output["enrichment_results_coef"]

# Placeholder for cell type name
ctype = "not known yet"

# Map real organ names to colors for plotting
org_color_real = {
    real_names_key: org_colors[real_names_dict_key]
    for real_names_dict_key, real_names_key in real_names.items()
}

# Set plot parameters
plt.rcParams.update({"font.size": 20})

# Load model results from pickle file
with open(organ_spec_dl, "rb") as fp:
    organ_spec = pickle.load(fp)

# Load gene sets for enrichment analysis
msigdb = dc.get_resource("MSigDB")
# Filter for hallmark and GO biological process terms
hallmark = msigdb[msigdb["collection"] == "hallmark"]
goterms = msigdb[msigdb["collection"] == "go_biological_process"]
# Remove duplicated entries for clean gene sets
hallmark = hallmark[~hallmark.duplicated(["geneset", "genesymbol"])]
goterms = goterms[~goterms.duplicated(["geneset", "genesymbol"])]

# Load pathway and TF activity inference scores
pathways = pd.read_csv(progeny_path)
pathways["collection"] = "progeny"
tfs = pd.read_csv(collectri_path)
tfs["collection"] = "collecttri"
goterms_mf = msigdb[msigdb["collection"] == "go_molecular_function"]

# Build network for enrichment analysis
net = msigdb[msigdb["collection"].isin(["go_biological_process", "hallmark"])].rename(
    columns={"geneset": "source", "genesymbol": "target"}
)
net["weight"] = 1
net = pd.concat([net, pathways, tfs])
net = net[~net.duplicated(["source", "target"])]
net = net[net["target"].notna()]
collections = net["collection"].unique()

ulm_dict = {}
pval_dict = {}

# Run enrichment analysis for each cell type and organ
for ctype in views:
    rel_info_dict = {}
    ulm_dict[ctype] = {}
    pval_dict[ctype] = {}
    for organ in organs:
        df = organ_spec[ctype][organ].reset_index().set_index("gene")
        # Filter genes with valid model weights
        genes_to_keep = df.groupby(df.index)["w_re"].apply(
            lambda x: (x.dropna() > 0).all()
        )
        print(f"filtered out {len(df.index.unique()) - len(genes_to_keep)} genes.")
        working_model = df.loc[genes_to_keep]
        # Keep only random effect rows
        working_model = working_model[working_model["summary_row"] == "random effect"]
        working_model.loc[:, 'z'] = working_model.loc[:, "eff"] / working_model.loc[:, "sd_eff"]
        working_model["organ"] = organ
        working_model = working_model.loc[:, ["organ", "z"]]
        rel_info_dict[organ] = working_model
    # Combine all organs into a single DataFrame
    all_df = pd.concat([rel_info_dict[i] for i in organs])
    all_df = all_df.pivot(columns="organ", values="z").fillna(0)
    # Run enrichment for each collection
    for collection in collections:
        net_subset = net[net["collection"] == collection]
        ulm_estimate, ulm_pva = dc.run_ulm(
            all_df.T,
            net=net_subset,
            source="source",
            target="target",
            verbose=True,
            weight="weight",
        )
        pval_dict[ctype][collection] = ulm_pva
        ulm_dict[ctype][collection] = ulm_estimate

# Save enrichment results as pickle files
with open(enrichment_results_pval, "wb") as outfile:
    pickle.dump(pval_dict, outfile)
with open(enrichment_results_coef, "wb") as outfile:
    pickle.dump(ulm_dict, outfile)


def plot_dotplot(
    df,
    x="study",
    y="index",
    color="FC",
    size="pval_log",
    color_range=(-1.5, 1.5),
    size_range=(0, 3),
    size_title="-log10(adj, pval)",
    color_title="log2(FC)",
    title=ctype,
):
    """
    This function plot the differential expression results of a list of genes in a certain cell type
    Input:
    genes - list of genes
    df - pandas dataframe with differntial gene expression results
        example: 	index	study	FC	pval_log
                0	ABCA3	Misharin_Budinger_2018	0.00000	0.000000
                1	AC002066.1	Misharin_Budinger_2018	0.000000	0.000000
    x = x axis column
    y = y axis columnqq
    color - column used for color of dots
    size- column used for size of dots
    color_range - min and max values of color as tuple
    size_range - min and max values of size as tuple
    """
    if df[y].nunique() > 0:
        fig, ax = plt.subplots(
            1,
            figsize=(
                round(len(df[x].unique()) / 4) + 0.5,
                (0.1 * len(df[y].unique()) + 1.5),
            ),
            tight_layout=True,
        )

        # Scatterplot with limited size range
        sns.scatterplot(
            ax=ax,
            data=df,
            x=x,
            y=y,
            size=size,
            hue=color,
            sizes=((size_range[0] * 60) + 20, (size_range[1] * 60) + 20),
            palette="RdBu_r",
            size_norm=size_range,
            hue_norm=color_range,
        )

        yrange = [-0.5, df[y].nunique() + 2]
        ax.set_ylim(yrange)
        xrange = [-1, df[x].nunique()]
        ax.set_xlim(xrange)

        norm = plt.Normalize(color_range[0], color_range[1])
        sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm)
        sm.set_array([])

        # Access the legend handles and labels
        handles, labels = ax.get_legend_handles_labels()
        # Filter out the handles and labels for the size legend
        size_legend_handles = [h for i, h in enumerate(handles) if i > 6]
        # Create a custom size legend
        size_legend_labels = np.linspace(
            size_range[0], size_range[1], num=5
        )  # Specify the desired size values
        size_legend_handles = [
            plt.scatter([], [], s=(val * 60) + 20, color="black")
            for val in size_legend_labels
        ]

        ax.legend(
            size_legend_handles,
            [f"{val:.1f}" for val in size_legend_labels],
            loc="center left",
            bbox_to_anchor=(1.8, 0.5),
            title=size_title,
        )

        # Add colorbar
        fig.subplots_adjust(right=0.95)  # create space on the right hand side
        sub_ax = plt.axes([1, 0.35, 0.1, 0.25])  # add a small custom axis
        ax.figure.colorbar(sm, label=color_title, cax=sub_ax)
        ax.tick_params(axis="x", labelsize=10)
        ax.tick_params(axis="y", labelsize=10)
        labels = ax.get_xticklabels()

        ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor")
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("")

        # Adding organ color bar or markers
        organ_colors_bar = [org_color_real.get(study.get_text()) for study in labels]
        y_max = df[y].nunique()
        ax.scatter(
            x=np.arange(len(df[x].unique())),
            y=[y_max + 1] * len(df[x].unique()),
            c=organ_colors_bar,
            s=100,
            marker="s",
            label="Organ",
        )


def check_if_gene_everywhere(
    stat_df, pval_df, val_cutoff_pos=2, val_cutoff_neg=-2, pcut=0.05, organ_nr=3
):
    """
    This function checks if a gene is common (checks across organ), meaning it is present in at least 2/3 or 3/4 organs 

    input:
    stat_df - dataframes with test statistics
    pval_df - datafram with p value

    cutoffs for testing
    organ_nr: number of organs to check for commonality

    returns:
    three lists: one with upregulated genes, one with downregulated genes, one with in general regulated genes
    """

    significant_df_down = pd.DataFrame((stat_df < val_cutoff_neg) & (pval_df < pcut))
    summed_weights_down = significant_df_down.sum(axis=1)
    summed_weights_down = list((summed_weights_down[summed_weights_down >= organ_nr]).index)

    significant_df_up = pd.DataFrame((stat_df > val_cutoff_pos) & (pval_df < pcut))
    summed_weights_up = significant_df_up.sum(axis=1)
    summed_weights_up = list((summed_weights_up[summed_weights_up >= organ_nr]).index)

    significant_df_all = pd.DataFrame(
        ((stat_df < val_cutoff_neg) | (stat_df > val_cutoff_pos)) & (pval_df < pcut)
    )
    summed_weights_all = significant_df_all.sum(axis=1)
    summed_weights_all = list((summed_weights_all[summed_weights_all >= organ_nr]).index)

    return summed_weights_down, summed_weights_up, summed_weights_all



def prep_collection_files(plot_df, genes_to_plot):
    subset = plot_df.loc[genes_to_plot, :].reset_index()
    subset["pval_log"] = np.log10(subset["pval"]) * -1
    subset.replace([np.inf, -np.inf], 0, inplace=True)
    return subset



# plot top based on average 
top=20
with PdfPages(enrichment_on_dl_path) as output_pdf_organs:
    for collection in ['collecttri', 'go_biological_process', 'hallmark']:
        for count,ctype in enumerate(views):

            df_col = ulm_dict[ctype][collection]
            pval_col = pval_dict[ctype][collection]

            #pos
            columns_with_positive_values = list(df_col.columns[(df_col > 0).all()])
            features = list(df_col.mean(axis = 0).sort_values()[-1*top:].index)
            overlapping_features_pos = list(set(pval_col.columns) & set(features) & set(columns_with_positive_values))
            # neg
            columns_with_neg_values = list(df_col.columns[(df_col < 0).all()])
            features = list(df_col.mean(axis = 0).sort_values()[:1*top].index)
            overlapping_feature_neg = list(set(pval_col.columns) & set(features) & set(columns_with_neg_values))
    
            for overlapping_features in [overlapping_feature_neg, overlapping_features_pos]:
    
                pval_df = pval_col.loc[:, overlapping_features].rename(index = real_names)
                stat_df = df_col.loc[:, overlapping_features].rename(index = real_names)
    
                stat = stat_df.melt(ignore_index = False, var_name = 'collection', value_name = 'FC')
                pval = pval_df.melt(ignore_index = False, var_name = 'collection', value_name = 'pval')
                
                stat = stat.reset_index()
                pval = pval.reset_index()
    
                plot = stat.merge(pval, on = ['index', 'collection'], how = 'outer').set_index('collection')
    
                subset = prep_collection_files(plot, genes_to_plot = list(plot.index))
                if subset.shape[0] > 0 :
                    if collection == 'go_biological_process':
                        subset['collection'] = subset['collection'].str.replace(
                            'GOBP_', ''
                            ).str.replace(
                            '_', ' '
                            ).str.replace(
                            'DIFFERENTIATION', 'DIFF.'
                            ).str.replace(
                            'MORPHOGENESIS', 'MORPHOGEN.'
                            ).str.replace(
                            'DEVELOPMENT', 'DEV.'
                            ).str.replace(
                            'REGULATION', 'REG.'
                            ).str.lower()
                    fig = plot_dotplot(subset, x = 'index', y = 'collection', color_range=(-2, 2), color_title='enrichment', title = f'{collection}, in {ctype}')
                    output_pdf_organs.savefig(bbox_inches="tight")
    
    



