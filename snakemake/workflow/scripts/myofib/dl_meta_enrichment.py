# Script for enrichment analysis of mixed effects model meta-analysis results
#
# This script loads organ-level model results and gene set resources, performs pathway and TF enrichment
# using decoupler, and generates dot plots for commonly up- and downregulated pathways/TFs across organs.
# Results are saved as pickle files and PDF plots for downstream analysis. Designed for Snakemake workflows.

import pickle

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# snakemake inputs
organ_spec_dl = snakemake.input["organ_spec_dl"]


# snakemake parameters
org_colors = snakemake.params["organ_colors"]
organs = org_colors.keys()
views = snakemake.params["views"]
real_names = snakemake.params["organ_names"]
enrichment = snakemake.params["enrichment"]

collectri_path = enrichment["collectri"]
progeny_path = enrichment["progeny"]

# snakemake outputs
enrichment_on_dl_path = snakemake.output["enrichment_on_dl_path"]
enrichment_results_pval = snakemake.output["enrichment_results_pval"]
enrichment_results_coef = snakemake.output["enrichment_results_coef"]


org_color_real = {
    real_names_key: org_colors[real_names_dict_key]
    for real_names_dict_key, real_names_key in real_names.items()
}

# Set plot parameters
plt.rcParams.update({"font.size": 20})

# load model results
with open(organ_spec_dl, "rb") as fp:
    organ_spec = pickle.load(fp)


# load gene sets for enrichment

msigdb = dc.get_resource("MSigDB")
# Filter by hallmark & go terms
hallmark = msigdb[msigdb["collection"] == "hallmark"]
goterms = msigdb[msigdb["collection"] == "go_biological_process"]

# Remove duplicated entries
hallmark = hallmark[~hallmark.duplicated(["geneset", "genesymbol"])]
goterms = goterms[~goterms.duplicated(["geneset", "genesymbol"])]


# Get pathway activity inference  scores
pathways = pd.read_csv(progeny_path)
pathways["collection"] = "progeny"
# Get TF activity inference score
tfs = pd.read_csv(collectri_path)
tfs["collection"] = "collecttri"
goterms_mf = msigdb[msigdb["collection"] == "go_molecular_function"]

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


rel_info_dict = {}

for organ in organs:
    df = organ_spec[organ].reset_index().set_index("gene")
    genes_to_keep = df.groupby(df.index)["w_re"].apply(lambda x: (x.dropna() > 0).all())
    print(f"filtered out {len(df.index.unique()) - len(genes_to_keep)} genes.")
    working_model = df.loc[genes_to_keep]

    working_model = working_model[working_model["summary_row"] == "random effect"]
    working_model["organ"] = organ
    working_model.loc[:, 'z'] = working_model.loc[:, "eff"] / working_model.loc[:, "sd_eff"]
    working_model = working_model.loc[:, ["organ", "z"]]
    rel_info_dict[organ] = working_model
all_df = pd.concat([rel_info_dict[i] for i in organs])
all_df = all_df.pivot(columns="organ", values="z").fillna(0)


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

    pval_dict[collection] = ulm_pva
    ulm_dict[collection] = ulm_estimate

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
    size_range=(0, 4),
    size_title="-log10(adj. pval)",
    color_title="log2(FC)",
    title="Disease fibroblasts",
    x_order=None,
    sig_threshold=1.301,
):
    """
    Plot differential expression results as a dot plot.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing differential expression results.

    x : str
        Column used for x-axis.

    y : str
        Column used for y-axis.

    color : str
        Column used for dot color.

    size : str
        Column used for dot size.

    color_range : tuple
        Min/max values for color normalization.

    size_range : tuple
        Min/max values for size normalization.

    size_title : str
        Title for size legend.

    color_title : str
        Title for colorbar.

    title : str
        Plot title.

    x_order : list or None
        Optional explicit ordering of x-axis categories.

    sig_threshold : float
        Threshold above which points receive a black outline.
        Default: 1.301 (~ -log10(0.05))
    """

    if df[y].nunique() == 0:
        return None

    df = df.copy()

    if x_order is not None:
        df[x] = pd.Categorical(
            df[x],
            categories=x_order,
            ordered=True,
        )

    fig, ax = plt.subplots(
        1,
        figsize=(
            round(len(df[x].dropna().unique()) / 4) + 0.5,
            (0.3 * len(df[y].unique()) + 1),
        ),
        tight_layout=True,
    )

    # Marker size range (points²)
    marker_sizes = (20, 200)

    sns.scatterplot(
        ax=ax,
        data=df,
        x=x,
        y=y,
        size=size,
        hue=color,
        sizes=marker_sizes,
        palette="RdBu_r",
        size_norm=size_range,
        hue_norm=color_range,
        legend=False,
    )

    # ------------------------------------------------------------------
    # Add black rings around significant points
    # ------------------------------------------------------------------
    min_s, max_s = marker_sizes

    size_norm_obj = plt.Normalize(
        size_range[0],
        size_range[1],
        clip=True,
    )

    ring_sizes = (
        min_s
        + size_norm_obj(df[size])
        * (max_s - min_s)
    )

    ring_colors = np.where(
        df[size] > sig_threshold,
        "black",
        "lightgrey",
    )

    ax.scatter(
        x=df[x],
        y=df[y],
        s=ring_sizes,
        facecolors="none",
        edgecolors=ring_colors,
        linewidths=1.2,
        zorder=10,
    )

    # ------------------------------------------------------------------
    # Axis limits
    # ------------------------------------------------------------------
    yrange = [-0.5, df[y].nunique() + 2]
    ax.set_ylim(yrange)

    xrange = [-1, df[x].nunique()]
    ax.set_xlim(xrange)

    # ------------------------------------------------------------------
    # Colorbar
    # ------------------------------------------------------------------
    norm = plt.Normalize(
        color_range[0],
        color_range[1],
    )

    sm = plt.cm.ScalarMappable(
        cmap="RdBu_r",
        norm=norm,
    )
    sm.set_array([])

    fig.subplots_adjust(right=0.95)

    sub_ax = plt.axes([1, 0.35, 0.1, 0.25])

    ax.figure.colorbar(
        sm,
        label=color_title,
        cax=sub_ax,
    )

    # ------------------------------------------------------------------
    # Custom size legend
    # ------------------------------------------------------------------
    size_norm_obj = plt.Normalize(
        size_range[0],
        size_range[1],
    )

    size_legend_values = np.linspace(
        size_range[0],
        size_range[1],
        num=5,
    )

    min_s, max_s = marker_sizes

    size_legend_handles = [
        plt.scatter(
            [],
            [],
            s=min_s + size_norm_obj(v) * (max_s - min_s),
            color="black",
        )
        for v in size_legend_values
    ]

    ax.legend(
        size_legend_handles,
        [f"{v:.1f}" for v in size_legend_values],
        loc="center left",
        bbox_to_anchor=(1.8, 0.5),
        title=size_title,
        frameon=False,
    )

    # ------------------------------------------------------------------
    # Formatting
    # ------------------------------------------------------------------
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)

    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        ha="right",
        rotation_mode="anchor",
    )

    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("")

    # ------------------------------------------------------------------
    # Organ annotation row
    # ------------------------------------------------------------------
    labels = ax.get_xticklabels()

    organ_colors_bar = [
        org_color_real.get(label.get_text())
        for label in labels
    ]

    y_max = df[y].nunique()

    ax.scatter(
        x=np.arange(len(labels)),
        y=[y_max + 1] * len(labels),
        c=organ_colors_bar,
        s=100,
        marker="s",
        zorder=5,
    )

    return fig, ax

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


        df_col = ulm_dict[collection]
        pval_col = pval_dict[collection]

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
                fig = plot_dotplot(subset, x = 'index', y = 'collection', color_range=(-4, 4),
                color_title='enrichment', title = f'{collection}, in disease fibroblasts',
                x_order = ['heart', 'lung','kidney', 'liver'])
                output_pdf_organs.savefig(bbox_inches="tight")

    

