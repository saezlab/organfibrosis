import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import yaml
import os
from skbio.stats.composition import clr, multiplicative_replacement
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from statsmodels.formula.api import mixedlm

# plot paraneters
plt.rcParams.update({'font.size':18})

# Snakemake bindings supply all run-time configuration, inputs, and desired outputs.
config = snakemake.config
organs = config['meta_organs']
fibroblast_path = str(snakemake.input["fibroblast"])
scvi_path = str(snakemake.input["scvi"])
marker_paths = list(map(str, snakemake.input["markers"]))
stacked_bar_path = str(snakemake.output["stacked_bar"])
umap_output = dict(zip(organs, map(str, snakemake.output["umap"])))
proportion_output = dict(zip(organs, map(str, snakemake.output["proportion"])))
top5_output = dict(zip(organs, map(str, snakemake.output["top5"])))
top5_std_output = dict(zip(organs, map(str, snakemake.output["top5_std"])))
pan_output = dict(zip(organs, map(str, snakemake.output["pan"])))
pan_std_output = dict(zip(organs, map(str, snakemake.output["pan_std"])))
marker_files = dict(zip(organs, marker_paths))

# Prepare the destination directories for every PDF figure that will be produced.
all_output_files = [stacked_bar_path]
all_output_files.extend(umap_output.values())
all_output_files.extend(proportion_output.values())
all_output_files.extend(top5_output.values())
all_output_files.extend(top5_std_output.values())
all_output_files.extend(pan_output.values())
all_output_files.extend(pan_std_output.values())

for path in all_output_files:
    directory = os.path.dirname(path)
    if directory:
        os.makedirs(directory, exist_ok=True)

all_studies = [config['datasets'][organ]['studies_to_analyze'] for organ in organs]
all_studies_for_cols = list(config['general_plotting']['study_colors'].keys())
organ_colors = config['general_plotting']['organ_colors']
real_names = config['general_plotting']['organ_names']
study_color_dict = config['general_plotting']['study_colors']

org_color_real = {
    real_names_key: organ_colors[real_names_dict_key] 
    for real_names_dict_key, real_names_key in real_names.items()
}
condition_colors = {'control':'#102542','fibrosis':'#828282'}
n_studies = np.array([len(studies) for studies in all_studies]).sum()

views = config['general_plotting']['views']
resolutions = config['preprocessing']['fibroblasts']['scvi_resolutions']
diseasefibs = config['preprocessing']['fibroblasts']['diseasefibs']
diseasefib_colors = config['general_plotting']['diseasefibs']


def read_and_add_harmony_results(adata_to_add, path, suffix, id_study = False, pcs = 50):
    # Attach precomputed harmony/UMAP embeddings and Leiden labels to the supplied AnnData object.
    """
    function to read in PCA, UMAP, and leiden clustering results and 
    add them to adata object with the suffix so that several results can be added
    input:
    adata_to_add - adata object where results will be added
    path - path at wich PCA, UMAP and leiden results are saved
    suffix - suffix is added to names, i.e. UMPA results is saved at X_umap_suffix
    id_study - if id column is just the cell id or a merged cell id + study
    pcs - no. of pcs

    output:
    adata object with adde results
    """
    harm_results = pd.read_csv(path, 
                     index_col = 0, 
                     dtype = {'id':'str','leiden_0.2':str, 
                              'leiden_0.3':str,'leiden_0.4':str,
                              'leiden_0.5':str,'leiden_0.7':str})

    if id_study == False:
        assert(np.all(adata_to_add.obs_names == harm_results['id']))
    elif id_study == True:
        harm_results = harm_results.set_index('Unnamed: 1')
        harm_results = harm_results.loc[(np.array(adata_to_add.obs.index) + '_' + np.array(adata_to_add.obs['study']))]
        print('success')
    else:
        assert(np.all((np.array(adata_to_add.obs.index) + '_' + np.array(adata_to_add.obs['study'])) == harm_results['id']))
        harm_results = harm_results.set_index('id',drop = False)
        print('success')

    adata_to_add.obsm[f'X_harmony_{suffix}'] = np.array(harm_results.iloc[:,:pcs])
    adata_to_add.obsm[f'X_umap_{suffix}'] = np.array(harm_results.iloc[:,pcs + 1:pcs + 3])

    adata_to_add.obs[[f'leiden_0.2_{suffix}',f'leiden_0.3_{suffix}',f'leiden_0.4_{suffix}',
                      f'leiden_0.5_{suffix}',f'leiden_0.7_{suffix}']] = harm_results[['leiden_0.2',
                                                                                        'leiden_0.3',
                                                                                        'leiden_0.4',
                                                                                        'leiden_0.5',
                                                                                        'leiden_0.7']].to_numpy()

    return adata_to_add


def get_cell_count(adata, clustering):
    # Compute cell type counts per sample
    cell_counts = (
        adata.obs
        .groupby(['sample_study', clustering])
        .size()
        .unstack(fill_value=0)
    )

    # Convert counts to proportions
    prop_data = cell_counts.div(cell_counts.sum(axis=1), axis=0)

    # CLR transformation
    prop_data_clr = clr(multiplicative_replacement(prop_data.values))
    prop_data_clr = pd.DataFrame(prop_data_clr, index=prop_data.index, columns=prop_data.columns)

    return(prop_data_clr)


def analyse_proportions(prop_data_clr, meta_data, clustering):
    # Melt for analysis
    long_df = prop_data_clr.reset_index().melt(id_vars='sample_study', var_name=clustering, value_name='value')
    long_df = long_df.merge(meta_data, left_on='sample_study', right_index = True)

    # T-tests by study and cell type
    t_results = []
    for (study, cell_type), group in long_df.groupby(['study', clustering]):
        hf_vals = group[group['cond_test'] == 'fibrosis']['value']
        nf_vals = group[group['cond_test'] == 'control']['value']
        t_stat, p_val = ttest_ind(hf_vals, nf_vals, equal_var=False)
        t_results.append({'study': study, clustering: cell_type, 'statistic': t_stat, 'p.value': p_val})

    study_diff_stats = pd.DataFrame(t_results)

    # T and P matrices
    t_matrix = study_diff_stats.pivot(index=clustering, columns='study', values='statistic')
    p_matrix = study_diff_stats.pivot(index=clustering, columns='study', values='p.value')

    # Adjusted p-values
    adj_p_matrix = p_matrix.apply(lambda col: multipletests(col, method='fdr_bh')[1], axis=0)
    adj_p_matrix = pd.DataFrame(adj_p_matrix, index=p_matrix.index, columns=p_matrix.columns)

    # Linear mixed models
    lmer_results = []
    for cell_type, group in long_df.groupby(clustering):
        group['cond_test'] = pd.Categorical(group['cond_test'], categories=['control', 'fibrosis'])
        model = mixedlm("value ~ cond_test", group, groups=group["study"])
        result = model.fit()
        estimate = result.params['cond_test[T.fibrosis]']
        p_val = result.pvalues['cond_test[T.fibrosis]']
        var_components = result.cov_re.iloc[0, 0] / (result.cov_re.iloc[0, 0] + result.scale)
        lmer_results.append({clustering: cell_type, 'Estimate': estimate, 'p_val': p_val, 'perc_studyvar': var_components})

    study_diff_stats_lmer = pd.DataFrame(lmer_results)
    study_diff_stats_lmer['adj_pval'] = multipletests(study_diff_stats_lmer['p_val'], method='fdr_bh')[1]

    # Star annotation matrix
    star_matrix = adj_p_matrix.applymap(lambda p: "*" if p < 0.05 else "")

    return(t_matrix, star_matrix, study_diff_stats_lmer)



# Functions for color logic
def color_fun(estimates, p_vals, up_color="blue", down_color="red", nonsig_color="grey"):
    return [
        up_color if est > 0 and p < 0.05
        else down_color if est < 0 and p < 0.05
        else nonsig_color
        for est, p in zip(estimates, p_vals)
    ]

def heatmap_color_fun(matrix, up_color="blue", down_color="red"):
    min_val, max_val = np.min(matrix), np.max(matrix)
    def color_mapper(x):
        if x < 0:
            return down_color
        elif x > 0:
            return up_color
        else:
            return "white"
    return np.vectorize(color_mapper)(matrix)

# read anndata and add harmonization results
adata = sc.read(fibroblast_path)
adata = read_and_add_harmony_results(adata, scvi_path, "scvi", id_study=True, pcs=30)

# annotate fibroblasts
# Tag disease-associated fibroblast clusters for each organ according to the config.
adata.obs["substate"] = "rest"
for organ in organs:
    resolution = resolutions[organ]
    clusters = diseasefibs[organ]
    adata.obs.loc[
        (adata.obs["organ"] == organ)
        & (adata.obs[f"leiden_{resolution}_scvi"].isin(clusters)),
        "substate",
    ] = "diseasefib"


adata.obs["sample_study"] = (
    adata.obs["sample"].astype("str") + "_" + adata.obs["study"].astype("str")
)



# UMAP overview per organ
for organ in organs:
    print(organ)
    resolution = resolutions[organ]
    plot = adata[adata.obs['organ'] == organ].copy()
    n_clusters = plot.obs[f'leiden_{resolution}_scvi'].nunique()
    print(n_clusters)
    fig = sc.pl.embedding(
            plot,
            basis = 'X_umap_scvi',
            color = [f'leiden_{resolution}_scvi','substate'], 
            #palette = condition_colors,
            show = False,
            title = [f'leiden clusters', f'cell states'],
            return_fig = True,
            palette = None,
            wspace = 0.2,
            )
    for ax in fig.axes:
        ax.set_xlabel("UMAP 1") 
        ax.set_ylabel("UMAP 2")  
    fig.suptitle(f"{real_names[organ]}", fontsize=20, y=1.05) 
    fig.savefig(umap_output[organ])
    plt.close(fig)



# Distribution of substate proportions per study
fig, axs = plt.subplots(2,2, tight_layout = True, figsize = (7,7))
ax = axs.ravel()
for count, organ in enumerate(organs):
    print(organ)
    resolution = resolutions[organ]
    df = adata[adata.obs['organ'] == organ].obs

    # Step 1: Create a crosstab of counts
    count_table = pd.crosstab(df['study'], df['substate'])

    # Step 2: Normalize to get percentages (row-wise)
    percentage_table = count_table.div(count_table.sum(axis=1), axis=0)

    # Step 3: Plot a stacked bar chart
    percentage_table.plot(kind='bar', stacked=True, ax = ax[count])

    ax[count].set_ylabel("proportion")
    ax[count].set_title(real_names[organ])
    ax[count].get_legend().remove()

ax[3].legend(loc='upper left', bbox_to_anchor=(1, 1))
fig.savefig(stacked_bar_path)
plt.close(fig)



# Differential abundance summary per organ
for organ in organs:
    resolution = resolutions[organ]
    adata_o = adata[adata.obs['organ'] == organ]
    adata_o.obs['sample_study'] = adata_o.obs['sample'].astype(str) + '_' + adata_o.obs['study'].astype(str)

    meta_data = adata_o.obs[['sample_study', 'study', 'cond_test']].drop_duplicates().set_index('sample_study')
    clustering = f'leiden_{resolution}_scvi'
    prop_data_clr = get_cell_count(adata_o, clustering)
    
    
    t_matrix, star_matrix, study_diff_stats_lmer = analyse_proportions(prop_data_clr, meta_data,clustering)
    fig, ax = plt.subplots(1, 2, width_ratios=[3, 6], figsize=(7, 8))

    sns.barplot(data = study_diff_stats_lmer, y = clustering, x = 'Estimate', ax = ax[0])
    ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation=45, ha='right')
    ax[0].set_ylabel("leiden cluster")
    sns.heatmap(
        t_matrix,
        annot=star_matrix,
        fmt="s",
        cmap=sns.diverging_palette(240, 10, n=9, as_cmap=True),
        linewidths=0.5,
        cbar_kws={"label": "t-value (dis. fibroblast vs. rest)"}, ax = ax[1], vmin = -4, vmax = 4
    )
    plt.title(f"{real_names[organ]}")
    ax[0].set_xlabel(f"estimate")
    plt.ylabel("leiden cluster")
    plt.xlabel("study")
    plt.tight_layout()
    fig.savefig(proportion_output[organ])
    plt.close(fig)



# Marker gene visualization per organ
pan_tissue_mf_markers = [
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "FN1",
    "TNC",
    "POSTN",
    "FAP",
    "PDGFRA",
    "SFRP2",
    "TIMP1"
]




# plot top 5 marker genes per organ & common MF markers
marker_dict = {}
for organ in organs:

    adata_o = adata[adata.obs['organ'] == organ]
    marker_output = marker_files[organ]
    markers = pd.read_csv(
        marker_output,
        index_col=0
    )
    top5 = markers.loc[0:5, f'{diseasefibs[organ][0]}_leiden_{resolutions[organ]}'] 
    
    dotplot_top5_std = sc.pl.dotplot(adata_o, top5, groupby="substate", standard_scale="var", title = real_names[organ], show=False)
    fig_top5_std = getattr(dotplot_top5_std, "figure", None) or getattr(dotplot_top5_std, "fig", None) or plt.gcf()
    fig_top5_std.savefig(top5_std_output[organ], bbox_inches='tight')
    plt.close(fig_top5_std)

    dotplot_top5 = sc.pl.dotplot(adata_o, top5, groupby="substate", title = real_names[organ], show=False)
    fig_top5 = getattr(dotplot_top5, "figure", None) or getattr(dotplot_top5, "fig", None) or plt.gcf()
    fig_top5.savefig(top5_output[organ], bbox_inches='tight')
    plt.close(fig_top5)

    dotplot_pan_std = sc.pl.dotplot(adata_o, pan_tissue_mf_markers, groupby="substate", standard_scale="var", title = real_names[organ], show=False)
    fig_pan_std = getattr(dotplot_pan_std, "figure", None) or getattr(dotplot_pan_std, "fig", None) or plt.gcf()
    fig_pan_std.savefig(pan_std_output[organ], bbox_inches='tight')
    plt.close(fig_pan_std)

    dotplot_pan = sc.pl.dotplot(adata_o, pan_tissue_mf_markers, groupby="substate", title = real_names[organ], show=False)
    fig_pan = getattr(dotplot_pan, "figure", None) or getattr(dotplot_pan, "fig", None) or plt.gcf()
    fig_pan.savefig(pan_output[organ], bbox_inches='tight')
    plt.close(fig_pan)
