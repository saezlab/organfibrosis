# Author: Leonie Küchenhoff
# This script performs differential gene expression meta-analysis across multiple organs
# using data from various studies. It combines results for each organ and then does a 
# mixed effects model.
# The results are saved in pickle files for further analysis.

import pandas as pd
from statsmodels.stats.meta_analysis import combine_effects
import pickle

# Snakemake inputs
deg_files = snakemake.input['deg_file']

# Snakemake parameters
org_colors = snakemake.params['organ_colors']
organs = org_colors.keys()
views = snakemake.params['views']
real_names = snakemake.params['organ_names']

# Snakemake outputs
organ_spec_path = snakemake.output['organ_spec']
cross_organ_path = snakemake.output['cross_organ']

# Prepare input DataFrame mapping organs, studies, and DEG files
input_df = pd.DataFrame(columns=['organ', 'study', 'deg_file'])
for count, deg_file in enumerate(deg_files):
    organ = deg_file.split("/")[-2]
    study = '_'.join(deg_file.split("/")[-1].split("_")[0:-3])
    input_df.loc[count, 'organ'] = organ
    input_df.loc[count, 'study'] = study
    input_df.loc[count, 'deg_file'] = deg_file

all_studies = input_df['study'].unique()

def prep_deg_files_allgenes(ctype):
    """
    Returns differential gene expression results for a cell type of interest for all genes.
    Args:
        ctype (str): Cell type of interest
    Returns:
        stat_df (pd.DataFrame): logFC results per study and gene
        se_df (pd.DataFrame): SE value results per study and gene
    """
    stat_col = f'log2FoldChange_{ctype}'
    se_col = f'lfcSE_{ctype}'
    stat_df = pd.DataFrame()
    se_df = pd.DataFrame()
    for index, row in input_df.iterrows():
        study = row['study']
        organ = row['organ']
        deg = pd.read_csv(row['deg_file'], index_col=0)
        cols = [col for col in deg.columns if ctype in col]
        if len(cols) > 0:
            deg_subset = deg.loc[:, cols].rename(columns={stat_col: study})
            stat = deg_subset[study]
            stat_df = stat_df.merge(stat, left_index=True, right_index=True, how='outer')
            se = deg_subset.loc[:, deg_subset.columns != study].rename(columns={se_col: study})
            se = se[study]
            se_df = se_df.merge(se, left_index=True, right_index=True, how='outer')
        else:
            stat_df[study] = -1
            se_df[study] = -1
    stat_df = stat_df.fillna(-1)
    se_df = se_df.fillna(-1)
    return stat_df, se_df

def combine_gene_summary(gene, ctype, studies):
    """
    For a given gene, extract effect sizes and SEs, filter out missing data, and perform meta-analysis if possible.
    Returns the full summary_frame DataFrame from combine_effects, or None if insufficient data.
    """
    logfc_row = deg_dict[ctype]['logFC'].loc[gene, studies]
    se_row = deg_dict[ctype]['se'].loc[gene, studies]
    valid_mask = (logfc_row != -1) & (se_row != -1)
    if valid_mask.sum() < 2:
        return None
    effect_sizes = logfc_row[valid_mask]
    variances = (se_row[valid_mask]) ** 2
    study_names = effect_sizes.index.tolist()
    result = combine_effects(effect_sizes, variances, method_re="iterated", row_names=study_names)
    return result.summary_frame()

def get_organ_info(df, organ_name):
    """
    Get info from one organ from organ-specific meta-analysis results.
    Returns effect and SD as Series.
    """
    org = df.reset_index().set_index('gene')
    org_effect = org[org['summary_row'] == 'random effect']['eff'].rename(organ_name).astype('float64')
    org_sd = org[org['summary_row'] == 'random effect']['sd_eff'].rename(organ_name).astype('float64')
    return org_effect, org_sd

def cross_organ_dl(gene, effect_df, se_df):
    """
    For the given gene, extract effect sizes and SEs across organs, filter out missing values, and perform meta-analysis.
    Returns the full summary_frame DataFrame from combine_effects, or None if insufficient data.
    """
    eff_row = effect_df.loc[gene]
    se_row = se_df.loc[gene]
    valid_mask = (~eff_row.isna()) & (~se_row.isna())
    if valid_mask.sum() < 2:
        return None
    valid_effects = eff_row[valid_mask]
    valid_variances = (se_row[valid_mask]) ** 2
    study_names = valid_effects.index.tolist()
    result = combine_effects(valid_effects, valid_variances, method_re='iterated', row_names=study_names)
    return result.summary_frame()

def flatten_concatenation(matrix):
    flat_list = []
    for row in matrix:
        flat_list += row
    return flat_list

# Load data and run models
deg_dict = {}
deg_count = pd.DataFrame(index=flatten_concatenation(all_studies))
for ctype in views:
    deg_dict[ctype] = {}
    stat_df, se_df = prep_deg_files_allgenes(ctype)
    deg_dict[ctype]['logFC'] = stat_df
    deg_dict[ctype]['se'] = se_df

# Run meta-analysis model per organ
organ_results = {}
for ctype in views:
    print(ctype)
    organ_results[ctype] = {}
    for organ in organs:
        print(organ)
        studies = input_df[input_df['organ'] == organ]['study'].unique()
        print(studies)
        if (organ == 'reheatHeart') & (ctype == 'epithelial'):
            organ_results[ctype][organ] = pd.DataFrame(
                columns=['gene', 'summary_row', 'eff', 'sd_eff', 'ci_low', 'ci_upp', 'w_fe', 'w_re']
            ).set_index(['gene', 'summary_row'])
        else:
            results_dict = {
                gene: combine_gene_summary(gene, ctype, studies)
                for gene in deg_dict[ctype]['se'].index
            }
            combined_results = pd.concat(results_dict, names=['gene', 'summary_row'])
            organ_results[ctype][organ] = combined_results

# Save organ-specific results
with open(organ_spec_path, 'wb') as fp:
    pickle.dump(organ_results, fp)

# Combine organ-specific results into one score (cross-organ meta-analysis)
co_results = {}
for ctype in views:
    lung_effect, lung_se = get_organ_info(organ_results[ctype]['HCAlung'], 'lung')
    heart_effect, heart_se = get_organ_info(organ_results[ctype]['reheatHeart'], 'heart')
    kidney_effect, kidney_se = get_organ_info(organ_results[ctype]['kidney'], 'kidney')
    liver_effect, liver_se = get_organ_info(organ_results[ctype]['liver'], 'liver')
    effect_df = pd.concat([lung_effect, heart_effect, kidney_effect, liver_effect], axis=1)
    se_df = pd.concat([lung_se, heart_se, kidney_se, liver_se], axis=1)
    results_dict_cross_org = {
        gene: cross_organ_dl(gene, effect_df, se_df)
        for gene in effect_df.index
    }
    combined_results_cross_org = pd.concat(results_dict_cross_org, names=['gene', 'summary_row'])
    co_results[ctype] = combined_results_cross_org

# Save cross-organ results
with open(cross_organ_path, 'wb') as fp:
    pickle.dump(co_results, fp)
