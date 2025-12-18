# This script processes patient-level scores 
# It computes per-study means and standard deviations for each gene set, then calculates z-scores
# and assigns immune/ECM status to each patient based on these standardized values.


import pandas as pd


est_path = snakemake.input['out_path_est']
pval_path = snakemake.input['out_path_pval']
save_results_path = snakemake.output[0]
patient_scoring_sets = snakemake.params[1]

full_genesets = pd.read_csv(patient_scoring_sets, index_col = 0)
all_sets = list(full_genesets['geneset'].unique())

# process patients scores and save processed results in dictionaries
result_dict = {}
pval_dict = {}
for est_p, pval_p in zip(est_path, pval_path):

    print(pval_p)
    print(est_p)

    organ = pval_p.split("/")[-2]
    study = '_'.join(pval_p.split("/")[-1].split("_")[0:-1])

    print(organ)
    print(study)

    ulm_est = (pd.read_csv(est_p)
               .rename(columns = {'Unnamed: 0':'patient'})
              )
    ulm_est['organ'] = organ
    ulm_est = ulm_est[all_sets  + ['cond_test','study','patient', 'organ','grouping']]
    
    pval = (pd.read_csv(pval_p)
            .rename(columns = {'Unnamed: 0':'patient'})
           )
    pval['organ'] = organ
    pval = pval[all_sets + ['study', 'patient', 'organ']]

    
    result_dict[study] = ulm_est
    pval_dict[study] = pval


# Concatenate all results and p-values from all studies into single DataFrames
full_results_df = pd.concat(result_dict, ignore_index=True)
full_pval_df = pd.concat(pval_dict, ignore_index=True)

# Define columns for study/condition grouping
study_cond_columns = list(all_sets + ['cond_test'] + ['study'])

# Calculate mean of each gene set per study and condition
grouped = full_results_df.loc[:, study_cond_columns].groupby(['study', 'cond_test']).mean().reset_index()
study_ctrl_means = grouped[grouped['cond_test'] == 'control']
# Merge mean values for control samples back to the full results
merged_df = pd.merge(full_results_df, study_ctrl_means, on='study', suffixes=('', '_mean'))

# Also calculate the standard deviation per study and condition
# Use absolute value to avoid negative stdevs (for robustness)
grouped_std = full_results_df.loc[:, study_cond_columns].groupby(['study', 'cond_test']).std().apply(abs).reset_index()
study_ctrl_std = grouped_std[grouped_std['cond_test'] == 'control']
# Merge stdev values for control samples back to the merged DataFrame
merged_df = pd.merge(merged_df, study_ctrl_std, on='study', suffixes=('', '_std'))

# Compute difference and standardized difference (z-score) for each gene set
# Store column names for later reference
diff_columns = []
stdv_columns = []
for col in all_sets:
    # Difference from control mean
    merged_df[col + '_diff'] = merged_df[col] - merged_df[col + '_mean']
    diff_columns.append(str(col + '_diff'))
    # Standardized difference (z-score)
    merged_df[col + '_stdev'] = (merged_df[col] - merged_df[col + '_mean']) / merged_df[col + '_std']
    stdv_columns.append(str(col + '_stdev'))

# Assign immune and ECM status based on z-score thresholds
# 'high' if z-score >= 1, 'low' otherwise
# Default to 'none' if not assigned

to_save_patients = merged_df
to_save_patients['immune'] = 'none'
to_save_patients.loc[to_save_patients['GOBP_IMMUNE_RESPONSE_stdev'] >= 1, 'immune'] = 'high'
to_save_patients.loc[to_save_patients['GOBP_IMMUNE_RESPONSE_stdev'] < 1, 'immune'] = 'low'

to_save_patients['ecm'] = 'none'
to_save_patients.loc[to_save_patients['NABA_CORE_MATRISOME_stdev'] >= 1, 'ecm'] = 'high'
to_save_patients.loc[to_save_patients['NABA_CORE_MATRISOME_stdev'] < 1, 'ecm'] = 'low'

# Save the final patient-level results to CSV

to_save_patients.to_csv(save_results_path)