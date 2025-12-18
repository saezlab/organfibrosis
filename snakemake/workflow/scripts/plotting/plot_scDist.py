import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


# Snakemake inputs
scdist_results = snakemake.input['scdist_results']
output_csv = snakemake.output['scdist_summary_csv']
study_results = snakemake.output['study_results']
all_results_pdf = snakemake.output['all_results']

organs = snakemake.params['organs']
views = snakemake.params['views']
org_colors = snakemake.params["org_colors"]
real_names = snakemake.params["real_names"]
datasets = snakemake.params["datasets"]
study_colors = snakemake.params["study_colors"]

org_color_real = {
    real_names_key: org_colors[real_names_dict_key]
    for real_names_dict_key, real_names_key in real_names.items()
}


all_studies = [datasets[organ]['studies_ownmodel'] for organ in organs]
studies_to_analyze = [datasets[organ]['studies_to_analyze'] for organ in organs]

# plot paraneters
plt.rcParams.update({'font.size':18})
# Read results for all studies/organs
results_dict = {}

for path in scdist_results:
    organ = path.split("/")[-2]
    study = '_'.join(path.split("/")[-1].split("_")[0:-1])
    results = pd.read_csv(path, delimiter=' ')
    results['study'] = study
    results['organ'] = real_names[organ]
    results = results.reset_index()
    results_dict[study] = results

all_results = pd.concat(results_dict, ignore_index=True)


with PdfPages(study_results) as output_pdf:
    # plot results
    for organ, studies_organ in zip(organs, all_studies):
        fig, ax = plt.subplots(1, len(studies_organ), figsize = ((len(studies_organ) * 3) + 2,5), sharey = True)
        for count, study in enumerate(studies_organ):

            subset = all_results[all_results['study'] == study]
            
            ax[count].errorbar(x=subset['index'], y=subset['Dist.'], 
                        yerr=[subset['Dist.'] - subset['95% CI (low)'], subset['95% CI (upper)'] - subset['Dist.']], 
                        fmt='o', capsize=5, label='Dist. with 95% CI')
            
            ax[count].set_xticklabels(ax[count].get_xticklabels(), rotation=45, ha='right')
            ax[count].set_ylabel('distance')
            ax[count].set_title(study, size = 18)
        
        fig.suptitle(f"{subset['organ'].iloc[0]}")
        fig.tight_layout()
        output_pdf.savefig()

def flatten_concatenation(matrix):
    flat_list = []
    for row in matrix:
        flat_list += row
    return flat_list

all_studies_to_analyze = flatten_concatenation(studies_to_analyze)


# Plot Distance per organ/cell type
fig, ax = plt.subplots(
    1, len(views), figsize=(20, 9), sharey=True, sharex=True, tight_layout=True
)
for count, ctype in enumerate(views):
    plot_deg_count_ctype = all_results[all_results['index'] == ctype]
    sns.boxplot(
        ax=ax[count],
        data=plot_deg_count_ctype[plot_deg_count_ctype['study'].isin(all_studies_to_analyze)],
        x="organ",
        y="Dist.",
        boxprops={"alpha": 0.4},
        palette=org_color_real,
    )
    sns.swarmplot(
        ax=ax[count],
        data=plot_deg_count_ctype[plot_deg_count_ctype['study'].isin(all_studies_to_analyze)],
        x="organ",
        y="Dist.",
        hue="study",
        palette=study_colors,
        s=8,
    )
    ax[count].legend().remove()
    ax[count].set_xticklabels(ax[count].get_xticklabels(), rotation=45, ha="right")
    ax[count].set_title(ctype)
ax[0].legend(bbox_to_anchor=(-.3, 1), loc="upper right", title="study", fontsize=14)
plt.suptitle("scDist results per organ/cell type", fontsize=20)
plt.savefig(all_results_pdf)


# Summarize mean and std
mean = all_results.groupby(['organ', 'index'])['Dist.'].mean().reset_index()
std = all_results.groupby(['organ', 'index'])['Dist.'].std().reset_index()
scdist_summ = mean.merge(std, on=['organ', 'index'], suffixes=('_mean', '_std'))

# Save summary as CSV
scdist_summ.to_csv(output_csv, index=False)