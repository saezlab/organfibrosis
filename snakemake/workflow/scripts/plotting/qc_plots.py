# Author: Leonie Küchenhoff
# This script generates quality control (QC) plots 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from skbio.stats.composition import clr
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from matplotlib.backends.backend_pdf import PdfPages
from statsmodels.formula.api import mixedlm

# File paths from snakemake
paths = snakemake.input["qc_df"]
mdata_paths = snakemake.input["metadata_df"]
count_paths = snakemake.input["ctypecount_df"]
count_patients_paths = snakemake.input["ctypecount_patient_df"]


# Output file paths from snakemake
output = snakemake.output["qc_pdf"]
output_metaplot = snakemake.output["metadata"]
ctype_count_pdf_abs = snakemake.output["ctype_count_abs"]
ctype_count_pdf_rel = snakemake.output["ctype_count_rel"]
plot_df_abs_csv = snakemake.output["plot_df_abs_csv"]
plot_df_rel_csv = snakemake.output["plot_df_rel_csv"]
ctype_count_change = snakemake.output["ctype_count_change"]
etiology_count = snakemake.output["etiology_count"]
etiology_count_rel = snakemake.output["etiology_count_rel"]
t_matrix_all_output = snakemake.output["t_matrix_all"]
adj_p_matrix_all_output = snakemake.output["adj_p_matrix_all"]
lm_all_output = snakemake.output["lm_all"]
sex_distribution = snakemake.output['sex_distribution']
age_distribution = snakemake.output['age_distribution']


# Organ and cell type colors, and other parameters from snakemake
organ_colors_light = snakemake.params["organ_colors_light"]
org_colors = snakemake.params["organ_colors"]
ctype_colors = snakemake.params["ctype_colors"]
views = snakemake.params["views"]
real_names = snakemake.params["real_names"]
organs = list(org_colors.keys())
study_sub = snakemake.params["study_sub"]
all_studies_for_cols = list(snakemake.config['general_plotting']['study_colors'].keys())

pairing = snakemake.params["pairing"]
study_organ_df = pd.DataFrame(pairing, columns=["organ", "study"])
study_color_dict = {
    row["study"]: org_colors[row["organ"]] for index, row in study_organ_df.iterrows()
}
category_colors = {'control':"tan",
                    'fibrosis':"darkgreen"}

# Plot parameters
plt.rcParams.update({"font.size": 21})


def analyse_proportions(prop_data_clr, meta_data, clustering):
    # Melt for analysis
    long_df = prop_data_clr.melt(
        id_vars="sample_study", var_name=clustering, value_name="value"
    )
    long_df = long_df.merge(meta_data, on="sample_study")

    # T-tests by study and cell type
    t_results = []
    for (study, cell_type), group in long_df.groupby(["study", clustering]):
        hf_vals = group[group["cond_test"] == "fibrosis"]["value"]
        nf_vals = group[group["cond_test"] == "control"]["value"]
        t_stat, p_val = ttest_ind(hf_vals, nf_vals, equal_var=False)
        t_results.append(
            {
                "study": study,
                clustering: cell_type,
                "statistic": t_stat,
                "p.value": p_val,
            }
        )

    study_diff_stats = pd.DataFrame(t_results)

    # T and P matrices
    t_matrix = study_diff_stats.pivot(
        index=clustering, columns="study", values="statistic"
    )
    p_matrix = study_diff_stats.pivot(
        index=clustering, columns="study", values="p.value"
    )

    # Adjusted p-values
    adj_p_matrix = p_matrix.apply(
        lambda col: multipletests(col, method="fdr_bh")[1], axis=0
    )
    adj_p_matrix = pd.DataFrame(
        adj_p_matrix, index=p_matrix.index, columns=p_matrix.columns
    )

    # Linear mixed models
    lmer_results = []
    for cell_type, group in long_df.groupby(clustering):
        group["cond_test"] = pd.Categorical(
            group["cond_test"], categories=["control", "fibrosis"]
        )
        model = mixedlm("value ~ cond_test", group, groups=group["study"])
        result = model.fit()
        estimate = result.params["cond_test[T.fibrosis]"]
        p_val = result.pvalues["cond_test[T.fibrosis]"]
        var_components = result.cov_re.iloc[0, 0] / (
            result.cov_re.iloc[0, 0] + result.scale
        )
        lmer_results.append(
            {
                clustering: cell_type,
                "Estimate": estimate,
                "p_val": p_val,
                "perc_studyvar": var_components,
            }
        )

    study_diff_stats_lmer = pd.DataFrame(lmer_results)
    study_diff_stats_lmer["adj_pval"] = multipletests(
        study_diff_stats_lmer["p_val"], method="fdr_bh"
    )[1]

    # Star annotation matrix
    star_matrix = adj_p_matrix.applymap(lambda p: "*" if p < 0.05 else "")

    return (t_matrix, star_matrix, study_diff_stats_lmer, adj_p_matrix)


# Prepare empty DataFrames to collect QC metrics for each study
counts_df = pd.DataFrame(columns=["study", "nCount_RNA"])
features_df = pd.DataFrame(columns=["study", "nFeature_RNA"])
pctmt_df = pd.DataFrame(columns=["study", "pct_counts_mt"])

# Loop through each QC file and append relevant columns to the DataFrames
for path in paths:
    qc_file = pd.read_csv(path, dtype={"nFeature_RNA": np.float32})
    study = qc_file.loc[:, "study"].iloc[0]
    counts_df = pd.concat([counts_df, qc_file[["study", "nCount_RNA"]]])
    features_df = pd.concat([features_df, qc_file[["study", "nFeature_RNA"]]])
    pctmt_df = pd.concat([pctmt_df, qc_file[["study", "pct_counts_mt"]]])

# Create violin plots for each QC metric across studies
fig, ax = plt.subplots(3, figsize=(11, 18), tight_layout=True)
sns.violinplot(
    data=counts_df, x="study", y="nCount_RNA", ax=ax[0], palette=study_color_dict
)
sns.violinplot(
    data=features_df, x="study", y="nFeature_RNA", ax=ax[1], palette=study_color_dict
)
sns.violinplot(
    data=pctmt_df, x="study", y="pct_counts_mt", ax=ax[2], palette=study_color_dict
)

# Set titles for each subplot
ax[0].set_title("total counts")
ax[0].set_ylim(-20, 15000)
ax[1].set_title("total features")
ax[1].set_ylim(-20, 10000)
ax[2].set_title("% mito. counts")

# Rotate x-axis labels for better readability
for i in ax.flat:
    i.set_xticklabels(i.get_xticklabels(), rotation=45, ha="right")

# Save the QC plots to file
plt.savefig(output, bbox_inches="tight")

# Metadata plots

# Count the number of reference and fibrosis samples per study
counts = {"Fibrosis": [], "Reference": [], "File": [], "Organ": []}
# Dictionary to store grouping information for each study
group_dict = {}
sex_dict = {}
full_meta_dict = {}
for file_path in mdata_paths:
    # Read metadata table for each study
    df = pd.read_csv(file_path)
    studyname = file_path.split("/")[-1].split(".")[0]
    # sex
    if 'sex' in df.columns:
        sex = df.groupby(['organ','study', 'sex']).size().reset_index().rename(columns = {0:'count'})
        sex['sex'] = sex['sex'].str.lower().fillna('unknown')
        sex['relative'] = sex['count'] / sex['count'].sum()
    else:
        sex = df.groupby(['organ','study']).size().reset_index().rename(columns = {0:'count'})
        sex['sex'] = 'unknown'
        sex['relative'] = 1
    sex_dict[studyname] = sex
    full_meta_dict[studyname] = df
    # Grouping overview for each study
    group = (
        df.groupby(["organ", "study", "grouping"])
        .size()
        .reset_index()
        .rename(columns={0: "count"})
    )
    group["grouping"] = group["grouping"].fillna("unknown")
    group['relative'] = (group['count'] / group['count'].sum()) * 100
    group_dict[studyname] = group
    # Count fibrosis and reference samples
    fibrosis_count = df[df["cond_test"] == "fibrosis"].shape[0]
    reference_count = df[df["cond_test"] == "control"].shape[0]
    organ = study_organ_df[study_organ_df["study"] == studyname]["organ"].values[0]
    # Store counts for plotting
    counts["Fibrosis"].append(fibrosis_count)
    counts["Reference"].append(reference_count)
    counts["File"].append(studyname)
    counts["Organ"].append(organ)

# Create DataFrame from sample counts
counts_df = pd.DataFrame(counts)
sex_df = pd.concat(sex_dict.values(), ignore_index=True)
meta_df = pd.concat(full_meta_dict.values(), ignore_index=True)

# Plot stacked bar plot for fibrosis and reference sample counts
fig, ax = plt.subplots(1, figsize=(11, 4))
bars = ax.bar(counts_df["File"], counts_df["Fibrosis"], label="Fibrosis", linewidth =0)
bars_reference = ax.bar(
    counts_df.index,
    counts_df["Reference"],
    label="Reference",
    bottom=counts_df["Fibrosis"],
    color="#828282",
)
# Color fibrosis bars by organ
for bar, study in zip(bars, counts_df["File"]):
    bar.set_color(study_color_dict.get(study, "grey"))
# Create legends for organs and reference
organ_legend = [
    Patch(facecolor=color, label=real_names[organ]) for organ, color in org_colors.items()
]
reference_patch = Patch(facecolor="#828282", label="reference")
plt.legend(
    handles=organ_legend + [reference_patch], loc="upper left", bbox_to_anchor=(1, 1)
)
#ax.set_title("Number of fibrosis and reference samples per study")
ax.set_xlabel("dataset")
ax.set_ylabel("sample count")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
plt.savefig(output_metaplot, bbox_inches="tight")

# Cell type number overview plots
# Read cell type count files and combine into a single DataFrame
dataframes = []
for file_path in count_paths:
    df = pd.read_csv(file_path)
    study_name = df.columns[1]
    df.columns = ["CellType", study_name]
    dataframes.append(df)
combined_df = pd.concat([df.set_index("CellType") for df in dataframes], axis=1).fillna(
    0
)
combined_df = combined_df.loc[views + ["other"]]
# Plot both absolute and relative cell type counts
plot_df_abs = combined_df.T
rel_values = combined_df / combined_df.sum(axis=0)
plot_df_rel = rel_values.T
ctype_colors["other"] = "lightgrey"
for save, plot_df in zip(
    [ctype_count_pdf_abs, ctype_count_pdf_rel], [plot_df_abs, plot_df_rel]
):
    # Adjust order of columns
    plot_df = plot_df.loc[:, views + ["other"]]
    # Plot stacked bar plot for cell type counts
    fig, ax = plt.subplots(figsize=(13,6), tight_layout=True)
    plot_df.plot(kind="bar", stacked=True, color=ctype_colors, ax=ax, width=0.9)
    # Add organ color bar above bars
    organ_colors_bar = [
        org_colors.get(
            study_organ_df[study_organ_df["study"] == study]["organ"].values[0],
            "lightgrey",
        )
        for study in plot_df.index
    ]
    y_max = ax.get_ylim()[1]
    ax.scatter(
        x=np.arange(len(plot_df.index)),
        y=[y_max] * len(plot_df.index),
        c=organ_colors_bar,
        s=100,
        marker="s",
        label="Organ",
    )
    # Create legends for cell types and organs
    celltype_legend = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=cell_type,
            markersize=10,
            markerfacecolor=color,
        )
        for cell_type, color in ctype_colors.items()
    ]
    organ_legend = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=organ,
            markersize=10,
            markerfacecolor=color,
        )
        for organ, color in org_colors.items()
    ]
    first_legend = ax.legend(
        handles=celltype_legend,
        title="cell types",
        loc="upper left",
        bbox_to_anchor=(1, 1),
    )
    ax.add_artist(first_legend)
    ax.legend(
        handles=organ_legend, title="organs", loc="upper left", bbox_to_anchor=(1, 0.55)
    )
    ax.set_title("Cell types across studies")
    ax.set_xlabel("study")
    ax.set_ylabel("relative cell count")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(save, bbox_inches="tight")

# Save cell type count DataFrames to CSV
plot_df_abs.to_csv(plot_df_abs_csv)
plot_df_rel.to_csv(plot_df_rel_csv)

# Summarize cell type count differences for each patient
# Read patient-level cell type count files and combine
### new ###
dataframes = []
for path in count_patients_paths:
    df = pd.read_csv(path).iloc[:, 1:]
    study_name = path.split("/")[-1].split(".")[0].replace("_ctypecount_patient", "")
    df["organ"] = path.split("/")[-2]
    df["study"] = study_name
    dataframes.append(df)
combined_df = pd.concat([df for df in dataframes], axis=0).fillna(0)


# drop kidney studies that were excluded from the analysis
studies_to_exclude =list(
    set(combined_df[combined_df["organ"] == "kidney"]['study'].unique()) - set(study_sub)
    )
combined_df = combined_df[~combined_df["study"].isin(studies_to_exclude)]


cell_sums = combined_df[views + ["other"]].sum(axis=1)
combined_df[views + ["other"]] = combined_df.loc[:, views + ["other"]].div(
    cell_sums, axis=0
)
colors = [ctype_colors.get(cell_type, "grey") for cell_type in combined_df.index]

t_matrix_dict = {}
adj_p_matrix_dict = {}
lm_dict = {}
# Plot per-organ cell type changes using barplots and heatmaps
with PdfPages(ctype_count_change) as pdf:
    for organ in organs:
        clr_patient_df = combined_df[combined_df["organ"] == organ].copy()
        views_to_use = views.copy()
        # heart does not have epithelial cells
        if organ == "reheatHeart":
            print('excluding epithelial cells for heart')
            clr_patient_df = clr_patient_df.drop("epithelial", axis=1)
            views_to_use.remove("epithelial")
        # Centered log-ratio (clr) transformation for compositional data
        clr_patient_df[clr_patient_df == 0] = 0.0000000000001
        clr_patient_df.loc[:, views_to_use] = clr(clr_patient_df.loc[:, views_to_use])


        clr_patient_df["sample_study"] = (
            clr_patient_df["sample"].astype(str)
            + "_"
            + clr_patient_df["study"].astype(str)
        )
        prop_data_clr = clr_patient_df.loc[:, views_to_use + ["other", "sample_study"]]
        meta_data = clr_patient_df.loc[
            :, ["sample", "sample_study", "cond_test", "study", "other"]
        ]
        t_matrix, star_matrix, study_diff_stats_lmer,adj_p_matrix = analyse_proportions(
            prop_data_clr, meta_data, "ctype"
        )
        fig, ax = plt.subplots(1, 2, width_ratios=[3, 6], figsize=(7, 8))
        sns.barplot(data=study_diff_stats_lmer, y="ctype", x="Estimate", ax=ax[0])
        ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation=45, ha="right")
        sns.heatmap(
            t_matrix,
            annot=star_matrix,
            fmt="s",
            cmap=sns.diverging_palette(240, 10, n=9, as_cmap=True),
            linewidths=0.5,
            cbar_kws={"label": "t-value (fibrosis - reference)"},
            ax=ax[1],
            vmin=-4,
            vmax=4,
        )
        plt.title(f"{real_names[organ]}")
        ax[0].set_ylabel("cell type")
        ax[0].set_xlabel("estimate")
        ax[1].get_yaxis().set_visible(False)
        plt.xlabel("study")
        plt.tight_layout()
        plt.show()
        pdf.savefig(fig)

        t_matrix = t_matrix.T.reset_index()
        adj_p_matrix = adj_p_matrix.T.reset_index()
        study_diff_stats_lmer = study_diff_stats_lmer
        t_matrix['organ'] = organ
        adj_p_matrix['organ'] = organ
        study_diff_stats_lmer['organ'] = organ

        t_matrix_dict[organ] = t_matrix
        adj_p_matrix_dict[organ] = adj_p_matrix
        lm_dict[organ] = study_diff_stats_lmer

t_matrix_all = pd.concat(t_matrix_dict.values(), ignore_index=True)
star_matrix_all = pd.concat(adj_p_matrix_dict.values(), ignore_index=True)
lm_all = pd.concat(lm_dict.values(), ignore_index=True)
# Save the combined DataFrames to CSV
t_matrix_all.to_csv(t_matrix_all_output, index=True)
star_matrix_all.to_csv(adj_p_matrix_all_output, index=True)
lm_all.to_csv(lm_all_output, index=True)    



# Plot disease etiology distributions per study

group_df = pd.concat(group_dict.values(), ignore_index=True)
fig, ax = plt.subplots(1, 4, figsize=(16, 3), sharey=True)
for count, organ in enumerate(group_df["organ"].unique()):
    groups_organ = group_df[group_df["organ"] == organ]
    pivot_df = groups_organ.pivot(index="study", columns="grouping", values="count")
    # Move control column to first position
    first_column = pivot_df.pop("control")
    pivot_df.insert(0, "control", first_column)
    s = pivot_df.plot(kind="bar", stacked=True, ax=ax[count])
    ax[count].set_xlabel("study")
    ax[count].set_ylabel("samples")
    ax[count].legend(title="etiologies", loc="lower left", bbox_to_anchor=(0, 1))
    ax[count].set_xticklabels(ax[count].get_xticklabels(), rotation=45, ha="right")
plt.savefig(etiology_count, bbox_inches="tight")



# Plot relative disease etiology distributions per study

group_df = pd.concat(group_dict.values(), ignore_index=True)
fig, ax = plt.subplots(1, 4, figsize=(16, 3), sharey=True, width_ratios=[5,5,8,4])
for count, organ in enumerate(group_df["organ"].unique()):
    groups_organ = group_df[group_df["organ"] == organ]
    pivot_df = groups_organ.pivot(index="study", columns="grouping", values="relative")
    # Move control column to first position
    first_column = pivot_df.pop("control")
    pivot_df.insert(0, "control", first_column)
    s = pivot_df.plot(kind="bar", stacked=True, ax=ax[count], width = 0.9)
    ax[count].set_xlabel("study")
    ax[count].set_ylabel("relative count [%]")
    ax[count].legend(title="etiologies", loc="lower left", bbox_to_anchor=(0, 1))
    ax[count].set_xticklabels(ax[count].get_xticklabels(), rotation=45, ha="right")
plt.savefig(etiology_count_rel, bbox_inches="tight")


# plot sex distribution
plt.rcParams.update({'font.size':20})
study_organ_df = sex_df[['organ','study']]
fig, ax = plt.subplots(figsize = (11, 4))
plot_df = sex_df.pivot(index = 'study', columns = 'sex', values='relative').loc[all_studies_for_cols, :]* 100
colors = ['#4D6504', '#EB753E', 'lightgrey'] ## Colors for the bars
plot_df.plot(kind='bar', stacked=True, color=colors, ax = ax, width = 0.9) ## Plot
# Adding organ color bar or markers
organ_colors_bar = [org_colors.get(study_organ_df[study_organ_df['study'] == study]['organ'].values[0], 'lightgrey') for study in plot_df.index]
y_max = ax.get_ylim()[1]
ax.scatter(x=np.arange(len(plot_df.index)), y=[y_max] * len(plot_df.index), c=organ_colors_bar, s=100, marker='s')
ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
ax.set_ylabel('relative count [%]')
plt.savefig(sex_distribution, bbox_inches="tight")


fig, ax = plt.subplots(1, figsize = (11,4))

sns.violinplot(
    data = meta_df,
    x = meta_df['study'],
    y = 'age',
    hue = 'cond_test',
    palette = category_colors,
    width = 1.2,
    split = True,
    dodge = True,
    hue_order = ['control', 'fibrosis'],
    inner = 'stick'
)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
y_max = ax.get_ylim()[1]
ax.scatter(x=np.arange(len(plot_df.index)), y=[y_max] * len(plot_df.index), c=organ_colors_bar, s=100, marker='s')

plt.savefig(age_distribution, bbox_inches="tight")