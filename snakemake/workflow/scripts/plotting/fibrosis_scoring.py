import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

# Load configuration and input files from Snakemake
organs = snakemake.params.organs
organ_colors = snakemake.params.organ_colors
real_names = snakemake.params.organ_names
condition_colors = snakemake.params.condition_colors
geneset_file = snakemake.input.geneset

# Set plot style and backend
plt.rcParams.update({'font.size': 20})
plt.switch_backend('Agg')

# Load genesets and metadata
full_genesets = pd.read_csv(geneset_file, index_col=0)
all_sets = list(full_genesets['geneset'].unique())
metadata_cols = ["identifier", "grouping", "study", 'organ', "cond_test", 'sample', 'region',
                 'sex', 'batch', 'tech', 'age', 'fibrosis score (interstitial fibrosis) in %',
                 'ischemia time in sec', 'LVEF', 'BMI', 'Trichrome % fibrotic', 'modality', 'Fibrosis',
                 'eGFR', 'IFTA_binary', 'IFTA', 'clin_A1C', 'clin_global_glomerulosclerosis',
                 'clin_arteriosclerosis', 'hypertension', 'SAF.Score', 'Steatosis', 'Ballooning',
                 'Inflammation', 'Diabetes.type.2']



def add_sig_bracket(ax, y1, y2, x, h, text="*"):
    # added function to add sig. brackets
    ax.plot(
        [x, x + h, x + h, x],
        [y1, y1, y2, y2],
        color="black",
        linewidth=1
    )
    ax.text(
        x + h * 1.1,
        (y1 + y2) / 2,
        text,
        va="center",
        ha="left",
        fontsize=12
    )


def get_non_overlapping_level(y1, y2, occupied):
    """
    Find the first available bracket level where this y-range
    does not overlap with another bracket on the same level.
    Helper function for plotting sig. levels.
    """
    level = 0

    while True:
        collision = False

        for oy1, oy2, olevel in occupied:
            same_level = olevel == level
            overlaps = not (y2 < oy1 or y1 > oy2)

            if same_level and overlaps:
                collision = True
                break

        if not collision:
            return level

        level += 1



# Load ULM estimation and p-value results for all organs/studies
result_dict = {}
pval_dict = {}
for ulmest_file, ulmpval_file in zip(snakemake.input.ulmest, snakemake.input.ulmpval):
    ulm_est = pd.read_csv(ulmest_file).rename(columns={'Unnamed: 0': 'patient'})
    organ = ulmest_file.split("/")[-2]
    cols_to_keep = list(set(all_sets + metadata_cols) & set(ulm_est.columns))
    ulm_est = ulm_est[cols_to_keep]
    ulm_est['organ'] = organ

    pval = pd.read_csv(ulmpval_file).rename(columns={'Unnamed: 0': 'patient'})
    result_dict[ulm_est['study'].iloc[0]] = ulm_est
    pval_dict[pval['study'].iloc[0]] = pval
    pval['organ'] = organ



full_results_df = pd.concat(result_dict, ignore_index=True)
full_pval_df = pd.concat(pval_dict, ignore_index=True)

# Compute control means and standard deviations for normalization
group_columns = list(set(all_sets) & set(full_results_df.columns))
grouped = full_results_df.loc[:, group_columns + ['cond_test', 'study']].groupby(['study', 'cond_test']).mean().reset_index()
study_ctrl_means = grouped[grouped['cond_test'] == 'control']
merged_df = pd.merge(full_results_df, study_ctrl_means, on='study', suffixes=('', '_mean'))

grouped_std = full_results_df.loc[:, group_columns + ['cond_test', 'study']].groupby(['study', 'cond_test']).std().apply(abs).reset_index()
study_ctrl_std = grouped_std[grouped_std['cond_test'] == 'control']
merged_df = pd.merge(merged_df, study_ctrl_std, on='study', suffixes=('', '_std'))

# Calculate normalized differences from control
diff_columns = []
stdv_columns = []
for col in group_columns:
    merged_df[col + '_diff'] = merged_df[col] - merged_df[col + '_mean']
    diff_columns.append(col + '_diff')
    merged_df[col + '_stdev'] = (merged_df[col] - merged_df[col + '_mean']) / merged_df[col + '_std']
    stdv_columns.append(col + '_stdev')

# Define genesets to plot
genesets = ['NABA_CORE_MATRISOME', 'HALLMARK_INFLAMMATORY_RESPONSE', 'kidney', 'heart', 'lung', 'liver']
genesets_std = [i + '_stdev' for i in genesets]

# ============ PLOT 1: Sample counts per study heatmap ============
for organ in organs:
    sub = merged_df[merged_df['organ'] == organ]
    counts = pd.crosstab(sub["study"], sub["grouping"])
    percentages = counts.div(counts.sum(axis=1), axis=0) * 100

    fig, ax = plt.subplots(figsize=(sub['grouping'].nunique(), sub['study'].nunique()))
    sns.heatmap(percentages, cmap="Blues", annot=counts, fmt="d",
                cbar_kws={"label": "% of samples per study"}, linewidths=0.5, linecolor="white", ax=ax)
    ax.set_xlabel("disease etiology")
    ax.set_ylabel("study")
    ax.set_title("Sample counts per study")
    plt.tight_layout()
    plt.savefig(snakemake.output.sample_counts, dpi=300, bbox_inches='tight')
    plt.close()

# ============ PLOT 2: Geneset boxplots (raw values) ============
hue_order = ["control", "fibrosis"]
for idx, value in enumerate(genesets):
    fig, ax = plt.subplots(1, 4, figsize=(35, 10), sharey=True, tight_layout=True)
    for count, organ in enumerate(organs):
        data_fibrosis = (
            merged_df[merged_df["organ"] == organ]
            .loc[:, [value, "study", "cond_test"]]
            .melt(id_vars=["study", "cond_test"], var_name="variable", value_name="value")
        )

        sns.boxplot(data=data_fibrosis, x="study", y="value", hue="cond_test", ax=ax[count],
                    hue_order=hue_order, palette=condition_colors)
        sns.stripplot(data=data_fibrosis, x="study", y="value", hue="cond_test", ax=ax[count],
                      dodge=True, edgecolor="black", linewidth=0.2, legend=False,
                      palette=condition_colors, hue_order=hue_order)

        ax[count].set_title(f'{value} score \n in {organ}')
        ax[count].set_xlabel("study")
        ax[count].set_ylabel("ULM enrichment")
        ax[count].legend(loc='upper left', bbox_to_anchor=(1, 1))
        ax[count].set_xticklabels(ax[count].get_xticklabels(), rotation=45, ha='right')

    fig.suptitle(f"{value} fibrosis compared to control")
    plt.savefig(snakemake.output.geneset_boxplots[idx], dpi=300, bbox_inches='tight')
    plt.close()

# ============ PLOT 3: Standardized geneset boxplots ============
for idx, value in enumerate(genesets_std):
    fig, ax = plt.subplots(1, 4, figsize=(35, 10), sharey=False, tight_layout=True)
    for count, organ in enumerate(organs):
        data_fibrosis = (
            merged_df[merged_df["organ"] == organ]
            .loc[:, [value, "study", "cond_test"]]
            .melt(id_vars=["study", "cond_test"], var_name="variable", value_name="value")
        )

        sns.boxplot(data=data_fibrosis, x="study", y="value", hue="cond_test", ax=ax[count],
                    hue_order=hue_order, palette=condition_colors)
        sns.stripplot(data=data_fibrosis, x="study", y="value", hue="cond_test", ax=ax[count],
                      dodge=True, edgecolor="black", linewidth=0.2, legend=False,
                      palette=condition_colors, hue_order=hue_order)

        ax[count].set_title(f'{value} score \n in {organ}')
        ax[count].set_xlabel("study")
        ax[count].set_ylabel("ULM enrichment (standardized)")
        ax[count].legend(loc='upper left', bbox_to_anchor=(1, 1))
        ax[count].set_xticklabels(ax[count].get_xticklabels(), rotation=45, ha='right')
        ax[count].hlines(y=0, xmin=-0.5, xmax=5, linewidth=2, color='black')

    fig.suptitle(f"{value} fibrosis compared to control")
    plt.savefig(snakemake.output.geneset_std_boxplots[idx], dpi=300, bbox_inches='tight')
    plt.close()

# ============ PLOT 4: Organ score vs ECM score scatter plots ============
for idx, organ in enumerate(organs):
    merged_df_sub = merged_df[merged_df['organ'] == organ]
    conditions = merged_df_sub["grouping"].unique()
    palette = {"control": "darkgrey"}
    other_conditions = [c for c in conditions if c != "control"]
    colors = sns.color_palette("tab10", len(other_conditions))
    palette.update(dict(zip(other_conditions, colors)))

    fig, ax = plt.subplots(1, figsize=[5, 5])
    sns.scatterplot(ax=ax, data=merged_df_sub, x=f'{real_names[organ]}_stdev', y='NABA_CORE_MATRISOME_stdev',
                    hue='grouping', palette=palette)

    ax.hlines(y=0, xmin=merged_df_sub[f'{real_names[organ]}_stdev'].min() - 0.5,
              xmax=merged_df_sub[f'{real_names[organ]}_stdev'].max() + 0.5, linewidth=2, color='grey')
    ax.vlines(x=0, ymin=merged_df_sub['NABA_CORE_MATRISOME_stdev'].min() - 0.5,
              ymax=merged_df_sub['NABA_CORE_MATRISOME_stdev'].max() + 0.5, linewidth=2, color='grey')
    ax.set_xlabel(f'{real_names[organ]} fibrosis score')
    ax.set_ylabel('ECM score')
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig(snakemake.output.scatter_plots[idx], dpi=300, bbox_inches='tight')
    plt.close()

# ============ PLOT 5-8: Study-specific phenotype boxplots ============
study_specs = [
    ('McCown_2025_sn', 'Fibrosis', 'Fibrosis score', [-1, 10, 25, 50, 100], ["0-10%", "11-25%", "26-50%", "51-100%"]),
    ('Gribben_2024', 'Fibrosis', 'Fibrosis stage', None, [0, 1, 2, 3, 4]),
    ('Wilson_2022', 'IFTA', 'IFTA', None, {"0\xa0%": "0-10%", "1-10%": "0-10%"}),
    ('Simonson_2023', 'Trichrome % fibrotic', 'Trichrome \n % fibrotic', None, ['low', 'high']),
]

all_plot_points = []
all_pvalues = []

for plot_idx, (study, phenotype, ylabel, bins, labels_or_mapping) in enumerate(study_specs):
    merged_df_sub = merged_df[merged_df['study'] == study].copy()

    if len(merged_df_sub) == 0:
        continue

    organ = merged_df_sub['organ'].iloc[0]
    x_col = f'{real_names[organ]}_stdev'

    # Apply phenotype transformation if needed
    if isinstance(labels_or_mapping, dict):
        merged_df_sub.loc[:, phenotype] = (
            merged_df_sub[phenotype]
            .replace(labels_or_mapping)
        )

    elif bins is not None:
        merged_df_sub.loc[:, phenotype] = pd.cut(
            merged_df_sub[phenotype],
            bins=bins,
            labels=labels_or_mapping
        )

    elif phenotype == 'Fibrosis' and study == 'Gribben_2024':
        merged_df_sub.loc[:, phenotype] = (
            merged_df_sub.loc[:, phenotype]
            .astype('int')
            .astype('category')
        )

    # Define y-axis order
    if study == 'Wilson_2022':
        order = list(reversed(merged_df_sub[phenotype].dropna().unique()))
    elif isinstance(labels_or_mapping, list):
        order = labels_or_mapping
    else:
        order = list(merged_df_sub[phenotype].dropna().unique())

    # Save plotted data points
    plot_points = merged_df_sub[[phenotype, x_col]].copy()
    plot_points["study"] = study
    plot_points["organ"] = real_names[organ]
    plot_points["phenotype"] = phenotype
    plot_points["category"] = plot_points[phenotype]
    plot_points["x_col"] = x_col
    plot_points["x_value"] = plot_points[x_col]
    plot_points["category_order"] = plot_points["category"].apply(
        lambda x: order.index(x) if x in order else pd.NA
    )

    all_plot_points.append(
        plot_points[
            [
                "study",
                "organ",
                "phenotype",
                "category",
                "category_order",
                "x_col",
                "x_value",
            ]
        ]
    )

    fig, ax = plt.subplots(1, figsize=[5, 3], tight_layout=True)

    sns.boxplot(
        ax=ax,
        data=merged_df_sub,
        x=x_col,
        y=phenotype,
        color=organ_colors[organ],
        order=order
    )

    # Pairwise one-sided Mann-Whitney U tests
    # cat1 is always the smaller / earlier y-axis category.
    pairs = [
        (order[i], order[j])
        for i in range(len(order))
        for j in range(i + 1, len(order))
    ]

    x_min, x_max = ax.get_xlim()
    x_range = x_max - x_min

    bracket_h = 0.03 * x_range
    level_gap = 0.08 * x_range
    occupied = []
    max_level = -1

    for cat1, cat2 in pairs:
        vals1 = merged_df_sub.loc[
            merged_df_sub[phenotype] == cat1, x_col
        ].dropna()

        vals2 = merged_df_sub.loc[
            merged_df_sub[phenotype] == cat2, x_col
        ].dropna()

        stat, pval = mannwhitneyu(
            vals1,
            vals2,
            alternative="less"
        )

        is_significant = pval < 0.05

        all_pvalues.append({
            "study": study,
            "organ": real_names[organ],
            "phenotype": phenotype,
            "x_col": x_col,
            "category_1": cat1,
            "category_2": cat2,
            "n_category_1": len(vals1),
            "n_category_2": len(vals2),
            "pvalue": pval,
            "significant_p_lt_0_05": is_significant,
        })

        if is_significant:
            y1 = order.index(cat1)
            y2 = order.index(cat2)

            level = get_non_overlapping_level(y1, y2, occupied)
            occupied.append((y1, y2, level))
            max_level = max(max_level, level)

            bracket_x = (
                x_max
                + 0.03 * x_range
                + level * level_gap
            )

            add_sig_bracket(
                ax=ax,
                y1=y1,
                y2=y2,
                x=bracket_x,
                h=bracket_h,
                text="*"
            )

    # Expand x-axis so all brackets are visible
    if max_level >= 0:
        ax.set_xlim(
            x_min,
            x_max + (0.15 + (max_level + 1) * 0.08) * x_range
        )

    ax.set_xlabel(f'{real_names[organ]} fibrosis score')
    ax.set_ylabel(ylabel)
    ax.set_title(study)

    plt.savefig(
        snakemake.output.study_plots[plot_idx],
        dpi=300,
        bbox_inches='tight'
    )




# Save all data points and all p-values
all_plot_points_df = pd.concat(all_plot_points, ignore_index=True)
all_pvalues_df = pd.DataFrame(all_pvalues)

all_plot_points_df.to_csv(snakemake.output.study_plots_data, index=False)
all_pvalues_df.to_csv(snakemake.output.study_plots_pval, index=False)





# ============ PLOT 9: Combined organ score boxplot by disease grouping ============
extra_column = "grouping"
n_categories = []
for organ in organs:
    categories = (
        merged_df.loc[merged_df["organ"] == organ, extra_column]
        .dropna().unique()
    )
    n_categories.append(len(categories))

fig, ax = plt.subplots(1, len(organs), figsize=(sum(n_categories) * 0.8, 8),
                        sharey=True, gridspec_kw={"width_ratios": n_categories})

for count, organ in enumerate(organs):
    diff_columns = f"{real_names[organ]}_stdev"
    data_fibrosis = merged_df[merged_df["organ"] == organ].loc[:, [diff_columns, "study", extra_column]]
    categories = data_fibrosis[extra_column].dropna().unique().tolist()
    order = ["control"] + sorted([x for x in categories if x != "control"])

    sns.boxplot(data=data_fibrosis, x=extra_column, y=diff_columns, order=order,
                ax=ax[count], color=organ_colors[organ])
    sns.boxplot(data=data_fibrosis[data_fibrosis[extra_column] == "control"],
                x=extra_column, y=diff_columns, order=order, ax=ax[count], color="darkgrey")

    ax[count].set_xticklabels(ax[count].get_xticklabels(), rotation=45, ha='right')
    ax[count].set_title(f"{real_names[organ]} score")
    ax[count].set_ylabel("normalized organ score")
    ax[count].set_xlabel("")
    ax[count].hlines(0, -0.5, len(order) - 0.5, color="black")

plt.tight_layout()
plt.savefig(snakemake.output.combined_boxplot, dpi=300, bbox_inches='tight')
plt.close()
