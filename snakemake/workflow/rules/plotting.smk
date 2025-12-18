import pandas as pd

wildcard_constraints:
    organ='[a-zA-Z]+'


dataset_combinations = config['datasets']
human_studies = ['HCAlung', 'reheatHeart', 'kidney', 'liver']
meta_organs = config['meta_organs']

expanded_list = []

for key, items in dataset_combinations.items():
    if key in human_studies:
        for item in items['studies_ownmodel']:
            expanded_list.append([key, item])

meta_expanded_list = [
    (organ, study)
    for organ in meta_organs
    for study in dataset_combinations[organ]['studies_to_analyze']
]


       
rule input_overlap:
    input:
        models = expand('results/models/{placeholder[0]}/{placeholder[1]}_mofa.hdf5', placeholder = expanded_list),
        markers = expand('results/mofa_input/{placeholder[0]}/mrkrs/{placeholder[1]}.csv', placeholder = expanded_list),
    output:
        jaccard_marker = 'plots/analysis/interorgan_comparison/marker_overlap_jaccard.pdf',
        abs_marker ='plots/analysis/interorgan_comparison/marker_overlap_absolute.pdf',
        jaccard_feat = 'plots/analysis/interorgan_comparison/model_genes_overlap_jaccard.pdf',
        abs_feat = 'plots/analysis/interorgan_comparison/model_genes_overlap_absolute.pdf',
        jaccard_marker_csv = 'plots/analysis/interorgan_comparison/marker_overlap_jaccard.csv',
        abs_marker_csv ='plots/analysis/interorgan_comparison/marker_overlap_absolute.csv',
        jaccard_feat_csv = 'plots/analysis/interorgan_comparison/model_genes_overlap_jaccard.csv',
        abs_feat_csv = 'plots/analysis/interorgan_comparison/model_genes_overlap_absolute.csv',
    params:
        celltype_colors = config['general_plotting']['celltype_colors_human'],
        organ_colors = config['general_plotting']['organ_colors'],
        organs = config['meta_organs'],
        pairing = expanded_list
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/plotting/input_overlap.py"
        

       
rule qc_plots:
    input:
        qc_df = expand('results/preprocessing/qc/{placeholder[0]}/{placeholder[1]}.csv', placeholder = expanded_list),
        metadata_df = expand('results/mofa_input/{placeholder[0]}/metadata/{placeholder[1]}.csv', placeholder = expanded_list),
        ctypecount_df = expand('results/preprocessing/qc/{placeholder[0]}/{placeholder[1]}_ctypecount.csv', placeholder = expanded_list),
        ctypecount_patient_df = expand('results/preprocessing/qc/{placeholder[0]}/{placeholder[1]}_ctypecount_patient.csv', placeholder = expanded_list),
    output:
        qc_pdf = 'plots/analysis/qc/qc_overview.pdf',
        metadata = 'plots/analysis/qc/sample_no.pdf',
        ctype_count_abs = 'plots/analysis/qc/ctype_count_abs.pdf',
        ctype_count_rel = 'plots/analysis/qc/ctype_count_rel.pdf',
        plot_df_abs_csv = 'plots/analysis/qc/ctype_count_abs.csv',
        plot_df_rel_csv = 'plots/analysis/qc/ctype_count_rel.csv',
        ctype_count_change = 'plots/analysis/qc/ctype_count_change.pdf',
        etiology_count = 'plots/analysis/qc/etiology_count.pdf',
        etiology_count_rel = 'plots/analysis/qc/etiology_count_rel.pdf',
        t_matrix_all = 'plots/analysis/qc/ctype_count_change_tvals.csv',
        adj_p_matrix_all = 'plots/analysis/qc/ctype_count_change_pvals.csv',
        lm_all = 'plots/analysis/qc/ctype_count_change_lm.csv',
        sex_distribution = 'plots/analysis/qc/gender_count.pdf',
        age_distribution = 'plots/analysis/qc/age_dist_studies.pdf'
    params:
        organ_colors = config['general_plotting']['organ_colors'],
        organ_colors_light = config['general_plotting']['organ_colors_light'],
        pairing = expanded_list,
        ctype_colors = config['general_plotting']['celltype_colors_human'],
        views = config['general_plotting']['views'],
        real_names = config['general_plotting']['organ_names'],
        study_sub = config['datasets']['kidney']['studies_to_analyze']
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/plotting/qc_plots.py"
        

rule myofib_identification_plots:
    input:
        fibroblast = 'results/preprocessing/fibroblasts/all_fibs.h5ad',
        scvi = 'results/integration/fibroblasts/integrated_latent_space.csv',
        markers = expand('results/preprocessing/fibroblasts/{organ}_markers_scvi.csv', organ=meta_organs),
    output:
        stacked_bar = 'plots/analysis/myofib_identification/substate_stacked_bar.pdf',
        umap = expand('plots/analysis/myofib_identification/{organ}_umap.pdf', organ=meta_organs),
        proportion = expand('plots/analysis/myofib_identification/{organ}_proportion_stats.pdf', organ=meta_organs),
        top5 = expand('plots/analysis/myofib_identification/{organ}_top5.pdf', organ=meta_organs),
        top5_std = expand('plots/analysis/myofib_identification/{organ}_top5_standard_scale.pdf', organ=meta_organs),
        pan = expand('plots/analysis/myofib_identification/{organ}_pan_markers.pdf', organ=meta_organs),
        pan_std = expand('plots/analysis/myofib_identification/{organ}_pan_markers_standard_scale.pdf', organ=meta_organs),
    resources:
        runtime=60,
        mem_mb = 120000
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/plotting/myofib_identification_plots.py"


# Plot predictions of scArches

rule plot_scArches:
    input:
        ref = 'results/integration/{organ}/with_embeddings.h5ad',
        query = 'results/integration/{organ}/post_label_transfer/{study}_scarches_model/adata.h5ad',
    output:
        confusion_matrix = 'results/integration/{organ}/post_label_transfer/{study}_scarches_model/confusion_matrix.pdf',
        scArches_ecdf = 'results/integration/{organ}/post_label_transfer/{study}_scarches_model/scArches_ecdf.pdf',
        scArches_umap = 'results/integration/{organ}/post_label_transfer/{study}_scarches_model/scArches_umap.pdf'
    resources:
        mem_mb = 30000
    conda:
        "../envs/scanpy.yaml"
    script:
        '../scripts/plotting/plot_integration.py'



rule plot_scDist:
    input:
        scdist_results = expand('results/scDist/{placeholder[0]}/{placeholder[1]}_results.csv', placeholder = expanded_list)
    output:
        scdist_summary_csv = 'plots/analysis/interorgan_comparison/scDist/summary.csv',
        study_results =  'plots/analysis/interorgan_comparison/scDist/scDist_results_studies_sep.pdf',
        all_results = 'plots/analysis/interorgan_comparison/scDist/scDist_results.pdf',
    params:
        organs = config['meta_organs'],
        org_colors = config['general_plotting']['organ_colors'],
        real_names = config['general_plotting']['organ_names'],
        views = config["general_plotting"]["views"],
        datasets = config['datasets'],
        study_colors = config['general_plotting']['study_colors']
    conda:
        "../envs/scanpy.yaml"
    script:
        '../scripts/plotting/plot_scDist.py'


rule plot_mofa_fibrosis:
    input:
        models = expand('results/models/{organ}/meta_model_mofa.hdf5', organ = meta_organs),
        metadata = expand('results/mofa_input/{placeholder[0]}/metadata/{placeholder[1]}.csv', placeholder = expanded_list)
    output:
        scatter = 'plots/analysis/interorgan_comparison/mofacell_fibrosis_factor_scatter.pdf',
        boxplot = 'plots/analysis/interorgan_comparison/mofacell_fibrosis_factor_boxplot_study.pdf',
        r2 = expand('plots/analysis/qc/r2/mofacell_R2_{organ}.csv', organ = meta_organs)
    params:
        organs = meta_organs,
        organ_names = config['general_plotting']['organ_names'],
        organ_colors = config['general_plotting']['organ_colors'],
        condition_colors = config['general_plotting']['condition_colors']
    conda:
        "../envs/scanpy.yaml"
    resources:
        mem_mb = 5000
    script:
        "../scripts/plotting/mofa_plot.py"


rule ccc_dl:
    input:
        ccc_results = "results/ccc/on_dl_allorgans.pckl"
    output:
        plots = "plots/analysis/ccc/ccc_dl_overview.pdf"
    conda:
        "../envs/scanpy.yaml"
    params:
        eff_cutoff_pos = 0.5,
        eff_cutoff_neg = -0.5
    resources:
        mem_mb = 5000,
        runtime=60
    script:
        "../scripts/plotting/ccc_dl.py"
