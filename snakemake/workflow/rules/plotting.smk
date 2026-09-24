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

extra_comparison_pairs = [
    ('mesenchymal', 'heart', 'lung'),
    ('endothelial', 'kidney', 'liver'),
    ('endothelial', 'lung', 'heart'),
    ('endothelial', 'liver', 'lung'),
    ('epithelial', 'lung', 'kidney'),
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
        umap2 = expand('plots/analysis/myofib_identification/{organ}_umap2.pdf', organ=meta_organs),
        proportion = expand('plots/analysis/myofib_identification/{organ}_proportion_stats.pdf', organ=meta_organs),
        top5 = expand('plots/analysis/myofib_identification/{organ}_top5.pdf', organ=meta_organs),
        top5_std = expand('plots/analysis/myofib_identification/{organ}_top5_standard_scale.pdf', organ=meta_organs),
        pan = expand('plots/analysis/myofib_identification/{organ}_pan_markers.pdf', organ=meta_organs),
        pan_std = expand('plots/analysis/myofib_identification/{organ}_pan_markers_standard_scale.pdf', organ=meta_organs),
        top5_all = expand('plots/analysis/myofib_identification/{organ}_top5_all.pdf', organ=meta_organs),
        top5_std_all = expand('plots/analysis/myofib_identification/{organ}_top5_standard_scale_all.pdf', organ=meta_organs),
        pan_all = expand('plots/analysis/myofib_identification/{organ}_pan_markers_all.pdf', organ=meta_organs),
        pan_std_all = expand('plots/analysis/myofib_identification/{organ}_pan_markers_standard_scale_all.pdf', organ=meta_organs),
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
        #scArches_ecdf = 'results/integration/{organ}/post_label_transfer/{study}_scarches_model/scArches_ecdf.pdf',
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
        r2 = expand('plots/analysis/qc/r2/mofacell_R2_{organ}.csv', organ = meta_organs),
        r2plot = 'plots/analysis/qc/r2/mofacell_summedR2.pdf',
        factor_corr_plot_path = 'plots/analysis/qc/r2/mofacell_weights_corr.pdf',
        topweight_plot_path = 'plots/analysis/qc/r2/top_genes_weights.pdf',
        topweight_csv_path = 'plots/analysis/qc/r2/top_genes_weights.csv',
        geneweights = expand("results/models/shared_factor_genes/{organ}_geneweights.csv", organ=meta_organs),
    params:
        organs = meta_organs,
        organ_names = config['general_plotting']['organ_names'],
        organ_colors = config['general_plotting']['organ_colors'],
        condition_colors = config['general_plotting']['condition_colors'],
        views = config["general_plotting"]["views"],
    conda:
        "../envs/scanpy.yaml"
    resources:
        mem_mb = 2000
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


rule circosplots:
    input:
        ccc_results = "results/ccc/on_dl_allorgans.pckl"
    output:
        plots = expand(
            "plots/analysis/ccc/circosplots/{organ}.pdf",
            organ=meta_organs
        ),
        ccc_table = "results/ccc/on_dl_allorgans.csv"
    params:
        organs = meta_organs,
        real_names = config["general_plotting"]["organ_names"],
        ctype_colors = config["general_plotting"]["celltype_colors_human"],
        eff_cutoff_pos = 0.5,
        top_n = 70
    conda:
        "../envs/ccc_cell2location.yaml"
    resources:
        mem_mb = 5000,
        runtime = 60
    script:
        "../scripts/plotting/circosplots.py"


rule fibrosis_scoring:
    input:
        ulmest = expand('results/patient_scoring/{placeholder[0]}/{placeholder[1]}_ulmest.csv', placeholder = expanded_list),
        ulmpval = expand('results/patient_scoring/{placeholder[0]}/{placeholder[1]}_ulmpval.csv', placeholder = expanded_list),
        geneset = "results/patient_scoring/pat_scoring_geneset.csv"
    output:
        sample_counts = 'plots/analysis/fibrosis_scoring/sample_counts.pdf',
        geneset_boxplots = expand('plots/analysis/fibrosis_scoring/geneset_{gs}_boxplot.pdf',
                                  gs=['NABA_CORE_MATRISOME', 'HALLMARK_INFLAMMATORY_RESPONSE', 'kidney', 'heart', 'lung', 'liver']),
        geneset_std_boxplots = expand('plots/analysis/fibrosis_scoring/geneset_{gs}_stdev_boxplot.pdf',
                                      gs=['NABA_CORE_MATRISOME', 'HALLMARK_INFLAMMATORY_RESPONSE', 'kidney', 'heart', 'lung', 'liver']),
        scatter_plots = expand('plots/analysis/fibrosis_scoring/{organ}_ecm_scatter.pdf', organ = meta_organs),
        study_plots = expand('plots/analysis/fibrosis_scoring/{study_name}.pdf',
                            study_name=['McCown_2025_sn', 'Gribben_2024', 'Wilson_2022', 'Simonson_2023']),
        combined_boxplot = 'plots/analysis/fibrosis_scoring/combined_organ_boxplot.pdf',
        study_plots_data ='plots/analysis/fibrosis_scoring/study_fibrosis_score_data.csv',
        study_plots_pval ='plots/analysis/fibrosis_scoring/study_fibrosis_score_pval.csv'
    params:
        organs = meta_organs,
        organ_colors = config['general_plotting']['organ_colors'],
        organ_names = config['general_plotting']['organ_names'],
        condition_colors = config['general_plotting']['condition_colors'],
        output_dir = 'plots/analysis/fibrosis_scoring'
    conda:
        "../envs/scanpy.yaml"
    resources:
        mem_mb = 8000,
        runtime = 60
    script:
        "../scripts/plotting/fibrosis_scoring.py"


rule pbulk_stats:
    input:
        pseudobulks = expand('results/zenodo/pseudobulks/{placeholder[0]}/{placeholder[1]}.h5ad', placeholder = meta_expanded_list)
    output:
        pdfs = expand('plots/analysis/pbulkstats/{organ}_pbulk.pdf', organ = meta_organs)
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/plotting/pbulk_stats.py"


rule compare_mofa_lmm:
    input:
        dl_results = "results/dl_meta/cross_organ_dl.pckl",
        mofa_weights = expand(
            "results/models/shared_factor_genes/{organ}_geneweights.csv",
            organ=meta_organs
        )
    output:
        scatterplots = "plots/analysis/interorgan_comparison/compare_mofa_lmm/mofa_lmm_scatterplots.pdf",
        heatmap = "plots/analysis/interorgan_comparison/compare_mofa_lmm/mofa_lmm_correlation_heatmap.pdf",
        corr_df = "plots/analysis/interorgan_comparison/compare_mofa_lmm/corr_df.csv"
    params:
        organs = meta_organs,
        views = config["general_plotting"]["views"],
        organ_names = config["general_plotting"]["organ_names"]
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/plotting/compare_mofa_lmm.py"


rule extra_comparisons:
    input:
        cross_organ_dl = "results/dl_meta/cross_organ_dl.pckl"
    output:
        plots = expand(
            "plots/analysis/interorgan_comparison/extra_comparisons/{celltype}_{organ1}_{organ2}.pdf",
            zip,
            celltype=[comparison[0] for comparison in extra_comparison_pairs],
            organ1=[comparison[1] for comparison in extra_comparison_pairs],
            organ2=[comparison[2] for comparison in extra_comparison_pairs]
        ),
        data = expand(
            "plots/analysis/interorgan_comparison/extra_comparisons/{celltype}_{organ1}_{organ2}.csv",
            zip,
            celltype=[comparison[0] for comparison in extra_comparison_pairs],
            organ1=[comparison[1] for comparison in extra_comparison_pairs],
            organ2=[comparison[2] for comparison in extra_comparison_pairs]
        )
    params:
        comparisons = extra_comparison_pairs,
        n_top = 10,
        effect_threshold = 0.5,
        organs = list(config["general_plotting"]["organ_names"].values()),
        organ_colors = config["general_plotting"]["organ_colors"],
        organ_names = config["general_plotting"]["organ_names"]
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/plotting/extra_comparisons.py"
