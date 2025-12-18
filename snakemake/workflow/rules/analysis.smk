import pandas as pd

wildcard_constraints:
    organ='[a-zA-Z]+'


script_dir = config['scriptdir']
dataset_combinations = config['datasets']

human_studies = ['HCAlung', 'reheatHeart', 'kidney', 'liver']


study_list = [(organ, study) for organ, studies_dict in dataset_combinations.items() for study in studies_dict['studies_to_analyze']]
all_studies = pd.DataFrame(study_list, columns=['organ', 'study'])


expanded_list = []

for key, items in dataset_combinations.items():
    if key in human_studies:
        for item in items['studies_to_analyze']:
            expanded_list.append([key, item])

all_human_studies = all_studies[all_studies['organ'].isin(human_studies)]

organs = all_studies['organ'].unique()

# one model per study
rule fit_mofa_model:
    input:
        meta = lambda w: expand('results/mofa_input/{{organ}}/metadata/{study}.csv', study=all_studies[all_studies['organ'] == w.organ]['study']),
        pbulk = lambda w: expand('results/mofa_input/{{organ}}/pbulk/{study}.csv',  study=all_studies[all_studies['organ'] == w.organ]['study']),
        coldata = lambda w: expand('results/mofa_input/{{organ}}/coldata/{study}.csv', study=all_studies[all_studies['organ'] == w.organ]['study']),
        markers = lambda w: expand('results/mofa_input/{{organ}}/mrkrs/{study}.csv', study=all_studies[all_studies['organ'] == w.organ]['study']),
    output:
        plot = 'results/models/{organ}/model_stats.pdf',
        file = 'results/models/{organ}/model_stats.csv'
    params:
        aesthestics_R = script_dir + '/analysis/aesthetics.R',
        mofa_dir = 'results/models/{organ}/'
    conda:
        "../envs/mofa.yaml"
    script:
        "../scripts/analysis/fit_mofa_model.R"

# one model per organ       
rule mofa_meta:
    input:
        meta = lambda w: expand('results/mofa_input/{{organ}}/metadata/{study}.csv', study = all_studies[all_studies['organ'] == w.organ]['study']),
        pbulk = lambda w: expand('results/mofa_input/{{organ}}/pbulk/{study}.csv', study =  all_studies[all_studies['organ'] == w.organ]['study']),
        coldata = lambda w: expand('results/mofa_input/{{organ}}/coldata/{study}.csv',study = all_studies[all_studies['organ'] == w.organ]['study']),
        markers = lambda w: expand('results/mofa_input/{{organ}}/mrkrs/{study}.csv',study = all_studies[all_studies['organ'] == w.organ]['study']),
    output:
        umap_pdf = 'results/models/{organ}/meta_model_umap.pdf',
        file = 'results/models/{organ}/meta_model_mofa.hdf5',
        mofa_pdf = 'results/models/{organ}/meta_model_mofa.pdf',
        stats_pdf = 'results/models/{organ}/meta_model_stats.pdf',
        factor_pdf = 'results/models/{organ}/meta_model_factors.pdf'
    params:
        aesthestics_R = script_dir + '/analysis/aesthetics.R',
        mofa_dir = 'results/models/{organ}/'
    conda:
        "../envs/mofa.yaml"
    script:
        "../scripts/analysis/mofa_meta.R"

  
rule fit_meta_organ_model:
    input:
        meta = expand('results/mofa_input/{placeholder[0]}/metadata/{placeholder[1]}.csv', placeholder = expanded_list),
        pbulk = expand('results/mofa_input/{placeholder[0]}/pbulk/{placeholder[1]}.csv',  placeholder = expanded_list),
        coldata = expand('results/mofa_input/{placeholder[0]}/coldata/{placeholder[1]}.csv', placeholder = expanded_list),
        markers = expand('results/mofa_input/{placeholder[0]}/mrkrs/{placeholder[1]}.csv', placeholder = expanded_list),
    output:
        umap_pdf = 'results/models/organ_comp/meta_model_umap.pdf',
        file = 'results/models/organ_comp/meta_model_mofa.hdf5',
        mofa_pdf = 'results/models/organ_comp/meta_model_mofa.pdf',
        stats_pdf = 'results/models/organ_comp/meta_model_stats.pdf',
        factor_pdf = 'results/models/organ_comp/meta_model_factors.pdf'
    params:
        aesthestics_R = script_dir + '/analysis/aesthetics.R',
        mofa_dir = 'results/models/organ_comp/'
    conda:
        "../envs/mofa.yaml"
    script:
        "../scripts/analysis/mofa_meta.R"
        
       
rule deg_perstudy:
    input:
        pbulk = 'results/mofa_input/{organ}/pbulk/{study}.csv',
        meta = 'results/mofa_input/{organ}/metadata/{study}.csv',
    output:
        deg_file = 'results/deg/{organ}/{study}_all_deg_tval.csv'
    params:
        config['preprocessing'].get('pseudobulk')
    conda:
        "../envs/scanpy.yaml"
    threads: 16
    script:
        "../scripts/analysis/deg_perstudy.py"
        
        
rule deg_analysis:
    input:
        deg_file = expand('results/deg/{placeholder[0]}/{placeholder[1]}_all_deg_tval.csv', placeholder = expanded_list),
    output:
        organ_spec_dot = 'plots/analysis/interorgan_comparison/organ_specific_degs.pdf',
        no_organ_consensus = 'plots/analysis/interorgan_comparison/no_organ_consensus.pdf',
        overall_consensus = 'plots/analysis/interorgan_comparison/overall_consensus.pdf',
        overlap_organ_consensus = 'plots/analysis/interorgan_comparison/upsetplot_deg_overlap.pdf',
        organ_consenus_up = 'results/deg/common/common_deg_up_perorgan.json',
        organ_consenus_down = 'results/deg/common/common_deg_down_perorgan.json',
        common_deg_file = 'results/deg/common/commonly_regulated_degs.json',
        common_deg_file_up = 'results/deg/common/commonly_regulated_degs_up.json',
        common_deg_file_down = 'results/deg/common/commonly_regulated_degs_down.json',
        combined_pval_organ = 'plots/analysis/interorgan_comparison/top500_degs_combinedpval_organ.pdf',
        deg_count = 'plots/analysis/interorgan_comparison/deg_nr.pdf',
        deg_sim_noclust = 'plots/analysis/interorgan_comparison/deg_dotproduct.pdf',
        deg_sim_clust = 'plots/analysis/interorgan_comparison/deg_dotproduct_clustered.pdf'
    params:
        celltype_colors = config['general_plotting']['celltype_colors_human'],
        organ_colors = config['general_plotting']['organ_colors'],
        study_colors = config['general_plotting']['study_colors'],
        organ_names = config['general_plotting']['organ_names'],
        views = config['general_plotting']['views']
    resources:
        runtime=180
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/analysis/deg_comparison.py"


rule deg_projection:
    input:
        deg_file = expand('results/deg/{placeholder[0]}/{placeholder[1]}_all_deg_tval.csv', placeholder = expanded_list),
        pb_file = expand('results/mofa_input/{placeholder[0]}/pbulk/{placeholder[1]}.csv', placeholder = expanded_list),
        metadata_file = expand('results/mofa_input/{placeholder[0]}/metadata/{placeholder[1]}.csv', placeholder = expanded_list),
        gene_dict_up_output = 'results/dl_meta/organ_spec_dl_up.pckl',
        gene_dict_down_output = 'results/dl_meta/organ_spec_dl_down.pckl'
    output:
        predicitions_auroc = 'results/deg/predicitions_auroc_{mode}.pkl',
        predicitions_clustered = 'plots/analysis/interorgan_comparison/prediction_deg_clustered_{mode}.pdf',
        predicitions_non_clustered = 'plots/analysis/interorgan_comparison/prediction_deg_{mode}.pdf',
        predicitions_box_withinorgan = 'plots/analysis/interorgan_comparison/prediction_deg_boxplot_within_organ_{mode}.pdf',
        predicitions_organ_agg = 'plots/analysis/interorgan_comparison/prediction_deg_aggregated_between_organ_{mode}.pdf',
        heatmap_output = 'plots/analysis/interorgan_comparison/prediction_deg_aggregated_between_organ_{mode}.pkl',
        stripplot_output = 'plots/analysis/interorgan_comparison/prediction_deg_boxplot_within_organ_{mode}.pkl'
    params:
        celltype_colors = config['general_plotting']['celltype_colors_human'],
        organ_colors = config['general_plotting']['organ_colors'],
        study_colors = config['general_plotting']['study_colors'],
        organ_names = config['general_plotting']['organ_names'],
        views = config['general_plotting']['views']
    wildcard_constraints:
        mode ="(all|unique)"
    conda:
        "../envs/scanpy.yaml"
    resources:
        runtime = 300
    script:
        "../scripts/analysis/deg_projection.py"


rule scDist_all:
    input:
        adata_files = 'results/preprocessing/{organ}/{study}.h5ad'
    output:
        feature_importance = 'results/scDist/{organ}/{study}_feature_importance.csv',
        results = 'results/scDist/{organ}/{study}_results.csv',
    resources:
         mem_mb=350000
    threads: 8
    conda:
        "../envs/scDist.yaml"
    script:
        "../scripts/analysis/scDist_all.R"

       
rule naba_score:
    input:
        pbulk = 'results/mofa_input/{organ}/pbulk/{study}.csv',
        meta = 'results/mofa_input/{organ}/metadata/{study}.csv',
    output:
        deg_file = 'results/naba/{organ}/{study}_naba_scores.csv'
    params:
        config['preprocessing'].get('pseudobulk'),
        'data/misc/genesets/naba/NABA_CORE_MATRISOME.v2023.2.Hs.gmt'
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/analysis/NABA_score.py"


rule patient_score:
    input:
        adata = 'results/preprocessing/{organ}/{study}.h5ad'
    output:
        out_path_est = 'results/patient_scoring/{organ}/{study}_ulmest.csv',
        out_path_pval = 'results/patient_scoring/{organ}/{study}_ulmpval.csv',
        pbulk = 'results/patient_scoring/{organ}/{study}_pbulk.csv'
    params:
        config['preprocessing'].get('pseudobulk'),
        'data/misc/genesets/patient_scoring_sets.csv'
    resources:
         mem_mb=120000
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/analysis/patient_scoring.py"


rule patient_score_part2:
    input:
        out_path_est = expand('results/patient_scoring/{placeholder[0]}/{placeholder[1]}_ulmest.csv', placeholder = expanded_list),
        out_path_pval = expand('results/patient_scoring/{placeholder[0]}/{placeholder[1]}_ulmpval.csv', placeholder = expanded_list)
    output:
        pbulk = 'results/patient_scoring/geneset_scores_full_pbulks.csv'
    params:
        config['preprocessing'].get('pseudobulk'),
        'data/misc/genesets/NABA_IMMUNE.csv'
    resources:
         mem_mb=120000
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/analysis/patient_scoring_part2.py"




rule fibroblasts:
    input:
        expand('results/preprocessing/{placeholder[0]}/{placeholder[1]}.h5ad', placeholder = expanded_list)
    output:
        adata = 'results/preprocessing/fibroblasts/all_fibs.h5ad'
    resources:
         mem_mb=120000
    threads: 8
    params:
        organs = config['meta_organs'],
        hvg = 2000,
        study_organ = all_human_studies
    conda:
        "../envs/harmony.yaml"
    script:
        "../scripts/analysis/fibroblast_pp.py"


rule organ_clusters_extract_markers:
    input:
        adata = 'results/preprocessing/fibroblasts/all_fibs.h5ad',
        scvi_results = 'results/integration/fibroblasts/integrated_latent_space.csv'
    output:
        marker_output_scvi = 'results/preprocessing/fibroblasts/{organ}_markers_scvi.csv',
        full_marker_output_scvi = 'results/preprocessing/fibroblasts/{organ}_full_markers_scvi.csv',
        pbulk_scvi = 'results/preprocessing/fibroblasts/{organ}_pbulk.csv',
    resources:
        mem_mb=180000,
        runtime=1000
    params:
        resolutions = config['preprocessing']['fibroblasts']['scvi_resolutions']
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/analysis/organ_clusters_extract_markers.py"



rule dl_meta:
    input:
        deg_file = expand('results/deg/{placeholder[0]}/{placeholder[1]}_all_deg_tval.csv', placeholder = expanded_list),
    output:
        organ_spec = 'results/dl_meta/organ_spec_dl.pckl',
        cross_organ = 'results/dl_meta/cross_organ_dl.pckl',
    params:
        organ_colors = config['general_plotting']['organ_colors'],
        organ_names = config['general_plotting']['organ_names'],
        views = config['general_plotting']['views']
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/analysis/dl_meta.py"


rule ccc_dl:
    input:
        organ_spec_dl = 'results/dl_meta/organ_spec_dl.pckl',
    output:
        ccc_result = 'results/ccc/on_dl_allorgans.pckl'
    params:
        organs = config['meta_organs'],
        views = config['general_plotting']['views']
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/analysis/ccc_on_dlresults.py"


rule analysis_dl:
    input:
        organ_spec_dl = 'results/dl_meta/organ_spec_dl.pckl',
        cross_organ_dl = 'results/dl_meta/cross_organ_dl.pckl',
        deg_file = expand('results/deg/{placeholder[0]}/{placeholder[1]}_all_deg_tval.csv', placeholder = expanded_list),
    output:
        cross_organ_per_ctype_dl = 'plots/analysis/interorgan_comparison/dl_meta/cross_organ_per_ctype_dl.pdf',
        upsetplot_dl_organspec_genes = 'plots/analysis/interorgan_comparison/dl_meta/upsetplot_dl_organspec_genes.pdf',
        gene_dict_up_output = 'results/dl_meta/organ_spec_dl_up.pckl',
        gene_dict_down_output = 'results/dl_meta/organ_spec_dl_down.pckl',
        count_organspec_genes = 'plots/analysis/interorgan_comparison/dl_meta/count_organspec_genes.pdf',
        counts_up = 'plots/analysis/interorgan_comparison/dl_meta/count_organspec_up.csv',
        counts_down = 'plots/analysis/interorgan_comparison/dl_meta/count_organspec_down.csv',
        topx = 'plots/analysis/interorgan_comparison/dl_meta/cross_organ_per_ctype_dl_topx.pdf',
    params:
        organ_colors = config['general_plotting']['organ_colors'],
        organ_names = config['general_plotting']['organ_names'],
        views = config['general_plotting']['views']
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/analysis/dl_meta_analysis_plot.py"


rule dl_meta_top_genes:
    input:
        organ_spec_dl = 'results/dl_meta/organ_spec_dl.pckl'
    output:
        top_genes_csv = 'results/dl_meta/organ_spec_top10_upregulated_allorgans.csv'
    params:
        organs = config['meta_organs'],
        views = config['general_plotting']['views'],
        organ_names = config['general_plotting']['organ_names']
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/analysis/dl_meta_top_genes.py"


rule dl_enrichment:
    input:
        organ_spec_dl = 'results/dl_meta/organ_spec_dl.pckl'
    output:
        enrichment_on_dl_path = 'plots/analysis/interorgan_comparison/dl_meta/enrichment_on_dl.pdf',
        enrichment_results_pval = 'results/dl_meta/enrichemnt_pval.pckl',
        enrichment_results_coef = 'results/dl_meta/enrichemnt_coef.pckl'
    params:
        organ_colors = config['general_plotting']['organ_colors'],
        organ_names = config['general_plotting']['organ_names'],
        views = config['general_plotting']['views'],
        enrichment = config['enrichment']
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/analysis/dl_meta_enrichment.py"
