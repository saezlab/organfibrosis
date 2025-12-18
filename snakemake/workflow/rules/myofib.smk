
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


rule pbulk_deg_myofib:
    input:
        adata_path = 'results/preprocessing/fibroblasts/all_fibs.h5ad',
        scvi_path = 'results/integration/fibroblasts/integrated_latent_space.csv'
    output:
        deg_paths = expand('results/myofib/deg_results/{placeholder[0]}/{placeholder[1]}.csv', placeholder = expanded_list),
        anno_path = 'results/myofib/deg_results/annotation.csv'
    resources:
         mem_mb=150000
    params:
        organs = config['meta_organs'],
        filters = config['preprocessing'].get('pseudobulk'),
        resolutions = config['preprocessing']['fibroblasts']['scvi_resolutions'],
        diseasefibs = config['preprocessing']['fibroblasts']['diseasefibs']
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/myofib/pbulk_myofib.py"


rule dl_meta:
    input:
        deg_file = expand('results/myofib/deg_results/{placeholder[0]}/{placeholder[1]}.csv', placeholder = expanded_list)
    output:
        organ_spec = 'results/myofib/dl_meta/organ_spec_dl.pckl',
        cross_organ = 'results/myofib/dl_meta/cross_organ_dl.pckl',
    params:
        organ_colors = config['general_plotting']['organ_colors'],
        organ_names = config['general_plotting']['organ_names'],
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/myofib/dl_meta.py"


rule dl_enrichment:
    input:
        organ_spec_dl = 'results/myofib/dl_meta/organ_spec_dl.pckl',
    output:
        enrichment_on_dl_path = 'plots/myofib/interorgan_comparison/dl_meta/enrichment_on_dl.pdf',
        enrichment_results_pval = 'results/myofib/dl_meta/enrichemnt_pval.pckl',
        enrichment_results_coef = 'results/myofib/dl_meta/enrichemnt_coef.pckl'
    params:
        organ_colors = config['general_plotting']['organ_colors'],
        organ_names = config['general_plotting']['organ_names'],
        views = config['general_plotting']['views'],
        enrichment = config['enrichment']
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/myofib/dl_meta_enrichment.py"



rule analysis_dl:
    input:
        organ_spec_dl = 'results/myofib/dl_meta/organ_spec_dl.pckl',
        cross_organ_dl = 'results/myofib/dl_meta/cross_organ_dl.pckl',
        deg_file = expand('results/myofib/deg_results/{placeholder[0]}/{placeholder[1]}.csv', placeholder = expanded_list),
    output:
        cross_organ_dl = 'plots/myofib/interorgan_comparison/dl_meta/cross_organ_dl.pdf',
        upsetplot_dl_organspec_genes = 'plots/myofib/interorgan_comparison/dl_meta/upsetplot_dl_organspec_genes.pdf',
        gene_dict_up_output = 'results/myofib/dl_meta/organ_spec_dl_up.pckl',
        gene_dict_down_output = 'results/myofib/dl_meta/organ_spec_dl_down.pckl',
        count_organspec_genes = 'plots/myofib/interorgan_comparison/dl_meta/count_organspec_genes.pdf',
        counts_up = 'plots/myofib/interorgan_comparison/dl_meta/count_organspec_up.csv',
        counts_down = 'plots/myofib/interorgan_comparison/dl_meta/count_organspec_down.csv',
        topx = 'plots/myofib/interorgan_comparison/dl_meta/cross_organ_dl_topx.pdf',
    params:
        organ_colors = config['general_plotting']['organ_colors'],
        organ_names = config['general_plotting']['organ_names'],
        views = config['general_plotting']['views']
    conda:
        "../envs/scanpy.yaml"
    resources:
        runtime = 300
    script:
        "../scripts/myofib/dl_meta_analysis_plot.py"




rule ccc_dl:
    input:
        organ_spec_dl = 'results/dl_meta/organ_spec_dl.pckl',
        myofib_organ_spec_dl = 'results/myofib/dl_meta/organ_spec_dl.pckl'
    output:
        ccc_result = 'results/myofib/ccc/on_dl_allorgans.pckl'
    params:
        organs = config['meta_organs'],
        views = config['general_plotting']['views']
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/myofib/ccc_on_dlresults.py"


rule percent_unique:
    input:
        organ_spec_dl_all = "results/dl_meta/organ_spec_dl.pckl",
        organ_spec_dl_myofib = "results/myofib/dl_meta/organ_spec_dl.pckl"
    params:
        organs = config['meta_organs'],
        all_studies = [config['datasets'][organ]['studies_to_analyze'] for organ in organs],
        org_colors = config['general_plotting']['organ_colors'],
        real_names = config['general_plotting']['organ_names'],
        ctypes = config["general_plotting"]["views"]
    conda:
        "../envs/scanpy.yaml"
    output:
        plot_pdf_all_mesench = "plots/myofib/percent_unique_all_mesenchymal.pdf",
        output_pdf_comparison = "plots/myofib/percent_unique_comparison.pdf",
        csv_all = "plots/myofib/percent_unique_all.csv",
        csv_myofib = "plots/myofib/percent_unique_diseasefib.csv",
    script:
        "../scripts/myofib/percent_unique.py"