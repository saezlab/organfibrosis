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

rule filter_co:
    input:
        data_path = "results/dl_meta/cross_organ_dl.pckl"
    output:
        out_path = "results/resource_website/co_data.pckl"
    conda:
        "../envs/scanpy.yaml"
    resources:
        runtime=60
    script:
        "../scripts/dataformat/filter_co_data.py"

rule filter_difib_co:
    input: 
        data_path = "results/myofib/dl_meta/cross_organ_dl.pckl"
    output:
        out_path = "results/resource_website/difib_co_data.pckl"
    conda:
        "../envs/scanpy.yaml"
    resources:
        runtime=60
    script:
        "../scripts/dataformat/filter_difib_co_data.py"

rule filter_spatial:
    input: 
        data_path = "results/spatial/all_organs.csv"
    output:
        out_path = "results/resource_website/spatial_data.pckl"
    params:
        organ_mapping = config["general_plotting"]["organ_names"]
    conda:
        "../envs/scanpy.yaml"
    resources:
        runtime=60
    script:
        "../scripts/dataformat/filter_spatial_data.py"

rule merge_ccc:
    input:
        data_path = "results/ccc/on_dl_allorgans.pckl"
    output:
        out_path = "results/resource_website/ccc_interactions.pckl"
    params:
        organs = config["meta_organs"],
        celltypes = config["general_plotting"]["views"],
        organ_mapping = config["general_plotting"]["organ_names"]
    conda:
        "../envs/scanpy.yaml"
    resources:
        runtime=60
    script:
        "../scripts/dataformat/merge_ccc_data.py"

rule merge_difib_deg:
    input:
        deg_file_list = expand('results/myofib/deg_results/{placeholder[0]}/{placeholder[1]}.csv', placeholder = expanded_list),
    output:
        out_path = "results/resource_website/difib_deg_data.parquet"
    params:
        real_names = config["general_plotting"]["organ_names"]
    conda:
        "../envs/parquet.yaml"
    resources:
        runtime=60
    script:
        "../scripts/dataformat/merge_difib_deg_data.py"

rule merge_deg:
    input:
        deg_file_list = expand('results/deg/{placeholder[0]}/{placeholder[1]}_all_deg_tval.csv', placeholder = expanded_list),
    output:
        out_path = "results/resource_website/deg_data.parquet"
    params:
        real_names = config["general_plotting"]["organ_names"],
        celltypes = config["general_plotting"]["views"]
    conda:
        "../envs/parquet.yaml"
    resources:
        runtime=60
    script:
        "../scripts/dataformat/merge_deg_files.py"


rule merge_enrich:
    input:
        pval_path = "results/dl_meta/enrichemnt_pval.pckl",
        coef_path = "results/dl_meta/enrichemnt_coef.pckl"
    params:
        real_names = config["general_plotting"]["organ_names"],
    output:
        out_path = "results/resource_website/enrich_collectri.pckl"
    conda:
        "../envs/scanpy.yaml"
    resources:
        runtime=60
    script:
        "../scripts/dataformat/merge_enrich_data.py"

rule merge_difib_enrich:
    input:
        pval_path = "results/myofib/dl_meta/enrichemnt_pval.pckl",
        coef_path = "results/myofib/dl_meta/enrichemnt_coef.pckl"
    params:
        real_names = config["general_plotting"]["organ_names"],
    output:
        out_path = "results/resource_website/enrich_difib_collectri.pckl"
    conda:
        "../envs/scanpy.yaml"
    resources:
        runtime=60
    script:
        "../scripts/dataformat/merge_difib_enrich_data.py"

rule merge_colocalization:
    input:
        data_path = expand('results/spatial/cosine_sim/{placeholder}.csv', placeholder = human_studies),
    output:
        out_path = "results/resource_website/ccc_colocalization.pckl"
    params:
        real_names = config["general_plotting"]["organ_names"],
    conda:
        "../envs/parquet.yaml"
    resources:
        runtime=60
    script:
        "../scripts/dataformat/merge_colocalization_data.py"