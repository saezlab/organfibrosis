
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
print(expanded_list)


all_human_studies = all_studies[all_studies['organ'].isin(human_studies)]
organs = all_studies['organ'].unique()


rule heart:
    input:
        ccc="results/ccc/on_dl_allorgans.pckl",
        genesets=[
            "data/misc/genesets/naba/NABA_PROTEOGLYCANS.v2023.2.Hs.gmt",
            "data/misc/genesets/naba/NABA_COLLAGENS.v2023.2.Hs.gmt",
        ],
        metadata="../KuppeRamirezLi_MI/visium/visium_meta.csv",
    output:
        cosine="results/spatial/cosine_sim/reheatHeart.csv",
        cosine_extra="results/spatial/cosine_sim/reheatHeart_extrapairs.csv",
        visium_plot="plots/spatial/visium_images/heart_COL1A1_TIMP1.pdf"
    params:
        organs=config['meta_organs'],
        visium_dir="../KuppeRamirezLi_MI/visium/visium_data_py",
        plot_genes=["COL1A1", "TIMP1"]
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/spatial/heart.py"

rule kidney:
    input:
        ccc="results/ccc/on_dl_allorgans.pckl",
        genesets=[
            "data/misc/genesets/naba/NABA_PROTEOGLYCANS.v2023.2.Hs.gmt",
            "data/misc/genesets/naba/NABA_COLLAGENS.v2023.2.Hs.gmt",
        ],
        adata="data/kidney/visium.h5ad",
    output:
        cosine="results/spatial/cosine_sim/kidney.csv",
        cosine_extra="results/spatial/cosine_sim/kidney_extrapairs.csv",
        visium_plot="plots/spatial/visium_images/kidney_COL1A1_TIMP1.pdf"
    params:
        organs=config['meta_organs'],
        sample_var="uid_slide",
        condition_var="condition_full",
        exclude_samples=[
            "V12D05-071_XY01_Ref",
            "V12D05-071_XY02_K21-00202-1",
            "V12D05-072_XY01_K21-198-3-C",
            "V12D05-072_XY02_K21-198-3-CM",
            "V12D05-072_XY03_K21-198-3-M",
            "V12D05-072_XY04_K22-265-5",
            "V12D05-074_XY04_Ref",
            "V12N14-057_XY01_Ref",
            "V12N14-082_XY03_Ref",
            "V12N14-084_XY01_K21-139-1-1",
            "V12N14-084_XY02_K21-139-1-2",
            "V12N16-374_XY01_Ref",
            "V12U06-350_XY03_Ref",
            "V12U21-008_XY03_K21-272-2",
            "V12U21-008_XY03_Ref",
            "V12U21-008_XY04_K21-272-3",
            "V12U21-009_XY01_Ref",
            "V12U21-009_XY04_K22-175-5",
            "V12U21-010_XY01_Ref",
            "V13F08-024_XY02_Ref-flip",
            "V13F08-026_XY04_Ref",
            "V13J17-332_XY02_Ref",
            "V13M27-060_XY04_Ref",
            "V13Y22-308_XY01_23-0156",
            "V13Y22-309_XY04_Ref",
            "V19S25-016_XY01_18-0006",
            "V19S25-019_XY02_M32",
            "V19S25-019_XY03_M61",
            "V19S25-019_XY04_F52",
            "V52L26-067_XY01_K21-198-3",
            "V52L26-067_XY02_K21-273-3",
            "V53M20-051_XY01_K22-284-1",
            "V53M20-103_XY01_K22-400-1",
        ],
        plot_genes=["COL1A1", "TIMP1"],
        plot_samples={"ref": "V13Y22-309_XY03_23-0160", "fib": "V12U06-350_XY01_22-0081"},
    conda:
        "../envs/scanpy.yaml"
    resources:
        mem_mb=200000
    script:
        "../scripts/spatial/kidney.py"

rule liver:
    input:
        ccc="results/ccc/on_dl_allorgans.pckl",
        genesets=[
            "data/misc/genesets/naba/NABA_PROTEOGLYCANS.v2023.2.Hs.gmt",
            "data/misc/genesets/naba/NABA_COLLAGENS.v2023.2.Hs.gmt",
        ],
        visium_files=expand(
            "data/liver/Andrews_2024/{slide}.h5ad",
            slide=[
                "Visium_C73_A1",
                "Visium_C73_B1",
                "Visium_C73_C1",
                "Visium_C73_D1",
                "Visium_PSC011_A1",
                "Visium_PSC011_B1",
                "Visium_PSC011_C1",
                "Visium_PSC011_D1",
            ],
        ),
    output:
        cosine="results/spatial/cosine_sim/liver.csv",
        cosine_extra="results/spatial/cosine_sim/liver_extrapairs.csv",
        visium_plot="plots/spatial/visium_images/liver_COL1A1_TIMP1.pdf"
    params:
        organs=config['meta_organs'],
        visium_dir="data/liver/Andrews_2024",
        slides=[
            "Visium_C73_A1",
            "Visium_C73_B1",
            "Visium_C73_C1",
            "Visium_C73_D1",
            "Visium_PSC011_A1",
            "Visium_PSC011_B1",
            "Visium_PSC011_C1",
            "Visium_PSC011_D1",
        ],
        plot_genes=["COL1A1", "TIMP1"],
        plot_samples={"ref": "Visium_C73_D1", "fib": "Visium_PSC011_A1"},
        condition_map={"primary sclerosing cholangitis": "fibrosis"},
    conda:
        "../envs/scanpy.yaml"
    resources:
        mem_mb=200000
    script:
        "../scripts/spatial/liver.py"

rule lung:
    input:
        ccc="results/ccc/on_dl_allorgans.pckl",
        genesets=[
            "data/misc/genesets/naba/NABA_PROTEOGLYCANS.v2023.2.Hs.gmt",
            "data/misc/genesets/naba/NABA_COLLAGENS.v2023.2.Hs.gmt",
        ],
        adata="data/HCAlung/spatial/adata_vis_human_spatial_paper.h5ad",
    output:
        cosine="results/spatial/cosine_sim/HCAlung.csv",
        cosine_extra="results/spatial/cosine_sim/HCAlung_extrapairs.csv",
        visium_plot="plots/spatial/visium_images/lung_COL1A1_TIMP1.pdf"
    params:
        organs=config['meta_organs'],
        sample_var="sample",
        condition_var="treatment",
        plot_genes=["COL1A1", "TIMP1"],
        plot_samples={
            "ref": "91_A1_RO-727_Healthy_processed_CM",
            "fib": "91_D1_24513-17_IPF_processed_CM",
        },
    conda:
        "../envs/scanpy.yaml"
    resources:
        mem_mb=200000
    script:
        "../scripts/spatial/lung.py"







rule ligand_receptor:
    input:
        genesets=[
            "data/misc/genesets/naba/NABA_PROTEOGLYCANS.v2023.2.Hs.gmt",
            "data/misc/genesets/naba/NABA_COLLAGENS.v2023.2.Hs.gmt"],
        cosine_sim= expand("results/spatial/cosine_sim/{organ}_extrapairs.csv", organ = config['meta_organs']),
        cross_organ_dl="results/myofib/dl_meta/cross_organ_dl.pckl",
        organ_spec_dl="results/myofib/dl_meta/organ_spec_dl_up.pckl",
    output:
        cosine_summary="results/spatial/all_organs.csv",
        scatter_org_plot="plots/spatial/scatter_org_spatial_diseasefibs.pdf",
        scatter_org_plot_csv="plots/spatial/scatter_org_spatial_diseasefibs.csv",
        upset_org_plot="plots/spatial/upset_org_spatial_diseasefibs.pdf",
        upset_org_plot_csv="plots/spatial/upset_org_spatial_diseasefibs.csv",
        scatter_crossorg_plot="plots/spatial/scatter_crossorg_spatial_diseasefibs.pdf",
        scatter_crossorg_plot_csv="plots/spatial/scatter_crossorg_spatial_diseasefibs.csv",
        upset_crossorg_plot="plots/spatial/upset_crossorg_spatial_diseasefibs.pdf",
        upset_crossorg_plot_csv="plots/spatial/upset_crossorg_spatial_diseasefibs.csv",
        top_genes_heatmap="plots/spatial/top_genes_heatmap_spatial_diseasefibs.pdf",
        top_genes_heatmap_csv="plots/spatial/top_genes_heatmap_spatial_diseasefibs.csv",
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/spatial/ligand_receptor.py"


rule ligand_receptor_coex:
    input:
        ccc="results/ccc/on_dl_allorgans.pckl",
        cosine_sim=expand("results/spatial/cosine_sim/{organ}.csv", organ=config['meta_organs'])
    output:
        per_interaction="results/spatial/ligand_receptor_coex/per_interaction.csv",
        merged_matrix="results/spatial/ligand_receptor_coex/merged_matrix.csv",
        top_genes="results/spatial/ligand_receptor_coex/top_genes.txt",
        top_interactors="results/spatial/ligand_receptor_coex/top_interactors.json",
        top_genes_heatmap="plots/spatial/ligand_receptor_coex/top_genes_heatmap.pdf",
        top_genes_heatmap_csv="plots/spatial/ligand_receptor_coex/top_genes_heatmap.csv",
        top_interactors_plot="plots/spatial/ligand_receptor_coex/top_interactors.pdf",
        top_interactors_plot_csv="plots/spatial/ligand_receptor_coex/top_interactors.csv",
        top_genes_boxplot="plots/spatial/ligand_receptor_coex/top_genes_boxplot.pdf",
        top_genes_boxplot_csv="plots/spatial/ligand_receptor_coex/top_genes_boxplot.csv"
    resources:
        runtime=60
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/spatial/ligand_receptor_coex.py"
