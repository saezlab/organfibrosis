wildcard_constraints:
    organ='[a-zA-Z]+'



rule convert_to_adata:
    input:
       'data/kidney/kpmp_v15_sc.h5Seurat',
       '../KiPfi_KPMPdata/data/KPMP/adjudication.csv'
    output:
        'data/kidney/kpmp_v15_sc.h5ad'
    params:
        slot = 'counts',
        reductions = 'umap'
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/preprocessing/seurat_to_h5ad.py"


rule prep_adata:
    input:
       'data/lung/{study}.h5ad'
    output:
       'results/preprocessing/lung/{study}.h5ad',
       'plots/preprocessing/lung/{study}_umap.pdf'
    params:
        config['general_plotting'].get('celltype_colors'),
        config['general_plotting'].get('condition_colors') 
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/preprocessing/prep_adata.py"



rule Reichart_extra_step:
    input:
       '../reheat2/raw/Reichart2022_DCM/scell_all.h5ad'
    output:
       'data/reheatHeart/Reichart_2022_pp.h5ad'
    conda:
        "../envs/scanpy.yaml"
    resources:
         mem_mb=180000
    script:
        "../scripts/preprocessing/RH_extra_pp.py"


rule Muto_extra_step:
    input:
        adata = ['data/kidney/Muto_2022/GSM5627690_cont1_filtered_feature_bc_matrix.h5', 
                'data/kidney/Muto_2022/GSM5627691_cont2_filtered_feature_bc_matrix.h5', 
                'data/kidney/Muto_2022/GSM5627692_cont3_filtered_feature_bc_matrix.h5', 
                'data/kidney/Muto_2022/GSM5627693_cont4_filtered_feature_bc_matrix.h5', 
                'data/kidney/Muto_2022/GSM5627694_cont5_filtered_feature_bc_matrix.h5',
                'data/kidney/Muto_2022/GSM5627695_PKD1_filtered_feature_bc_matrix.h5', 
                'data/kidney/Muto_2022/GSM5627696_PKD2_filtered_feature_bc_matrix.h5', 
                'data/kidney/Muto_2022/GSM5627697_PKD3_filtered_feature_bc_matrix.h5',
                'data/kidney/Muto_2022/GSM5627698_PKD4_filtered_feature_bc_matrix.h5',
                'data/kidney/Muto_2022/GSM5627699_PKD5_filtered_feature_bc_matrix.h5',
                'data/kidney/Muto_2022/GSM5627700_PKD6_filtered_feature_bc_matrix.h5',
                'data/kidney/Muto_2022/GSM5627701_PKD7_filtered_feature_bc_matrix.h5',
                'data/kidney/Muto_2022/GSM5627702_PKD8_filtered_feature_bc_matrix.h5'],
        meta = 'data/kidney/Muto_2022/GSE185948_metadata_RNA.csv'
    output:
       'data/kidney/Muto_2022_pp.h5ad'
    conda:
        "../envs/scanpy.yaml"
    resources:
         mem_mb=180000
    script:
        "../scripts/preprocessing/Muto_extra_pp.py"

rule Li_extra_step:
    input:
       meta = 'data/kidney/Li_2024/GSE234788_GEO_RNA_446267cells_meta_clinical_data.csv',
       cell_meta = 'data/kidney/Li_2024/GSE234788_GEO_RNA_446267cells_meta.csv'
    output:
       'data/kidney/Li_2024_pp.h5ad'
    conda:
        "../envs/scanpy.yaml"
    resources:
         mem_mb=180000
    script:
        "../scripts/preprocessing/Li_extra_pp.py"


rule prep_adata_heartreheat:
    input:
        lambda w: config['datasets']['reheatHeart']['paths'][w.study],
        config['datasets']['reheatHeart']['metadata'],
    output:
       adata = 'results/preprocessing/reheatHeart/{study}.h5ad',
       umap = 'plots/preprocessing/reheatHeart/{study}_umap.pdf',
       qc_df = 'results/preprocessing/qc/reheatHeart/{study}.csv',
       ctype_count = 'results/preprocessing/qc/reheatHeart/{study}_ctypecount.csv',
       ctype_count_patient = 'results/preprocessing/qc/reheatHeart/{study}_ctypecount_patient.csv',
    params:
        config['general_plotting'].get('celltype_colors'),
        config['general_plotting'].get('condition_colors'),
        lambda w: config['datasets']['reheatHeart']['technology'][w.study],
        config['preprocessing']['metacolumns']
    resources:
         mem_mb=170000
    conda:
        "../envs/scanpy.yaml"
    script:
        "../scripts/preprocessing/prep_adata_reheat.py"


rule prep_adata_lungHCA:
    input:
       'data/HCAlung/scRNA_unprocessed.h5ad'
    output:
        adata = expand('results/preprocessing/HCAlung/{study}.h5ad', study = config['datasets']['HCAlung']['studies_ownmodel']),
        umap = expand('plots/preprocessing/HCAlung/{study}_umap.pdf', study = config['datasets']['HCAlung']['studies_ownmodel']),
        qc_df = expand('results/preprocessing/qc/HCAlung/{study}.csv', study = config['datasets']['HCAlung']['studies_ownmodel']),
        ctype_count = expand('results/preprocessing/qc/HCAlung/{study}_ctypecount.csv', study = config['datasets']['HCAlung']['studies_ownmodel']),
        ctype_count_patients  = expand('results/preprocessing/qc/HCAlung/{study}_ctypecount_patient.csv', study = config['datasets']['HCAlung']['studies_ownmodel'])
    params:
        config['general_plotting'].get('celltype_colors'),
        config['general_plotting'].get('condition_colors'),
        config['preprocessing']['metacolumns']
    conda:
        "../envs/scanpy.yaml"
    resources:
         mem_mb=270000
    script:
        "../scripts/preprocessing/prep_adata_hca.py"
        
        
rule prep_adata_kidney:
    input:
        lambda w: config['datasets']['kidney']['paths'][w.study],
    output:
        adata = 'results/preprocessing/kidney/{study}.h5ad',
        umap = 'plots/preprocessing/kidney/{study}_umap.pdf',
        qc_df = 'results/preprocessing/qc/kidney/{study}.csv',
        ctype_count = 'results/preprocessing/qc/kidney/{study}_ctypecount.csv',
        ctype_count_patient = 'results/preprocessing/qc/kidney/{study}_ctypecount_patient.csv',
    params:
        config['general_plotting'].get('celltype_colors'),
        config['general_plotting'].get('condition_colors'),
        lambda w: config['datasets']['kidney']['metadata'][w.study],
        config['preprocessing']['metacolumns']
    conda:
        "../envs/scanpy.yaml"
    resources:
         mem_mb=300000
    script:
        "../scripts/preprocessing/prep_adata_kidney.py"



rule prep_adata_liver:
    input:
        lambda w: config['datasets']['liver']['paths'][w.study],
    output:
        adata = 'results/preprocessing/liver/{study}.h5ad',
        umap = 'plots/preprocessing/liver/{study}_umap.pdf',
        qc_df = 'results/preprocessing/qc/liver/{study}.csv',
        ctype_count = 'results/preprocessing/qc/liver/{study}_ctypecount.csv',
        ctype_count_patient = 'results/preprocessing/qc/liver/{study}_ctypecount_patient.csv',
    params:
        config['general_plotting'].get('celltype_colors'),
        config['general_plotting'].get('condition_colors'),
        lambda w: config['datasets']['liver']['metadata'][w.study],
        lambda w: config['datasets']['liver']['anno_file'][w.study],
        config['preprocessing']['metacolumns']
    conda:
        "../envs/scanpy.yaml"
    resources:
         mem_mb=200000
    script:
        "../scripts/preprocessing/prep_adata_liver.py"


rule MOFA_prep:
    input:
        'results/preprocessing/{organ}/{study}.h5ad'
    output:
        'results/mofa_input/{organ}/metadata/{study}.csv',
        'results/mofa_input/{organ}/pbulk/{study}.csv',
        'results/mofa_input/{organ}/coldata/{study}.csv',
    params:
        lambda w: config['datasets'][w.organ].get('days_post_treatment'),
        config['preprocessing']['metacolumns_unique']
    conda:
        "../envs/scanpy.yaml"
    resources:
         mem_mb=100000
    script:
        "../scripts/preprocessing/MOFA_prep.py"


rule get_markers:
    input:
        pbulk = 'results/mofa_input/{organ}/pbulk/{study}.csv',
        coldata = 'results/mofa_input/{organ}/coldata/{study}.csv'  
    output:
        marker = 'results/mofa_input/{organ}/mrkrs/{study}.csv',
    conda:
        "../envs/mofa.yaml"
    script:
        "../scripts/preprocessing/get_markers.R"



