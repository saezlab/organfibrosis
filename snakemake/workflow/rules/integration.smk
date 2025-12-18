import pandas as pd

preprocess_combinations = config['preprocessing']['integration']

study_list = [(organ, study) for organ, vars in preprocess_combinations.items() for study in vars['ground_truth']]
groundtruth_studies = pd.DataFrame(study_list, columns=['organ', 'study'])

# Run scANVI

rule scanvi:
    input:
        adata = lambda w: config['preprocessing']['integration'][w.organ]['path']
    output:
        scanvi_model = 'results/integration/{organ}/scanvi_model/model.pt',
        scvi_model = temp('results/integration/{organ}/scvi_model/model.pt'),
        adata = 'results/integration/{organ}/with_embeddings.h5ad'
    params:
        script = "scripts/snakemake/workflow/scripts/integration/scanvi.py",
        label_col = lambda w: config['preprocessing']['integration'][w.organ]['ctype_ground_truth'],
        batch_col = lambda w: config['preprocessing']['integration'][w.organ]['batch_ground_truth'],
        fake_param = 1
    resources:
        mem_mb = 40000,
        slurm = "gres=gpu:1"
    shell:
        'set +eu '
        ' && . $(conda info --base)/etc/profile.d/conda.sh '
        ' && conda activate /home/hd/hd_hd/hd_lx260/SOFTWARE/miniconda3/envs/scvi '
        ' && $CONDA_PREFIX/bin/python {params.script} '
        ' -i {input.adata} -oscan {output.scanvi_model} -oscvi {output.scvi_model} -adata {output.adata} -l {params.label_col} -b {params.batch_col}'


#Run scArches

rule scArches:
    input:
        ref_model = 'results/integration/{organ}/scanvi_model/model.pt',
        query_data = 'data/{organ}/{study}_pp.h5ad'
    output:
        #adata = 'results/integration/{organ}/post_label_transfer/{study}_scarches_model/adata.h5ad',
        model = 'results/integration/{organ}/post_label_transfer/{study}_scarches_model/model.pt',
    resources:
        mem_mb = 20000,
        slurm = "gres=gpu:1"
    params:
        label_col = lambda w: config['preprocessing']['integration'][w.organ]['ctype_ground_truth'],
        old_ctype = lambda w: config['preprocessing']['integration'][w.organ]['old_ctype'].get(w.study),
        script = "scripts/snakemake/workflow/scripts/integration/scarches.py",
        weight_decay = 0,
        fake_param = 1
    shell:
        'set +eu '
        ' && . $(conda info --base)/etc/profile.d/conda.sh '
        ' && conda activate /home/hd/hd_hd/hd_lx260/SOFTWARE/miniconda3/envs/scvi '
        " && $CONDA_PREFIX/bin/python {params.script} --model {input.ref_model} --query {input.query_data}"
        " -l {params.label_col} -ol {params.old_ctype} --wd {params.weight_decay} -o {output.adata}"



rule scvi_fibroblast_train:
    input:
        adata = 'results/preprocessing/fibroblasts/all_fibs.h5ad'
    output:
        trained = directory('results/integration/fibroblasts/scvi_trained')
    params:
        script = "scripts/snakemake/workflow/scripts/integration/scvi_fibroblasts.py",
        hvg_dir = 'results/preprocessing/fibroblasts/',
        seed = 7,
        max_epochs = 50
    resources:
        mem_mb = 100000,
        slurm = "gres=gpu:1"
    shell:
        'set +eu '
        ' && . $(conda info --base)/etc/profile.d/conda.sh '
        ' && conda activate /home/hd/hd_hd/hd_lx260/SOFTWARE/miniconda3/envs/scvi '
        ' && mkdir -p {output.trained} '
        ' && $CONDA_PREFIX/bin/python {params.script} '
        ' -i {input.adata} -v {params.hvg_dir} --intermediate {output.trained} '
        ' --seed {params.seed} --max-epochs {params.max_epochs}'


rule scvi_fibroblast_cluster:
    input:
        trained = rules.scvi_fibroblast_train.output.trained
    output:
        latent_space = 'results/integration/fibroblasts/integrated_latent_space.csv'
    params:
        script = "scripts/snakemake/workflow/scripts/integration/scvi_fibroblasts_cluster.py",
        plot_dir = "results/integration/fibroblasts/",
        seed = 30
    resources:
        mem_mb = 60000,
        slurm = "gres=gpu:1"
    shell:
        'set +eu '
        ' && . $(conda info --base)/etc/profile.d/conda.sh '
        ' && conda activate /home/hd/hd_hd/hd_lx260/SOFTWARE/miniconda3/envs/scvi '
        ' && mkdir -p {params.plot_dir} '
        ' && $CONDA_PREFIX/bin/python {params.script} '
        ' --intermediate {input.trained} -o {output.latent_space} -p {params.plot_dir} '
        ' --seed {params.seed}'




