#!env bash
echo $CONDA_PREFIX

mkdir -p logs

if [ -f logs/mofa.post-deploy.log ]; then
    rm logs/mofa.post-deploy.log
fi

$CONDA_PREFIX/bin/Rscript scripts/snakemake/workflow/envs/mofa.R >> logs/mofa.post-deploy.log 2>&1
