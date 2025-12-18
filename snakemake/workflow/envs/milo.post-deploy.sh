#!env bash
echo $CONDA_PREFIX

mkdir -p logs

if [ -f logs/milo.post-deploy.log ]; then
    rm logs/milo.post-deploy.log
fi

$CONDA_PREFIX/bin/Rscript scripts/snakemake/workflow/envs/milo.R >> logs/milo.post-deploy.log 2>&1
