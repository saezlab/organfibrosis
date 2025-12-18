#!env bash
echo $CONDA_PREFIX

mkdir -p logs

if [ -f logs/scDist.post-deploy.log ]; then
    rm logs/scDist.post-deploy.log
fi

$CONDA_PREFIX/bin/Rscript scripts/snakemake/workflow/envs/scDist.R >> logs/scDist.post-deploy.log 2>&1
