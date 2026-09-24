import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from pathlib import Path

plt.rcParams.update({'font.size': 16})

config = snakemake.config
real_names = config['general_plotting']['organ_names']
ctype_colors = config['general_plotting']['celltype_colors']

# Read inputs/outputs provided by the rule
input_paths = list(snakemake.input.pseudobulks) if hasattr(snakemake.input, 'pseudobulks') else list(snakemake.input)
output_paths = list(snakemake.output.pdfs) if hasattr(snakemake.output, 'pdfs') else list(snakemake.output)

# Build organ -> input file list mapping using the real name (directory) mapping
organ_to_inputs = {}
for out in output_paths:
    organ = Path(out).stem.replace('_pbulk', '')
    organ_to_inputs[organ] = []

for p in input_paths:
    matched = False
    for organ, real in real_names.items():
        if str(real) in p:
            organ_to_inputs.setdefault(organ, []).append(p)
            matched = True
            break
    if not matched:
        # fallback: match by organ slug appearing in path
        for organ in organ_to_inputs.keys():
            if organ in p:
                organ_to_inputs[organ].append(p)
                matched = True
                break

for organ, paths in organ_to_inputs.items():
    if not paths:
        continue

    studies = [Path(pp).stem for pp in paths]
    fig, ax = plt.subplots(len(studies), 2, figsize=(7, len(studies) * 3.5), tight_layout=True, sharey=True)
    if len(studies) == 1:
        ax = np.array([ax])

    for study_count, pp in enumerate(paths):
        study = Path(pp).stem
        print(f"Reading {pp}")
        adata = sc.read(pp)

        sns.scatterplot(y=np.log10(adata.obs['psbulk_counts']), x=np.log10(adata.obs['psbulk_n_cells']),
                        ax=ax[study_count, 0], hue=adata.obs['cond_test'],
                        palette={'fibrosis': 'darkslategrey', 'control': 'lightsalmon'})

        sns.scatterplot(y=np.log10(adata.obs['psbulk_counts']), x=np.log10(adata.obs['psbulk_n_cells']),
                        ax=ax[study_count, 1], hue=adata.obs['annotation_MOFA'],
                        palette=ctype_colors)

        ax[study_count, 0].get_legend().remove()
        ax[study_count, 1].get_legend().remove()

        ax[study_count, 0].set_ylim(0.5, 9)
        ax[study_count, 0].set_xlim(0.5, 5)
        ax[study_count, 1].set_xlim(0.5, 5)

        ax[study_count, 0].set_ylabel('log10(counts)')
        ax[study_count, 0].set_xlabel('log10(n cells)')
        ax[study_count, 1].set_xlabel('log10(n cells)')

        ax[study_count, 0].grid(True)
        ax[study_count, 1].grid(True)

        ax[study_count, 0].set_title(f'                                   {study}')

    fig.suptitle(real_names.get(organ, organ))

    # Find intended output path for this organ from provided outputs
    outpath = None
    for op in output_paths:
        if Path(op).stem.replace('_pbulk', '') == organ:
            outpath = Path(op)
            break
    if outpath is None:
        outpath = Path('plots/analysis/pbulkstats') / f"{organ}_pbulk.pdf"

    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(outpath), bbox_inches='tight')
    plt.close(fig)
