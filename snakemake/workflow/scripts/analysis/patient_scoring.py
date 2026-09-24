# Author: Leonie Küchenhoff

# In this script, different gene sets are enriched in each sample
# This is supposed to serve as an estimation for how fibrotic/immune active 
# each of the samples is

import scanpy as sc
import decoupler as dc
import pandas as pd
import numpy as np

# Get input and output paths from snakemake
adata_path = snakemake.input[0]
filter_by_expr = snakemake.params[0].get('filter_by_expr')
filter_by_prop = snakemake.params[0].get('filter_by_prop')
geneset_path = snakemake.input.genesets
out_path_est = snakemake.output[0]
out_path_pval = snakemake.output[1]
pbulk_path = snakemake.output[2]

# Read in anndata object
adata = sc.read(adata_path)
# exclude non-anotated cell types
adata = adata[~adata.obs['annotation_MOFA'].isna()].copy()
# Extract study name from metadata
study = adata.obs.iloc[0, :].loc['study']
# Extract number of cell types per sample
cell_types = adata.obs.groupby('sample')['annotation_MOFA'].nunique()
valid_samples = list(cell_types[cell_types > 2].index)

# OCEAN data was SoupX corrected, following step cannot handle floats
if study == 'McCown_2025_sn':
    adata.layers['counts'] = np.round(adata.X)

# make pseudobulk of all cells
full_pb = dc.get_pseudobulk(adata, sample_col='sample', 
                            groups_col = None, layer='counts', 
                            min_prop=0, min_smpls=0)
# filter for valid samples
full_pb= full_pb[full_pb.obs["sample"].isin(valid_samples)]

# Save pseudobulk data to CSV
pb_dat = pd.DataFrame(full_pb.X)
pb_dat.columns = full_pb.var.index.values
pb_dat.index = full_pb.obs.index.values
pb_dat.to_csv(pbulk_path)

# Free memory
del(adata)
# Normalize pseudobulk data
sc.pp.normalize_total(full_pb, target_sum=1e4)
sc.pp.log1p(full_pb)

# Read in curated gene sets
full_genesets = pd.read_csv(geneset_path, index_col=0)
full_genesets = full_genesets.fillna(1)

# Run enrichment analysis using decoupler ULM
# This estimates gene set activity per sample
acts = None
try:
    dc.run_ulm(
        full_pb,
        net=full_genesets,
        source='geneset',
        target='genesymbol',
        use_raw=False,
        weight='weight',
        verbose=True
    )
    acts = dc.get_acts(full_pb, obsm_key='ulm_estimate')
except Exception as e:
    print(f"Error running decoupler ULM: {e}")
    raise

# Extract results from enrichment and combine with metadata
estimate = acts.obsm['ulm_estimate']
pvals = acts.obsm['ulm_pvals']

# Define common columns to keep from metadata
columns = list(set(acts.obs.columns) & set([
    'cond_test', 'Fibrosis', 'eGFR', 'study', 'sex', 'batch', 'tech', 'age',
    'fibrosis score (interstitial fibrosis) in %', 'ischemia time in sec',
    'LVEF', 'BMI', 'Trichrome % fibrotic', 'modality', 'grouping',
    'region', 'sample', 'IFTA', 'clin_global_glomerulosclerosis'
]))

# Add metadata to results
estimate.loc[:, columns] = acts.obs.loc[:, columns]
pvals['study'] = study

# Save results to CSV
estimate.to_csv(out_path_est)
pvals.to_csv(out_path_pval)





