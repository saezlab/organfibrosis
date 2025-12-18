import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages

# Get input and output paths from snakemake
input_file = snakemake.input[0]
meta_file = snakemake.input[1]

cmap_cells = snakemake.params[0]
cmap_conditions = snakemake.params[1]
technology = snakemake.params[2]
final_meta_columns = snakemake.params[3]


output = snakemake.output.adata
output_pdf = snakemake.output.umap
qc_output = snakemake.output.qc_df
ctype_count = snakemake.output.ctype_count
ctype_count_patient = snakemake.output.ctype_count_patient

# Read metadata
metadata = pd.read_csv(
    meta_file,
    index_col=0,
    dtype={'race': 'category', 'Trichrome % fibrotic': 'category'}
)

# Extract study name from output file name
study = output.split('/')[-1][:-5]
print(study)

# Read AnnData object
adata = sc.read_h5ad(input_file)

# Initialize and standardize AnnData obs columns
adata.obs['annotation_MOFA'] = None
adata.obs['cond_test'] = 'control'
adata.obs['tech'] = technology
adata.obs = adata.obs.rename(
    columns={
        'sample_id': 'sample',
        'disease_code': 'grouping',
        'cell_type': 'cell_type1'
    }
)

# Determine which meta columns to use for merging
meta_columns = list(set(adata.obs.columns) & set(['age', 'sex']))

# Merge metadata and handle study-specific logic
if study == 'Reichart_2022':
    metadata = metadata[metadata['study'] == 'Reichart2022_DCM']
    # Set gene names as var index instead of ensembl ids
    adata.var = adata.var.reset_index().rename(columns={'index': 'ensembl_gene_id'})
    adata.var.index = adata.var['feature_name'].tolist()
    adata.obs = adata.obs.reset_index().merge(
        metadata,
        left_on=['donor_id', 'sex'],
        right_on=['sample_id', 'sex'],
        how='left'
    ).set_index('index')
    adata.obs['study'] = study
elif study == 'Kuppe_2022':
    # adata object already contains info
    pass
else:
    metadata = metadata[metadata['study'] == study]
    adata.obs = adata.obs.reset_index().merge(
        metadata,
        left_on=['sample', 'grouping'] + meta_columns,
        right_on=['sample_id', 'disease_code'] + meta_columns,
        how='left'
    ).set_index('index')
    adata.obs['region'] = 'LV'

# Add cell state info from authors for Chaffin_2022
if study == 'Chaffin_2022':
    cell_states = pd.read_csv(
        'data/reheatHeart/cell_state_annos/Chaffin2022_DCM.csv',
        index_col=0
    ).rename(columns={'cellstate': 'cell_type2'})
    adata.obs = adata.obs.reset_index().merge(
        cell_states,
        left_on='index',
        right_index=True,
        how='left'
    ).set_index('index')

# Set modality and calculate QC metrics
adata.obs['modality'] = 'sn'
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Rename columns for uniform nomenclature across atlases
adata.obs['nCount_RNA'] = adata.X.sum(axis=1)
adata.obs['nFeature_RNA'] = np.count_nonzero(adata.X.astype('int16').toarray(), axis=1)

# Save QC metrics to CSV
qc_df = adata.obs[['pct_counts_mt', 'nFeature_RNA', 'nCount_RNA']]
qc_df['study'] = study
qc_df.to_csv(qc_output)

# Filter cells and genes for quality
adata_sub = adata[adata.obs['pct_counts_mt'] < 25, :]
sc.pp.filter_cells(adata_sub, min_genes=200)
sc.pp.filter_cells(adata_sub, min_counts=500)
sc.pp.filter_cells(adata_sub, max_counts=40000)
sc.pp.filter_genes(adata_sub, min_cells=3)

# Ensure UMAP coordinates are present
if ('X_umap' not in adata_sub.obsm_keys()) and ('UMAP' in adata_sub.obsm_keys()):
    adata_sub.obsm['X_umap'] = adata_sub.obsm['UMAP'].to_numpy()
elif ('X_umap' not in adata_sub.obsm_keys()) and ('UMAP_HARMONY' in adata_sub.obsm_keys()):
    adata_sub.obsm['X_umap'] = adata_sub.obsm['UMAP_HARMONY']

# Unify annotations to common annotation across tissues
if study != 'Kuppe_2022':
    adata_sub.obs.loc[
        adata_sub.obs['cell_type1'].str.contains('vSMCs|PC', case=False),
        'annotation_MOFA'
    ] = 'pericytesSMCs'
    adata_sub.obs.loc[
        adata_sub.obs['cell_type1'].str.contains('Fib', case=False),
        'annotation_MOFA'
    ] = 'fibroblast'
    adata_sub.obs.loc[
        adata_sub.obs['cell_type1'].isin(['Myeloid', 'Mast']),
        'annotation_MOFA'
    ] = 'myeloid'
    adata_sub.obs.loc[
        adata_sub.obs['cell_type1'].str.contains('Lymphoid', case=False),
        'annotation_MOFA'
    ] = 'lymphoid'
    adata_sub.obs.loc[
        adata_sub.obs['cell_type1'].str.contains('Endo', case=False),
        'annotation_MOFA'
    ] = 'endothelial'
    adata_sub.obs['cond_test'] = 'control'
    adata_sub.obs.loc[
        adata_sub.obs['grouping'] != 'NF',
        'cond_test'
    ] = 'fibrosis'
    adata_sub.obs['grouping'] = adata_sub.obs['grouping'].replace({'NF': 'control'})
    adata_sub.obs['sex'] = adata_sub.obs['sex'].replace({'male': 'Male', 'female': 'Female'})
    adata_sub.layers['counts'] = adata_sub.X.copy()
    adata_sub.obs['study'] = study
else:
    # Remove specific regions for Kuppe_2022
    adata_sub = adata_sub[adata_sub.obs['major_labl'] != 'BZ']
    adata_sub = adata_sub[adata_sub.obs['major_labl'] != 'RZ']
    adata_sub = adata_sub[adata_sub.obs['major_labl'] != 'IZ']
    adata_sub.obs.loc[
        adata_sub.obs['cell_type1'].isin(['vSMCs', 'PC']),
        'annotation_MOFA'
    ] = 'pericytesSMCs'
    adata_sub.obs.loc[
        adata_sub.obs['cell_type1'].str.contains('Fib', case=False),
        'annotation_MOFA'
    ] = 'fibroblast'
    adata_sub.obs.loc[
        adata_sub.obs['cell_type1'].isin(['Myeloid', 'Mast']),
        'annotation_MOFA'
    ] = 'myeloid'
    adata_sub.obs.loc[
        adata_sub.obs['cell_type1'].str.contains('Lymphoid', case=False),
        'annotation_MOFA'
    ] = 'lymphoid'
    adata_sub.obs.loc[
        adata_sub.obs['cell_type1'].str.contains('Endo', case=False),
        'annotation_MOFA'
    ] = 'endothelial'
    adata_sub.obs['cond_test'] = 'control'
    adata_sub.obs.loc[
        adata_sub.obs['major_labl'] != 'CTRL',
        'cond_test'
    ] = 'fibrosis'
    adata_sub.obs['grouping'] = adata_sub.obs['major_labl'].replace(
        {'CTRL': 'control', 'FZ': 'MI-FZ'}
    )
    adata_sub.layers['counts'] = adata_sub.X.copy()
    adata_sub.obs['study'] = study
    adata_sub.obs = adata_sub.obs.rename(columns={'patient_region_id': 'sample'})
    del adata_sub.obsm['PCA']
    del adata_sub.obsm['HARMONY']

# Set ctype and mesenchymal annotation
adata_sub.obs['ctype'] = adata_sub.obs['annotation_MOFA']
adata_sub.obs.loc[
    adata_sub.obs['annotation_MOFA'].isin(['pericytesSMCs', 'fibroblast']),
    'annotation_MOFA'
] = 'mesenchymal'

# Keep only cells with annotation
adata_new = adata_sub[adata_sub.obs['annotation_MOFA'].notna()].copy()

# Summarize cell type counts
counts = {'endothelial': [], 'myeloid': [], 'lymphoid': [], 'mesenchymal': [], 'epithelial': []}
ctypes = ['endothelial', 'myeloid', 'lymphoid', 'epithelial', 'mesenchymal']
for i in ctypes:
    count = adata_sub.obs[adata_sub.obs['annotation_MOFA'] == i].shape[0]
    counts[i] = count
counts['other'] = adata_sub.obs[~adata_sub.obs['annotation_MOFA'].isin(ctypes)].shape[0]
counts['pericytesSMCs'] = adata_sub.obs[adata_sub.obs['ctype'] == 'pericytesSMCs'].shape[0]
counts['fibroblast'] = adata_sub.obs[adata_sub.obs['ctype'] == 'fibroblast'].shape[0]
# Save cell type counts
counts_df = pd.DataFrame({study: counts})
counts_df.to_csv(ctype_count)

# Also get counts per patient
if pd.api.types.is_categorical_dtype(adata_sub.obs['annotation_MOFA']):
    adata_sub.obs['annotation_MOFA'] = adata_sub.obs['annotation_MOFA'].cat.add_categories(['other'])
obs_data = adata_sub.obs.copy()
obs_data['annotation_MOFA'].fillna('other', inplace=True)
counts = obs_data.groupby(['sample', 'annotation_MOFA']).size().reset_index()
pivoted_df = counts.pivot(index='sample', columns='annotation_MOFA', values=0)
counts_withcond = pivoted_df.merge(
    obs_data[['sample', 'grouping', 'cond_test']],
    on='sample', how='left'
).drop_duplicates()
counts_withcond.to_csv(ctype_count_patient)

# Plot UMAPs and save to PDF
with PdfPages(output_pdf) as output_pdf:
    for i in [adata_sub, adata_new]:
        fig, axs = plt.subplots(1, 3, tight_layout=True, figsize=(15, 4))
        sc.pl.umap(i, color=['annotation_MOFA'], ax=axs[0], show=False, palette=cmap_cells)
        sc.pl.umap(i, color=['cond_test'], ax=axs[1], show=False, palette=cmap_conditions)
        sc.pl.umap(i, color=['sample'], ax=axs[2], show=False)
        output_pdf.savefig(fig)

# Print summary info
print(adata_sub.obs.head())
print(adata_sub.obs.dtypes)

# Keep only relevant metadata columns
cols_to_keep = list(set(adata_sub.obs.columns) & set(final_meta_columns))
adata_sub.obs = adata_sub.obs[cols_to_keep]

# Save processed AnnData object
adata_sub.write(output)
