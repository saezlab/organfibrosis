# Author: Leonie Küchenhoff

"""
This script preprocesses the study by Li et al. to make the format comparable
to other studies. A large part of this script is taken from the authors of
the paper and can be found here:
https://github.com/TheHumphreysLab/SHARE-seq-kidney/blob/main/Figure1/fig1_RNA_qc_clustering.py

This includes the preprocessing of the barcodes that need to be adapted to
merge with metadata. Another preprocessing step follows after this script.
"""

from typing import Dict, List

import anndata
import pandas as pd
import scanpy as sc

# Define snakemake inputs and outputs
meta_path = snakemake.input['meta']
cell_meta_path = snakemake.input['cell_meta']
output = snakemake.output[0]

# Read metadata files
meta = pd.read_csv(meta_path, index_col=0)
cell_meta = pd.read_csv(cell_meta_path, index_col=0)

def append_batch_id(adata: anndata.AnnData, batch_id: str) -> List[str]:
    """Add batch identifier to cell barcodes.
    
    Args:
        adata: AnnData object containing single-cell data
        batch_id: Identifier for the batch (e.g., 'B1')
        
    Returns:
        List of modified cell barcodes with batch information
    """
    return [f"{i},{batch_id}" for i in adata.obs.index]

# Load data from different sequencing batches
adata1 = sc.read('data/kidney/Li_2024/GSM7474780_S2.rna..h5ad')
adata1.obs.index = append_batch_id(adata1, 'B1')

# Load NovaSeq S4 files
s4_files = {
    'B2': "data/kidney/Li_2024/GSM7474782_S4_1.RNA.hg19.gene.bc.matrices.h5",
    'B3': "data/kidney/Li_2024/GSM7474784_S4_2.RNAlane3.hg19.gene.bc.matrices.h5",
    'B4': "data/kidney/Li_2024/GSM7474785_S4_2.RNAlane4.hg19.gene.bc.matrices.h5"
}

adatas = [adata1]
for batch_id, file_path in s4_files.items():
    adata = sc.read_10x_h5(file_path)
    adata.obs.index = append_batch_id(adata, batch_id)
    adatas.append(adata)

# Combine all batches
adata_merge = anndata.concat(adatas, join="outer")

# Process cell barcodes for batch B1
new_index = []
for idx in adata_merge.obs.index:
    if idx[-2:] == 'B1':
        parts = idx.split('.')
        new_idx = (f"{parts[0]}.{parts[1]},{parts[2]}.{parts[3]},"
                  f"{parts[4]}.{parts[5]},{parts[6]}.{parts[7]}")
        new_index.append(new_idx)
    else:
        new_index.append(idx)
adata_merge.obs.index = new_index

def extract_barcode_ids(barcode: str) -> Dict[str, int]:
    """Extract R1, R2, R3, P1, and batch IDs from barcode string.
    
    Args:
        barcode: String containing the cell barcode
        
    Returns:
        Dictionary with R1, R2, R3, P1, and batch component values
    """
    parts = barcode.split(',')
    return {
        'R1': int(parts[0].split('.')[-1]),
        'R2': int(parts[1].split('.')[-1]),
        'R3': int(parts[2].split('.')[-1]),
        'P1': int(parts[3].split('.')[-1]),
        'B': int(parts[4][-1])
    }

# Process all barcodes
barcode_components = [extract_barcode_ids(i) for i in adata_merge.obs.index]
for component in ['R1', 'R2', 'R3', 'P1', 'B']:
    adata_merge.obs[f'{component}_id'] = [bc[component] for bc in barcode_components]

# Flag invalid barcodes (R1, R2, or R3 > 96)
adata_merge.obs['barcode_error'] = [
    1 if any(adata_merge.obs[f'{r}_id'][i] > 96 for r in ['R1', 'R2', 'R3'])
    else 0
    for i in range(len(adata_merge.obs))
]

# Filter out cells with barcode errors
adata = adata_merge[adata_merge.obs['barcode_error'] == 0, :]

# Merge with metadata while preserving original index order
index_before = adata.obs.index
adata.obs = (
    adata.obs.merge(cell_meta, how='left', left_index=True, right_index=True)
    .merge(meta, left_index=True, right_index=True, how='left')
    .loc[index_before, :]
)

# Replace underscores with hyphens in gene names (required for downstream processing)
adata.var_names = adata.var_names.str.replace('_', '-', regex=True)
adata.var_names_make_unique()

# Save processed data
adata.write(output)

