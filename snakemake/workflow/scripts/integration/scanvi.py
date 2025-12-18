import os
import scvi
import scanpy as sc
import torch
from rich import print
import argparse

assert torch.cuda.is_available(), "No GPU found!"

# define settings for scvi
scvi.settings.seed = 7
torch.set_float32_matmul_precision("high")


parser = argparse.ArgumentParser()
# Define command line arguments
parser.add_argument('-i' , "--input", help="Path to adata", required=True)
parser.add_argument('-oscan', '--outputscan', help="Path to output scanvi_model", required=True)
parser.add_argument('-oscvi', '--outputscvi', help="Path to output scvi model", required=True)
parser.add_argument('-adata', '--adata', help="Path to output anndata with embeddings and umap", required=True)
parser.add_argument('-l', '--label_col', help="Column name of cell type labels in adata.obs", required=True)
parser.add_argument('-b', '--batch_col', help="Column name of bacthlabels in adata.obs", required=True)

# Parse the command line arguments
args = parser.parse_args()

input_path = args.input
scanvi_model_dir = args.outputscan
scvi_model_dir = args.outputscvi
adata_dir = args.adata
label_col = args.label_col
batch_col = args.batch_col

#  Load data
print('INFO: Loading data')
adata = sc.read_h5ad(input_path)
adata.layers['counts'] = adata.X.copy()
print(adata)

# Preprocess data
sc.pp.normalize_total(adata, target_sum=1e5)
sc.pp.log1p(adata)

adata.raw = adata

# select highly variable genes
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0075,
    max_mean=4,
    min_disp=0.1,
    batch_key=batch_col
)
keep = adata.var['highly_variable']
adata = adata[:, keep].copy()

# Setup scVI model
scvi.model.SCVI.setup_anndata(adata, 
                              layer="counts", 
                              batch_key=batch_col
                              )

model = scvi.model.SCVI(adata)

# Train model
model.train()

SCVI_LATENT_KEY = "X_scVI"
model.adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation(model.adata)

# Save model
print('INFO: Saving model')
model.save(os.path.dirname(scvi_model_dir), overwrite=True, save_anndata=False)

# make scANVI model from scVI model
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    labels_key=label_col,
    unlabeled_category="Unknown",
)

print('INFO: Training scANVI model')
scanvi_model.train()

SCANVI_LATENT_KEY = "X_scANVI"
scanvi_model.adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(scanvi_model.adata)

# Save model
print('INFO: Saving model')
scanvi_model.save(os.path.dirname(scanvi_model_dir), overwrite = True, save_anndata = True)


scanvi_model.adata.obsm[SCVI_LATENT_KEY] = model.adata.obsm[SCVI_LATENT_KEY].copy()

# Make embeddings and clustering

for model in ["scVI", "scANVI"]:
    sc.pp.neighbors(adata, use_rep="X_" + model)
    sc.tl.leiden(adata, key_added='leiden_' + model)
    adata.obsm['X_' + model + '_MDE'] = scvi.model.utils.mde(adata.obsm['X_' + model])

#  Save anndata
adata.write_h5ad(adata_dir)