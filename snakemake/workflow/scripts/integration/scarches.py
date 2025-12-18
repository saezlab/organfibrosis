import os
import scvi
import scanpy as sc
import torch
from rich import print
import argparse
import pandas as pd

assert torch.cuda.is_available(), "No GPU found!"

# define settings for scvi
scvi.settings.seed = 7
torch.set_float32_matmul_precision("high")

    

parser = argparse.ArgumentParser()
# Define command line arguments
parser.add_argument('-m' , "--model", help="Path to reference model", required=True)
parser.add_argument('-q', '--query', help="Path to query dataset (.h5ad)", required=True)
parser.add_argument('-o', '--output', help="Path to output directory", required=True)
parser.add_argument('-l', '--label_col', help="Name of column containing labels", required=True)
parser.add_argument('-ol', '--old_label', help="Name of column containing old ctype labels", required=True)
parser.add_argument('-w', '--wd', help="Weight decay for mapping", required=True)

# Parse the command line arguments
args = parser.parse_args()

model_path = args.model
query_adata_path = args.query

label_col = args.label_col
old_label = args.old_label
weight_decay = float(args.wd)
print('INFO: Weight decay set to {}'.format(weight_decay))
print('INFO: weight decay type: {}'.format(type(weight_decay)))

model_dir = args.output

#reference model
model_path = os.path.dirname(model_path)

#where to store the query model
model_dir = os.path.dirname(model_dir)

# check if paths exist
for path in [model_path, query_adata_path]:
    assert os.path.exists(path), "Path does not exist: {}".format(path)

# Load data


# load ref model
ref_model = scvi.model.SCANVI.load(model_path)

#load query adata
query_adata = sc.read_h5ad(query_adata_path)
query_adata.layers['counts'] = query_adata.X.copy()
query_adata.obs['grouping'] = 'Unknown'

# reformat label column and save old annotations
query_adata.obs['label_before'] = query_adata.obs[old_label].copy()
query_adata.obs[label_col] = 'Unknown'

# prepare the query for annotation
# will also output some info about common vars
scvi.model.SCANVI.prepare_query_anndata(query_adata, ref_model)

query = scvi.model.SCANVI.load_query_data(
    query_adata,
    ref_model,
)

#  Train model on query

query.train(max_epochs=200,
            plan_kwargs={"weight_decay": weight_decay}
            )



SCANVI_PREDICTIONS_KEY = label_col + '_pred'

# get latent representation of query
query.adata.obsm['X_scANVI'] = query.get_latent_representation()

#  Hard prediction of cell types

# predict labels for query cells
query.adata.obs[SCANVI_PREDICTIONS_KEY] = query.predict()

# get soft predictions for each cell type for each cell
# query.adata.obsm['scArches_soft_' + label_col + '_pred'] = query.predict(soft = True)

#  Save results
print('INFO: Saving model')
query.save(model_dir, overwrite = True, save_anndata = True)