# %%
from scipy.sparse import csr_matrix
import h5py
import pandas as pd

import platform
import os
# %%
import scanpy as sc
import numpy as np

# %%
def group_to_array(file, group, as_array=False, as_sparse=False, as_dataframe=False):
    """
    Converts a HDF5 group to a numpy array, sparse matrix or pandas dataframe.
    """
    arrays = []

    if isinstance(file[group], h5py.Group):
        for key, value in file[group].items():
            if isinstance(value, h5py.Dataset):
                if file[group][key].dtype == 'object':
                    value = value.asstr()
                
                if as_array or as_sparse:
                    arrays.append(np.array(value))
                elif as_dataframe:
                    arrays.append(pd.Series(value[:]))

            elif isinstance(value, h5py.Group) and 'levels' in value.keys() and 'levels' in value.keys():
                levels = np.array(value['levels'].asstr())
                if as_array or as_sparse:
                    out = [levels[val - 1] for val in value['values']]
                    arrays.append(np.array(out))
                    
                elif as_dataframe:
                    new_type = pd.CategoricalDtype(levels)
                    codes = [val -1 for val in value['values']]
                    out = pd.Categorical.from_codes(codes, dtype=new_type)
                    arrays.append(pd.Series(out))

        if as_array:
            return np.vstack(arrays).T
        elif as_sparse:
            return csr_matrix(np.vstack(arrays).T)
        elif as_dataframe:
            df = pd.concat(arrays, axis=1)
            df.columns = [key for key in file[group].keys()]
            if '_index' in df.columns:
                df = df.set_index('_index')
                df.index.name = None
            return df
    
    #get embeddings
    elif isinstance(file[group], h5py.Dataset) and 'cell.embeddings' in group:
            return np.array(file[group]).T

def assay_to_sparse_matrix(file, assay = 'RNA', slot = 'data'):
    """
    Converts a HDF5 assay to a csr_matrix.
    """
    if assay in file['assays']:
        nrows = file['meta.data/_index'][:].shape[0]
        features = file['assays/{0}/features'.format(assay)].asstr()[:]
        ncols = features.shape[0]
        slot = 'assays/{0}/{1}'.format(assay, slot)
        
        mat = csr_matrix((file[slot + '/data'][:],
                          file[slot + '/indices'][:],
                          file[slot + '/indptr'][:]),
                         shape=(nrows, ncols), dtype=file[slot + '/data'].dtype)
        return mat, features
    else:
        raise ValueError('Assay {0} not found'.format(assay))

def readh5Seurat_to_anndata(file, assay = 'RNA', slots = ['counts'], reductions=None):
    """ Reads from h5Seurat file and returns anndata object

    """
    with h5py.File(file, "r") as f:

        # check that main assay is in f
        if assay not in list(f['assays'].keys()):
            raise ValueError('Main assay %s not in h5Seurat file' % assay)
        
        # check that obs is in f
        if 'meta.data' not in list(f.keys()):
            raise ValueError('meta.data not in h5Seurat file')

        adata = None
        for slot in slots:
            if slot not in ['counts', 'data']:
                raise ValueError('Slot %s not in h5Seurat file' % slot)
            if adata is None:
                print('Loading main assay ({0}) as {1}'.format(assay, slot))
                mat, feat = assay_to_sparse_matrix(f, assay = assay, slot = 'counts')
                adata = sc.AnnData(X=mat, dtype=mat.dtype)
            else:
                print('Loading layer: {0}_{1}'.format(assay, slot))
                mat, _ = assay_to_sparse_matrix(f, assay = assay, slot = slot)
                adata.layers['{0}_{1}'.format(assay, slot)] = mat
        
        
        adata.var_names = feat

        # add obs
        print('Loading obs')
        adata.obs = group_to_array(f, 'meta.data', as_dataframe=True)
        
        if isinstance(reductions, dict):
            for red in reductions.keys():
                assert red in f['reductions'].keys(), 'Reduction %s not in h5Seurat file' % red
                print('Loading reduction: {0} as {1}'.format(red, 'X_' + reductions[red]))
                embedding = group_to_array(f, 'reductions/{0}/cell.embeddings'.format(red), as_array=True)
                adata.obsm['X_' + reductions[red]] = embedding
        
        if isinstance(reductions, list):

            for red in reductions:
                assert red in f['reductions'].keys(), 'Reduction %s not in h5Seurat file' % red
                print('Loading reduction: %s' % red)
                adata.obsm['X_' + red] = group_to_array(f, 'reductions/{0}/cell.embeddings'.format(red), as_array=True)
            

    return adata

# %%
if 'snakemake' in locals():
    input_file = snakemake.input[0]
    slot = snakemake.params['slot']
    reductions = snakemake.params['reductions']
    metadata = snakemake.input[1]



# %%
if not os.path.exists(input_file):
    raise FileNotFoundError(input_file + ' does not exist')

print('INFO: converting h5Seurat file {0} to anndata'.format(input_file))
# read in h5Seurat file to anndata
adata = readh5Seurat_to_anndata(input_file, reductions=reductions)
metadata_df = pd.read_csv(metadata)




# %%
if 'snakemake' in locals():
    adata.write(snakemake.output[0])