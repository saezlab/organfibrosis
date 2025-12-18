# script to extract marker genes from organ-integrated clusters

import scanpy as sc
import pandas as pd
import numpy as np
import decoupler as dc


organ = snakemake.wildcards.get('organ')
fibroblast_path = snakemake.input.adata
scvi_results_path = snakemake.input.scvi_results
marker_output_path_scvi = snakemake.output.marker_output_scvi
full_marker_output_path_scvi = snakemake.output.full_marker_output_scvi
pb_data = snakemake.output.pbulk_scvi
resolutions = snakemake.params['resolutions']



adata = sc.read(fibroblast_path)
org_adata = adata[adata.obs['organ'] == organ]


def read_and_add_harmony_results(adata_to_add, path, suffix, id_study = False, pcs = 50):
    # Attach precomputed harmony/UMAP embeddings and Leiden labels to the supplied AnnData object.
    """
    function to read in PCA, UMAP, and leiden clustering results and 
    add them to adata object with the suffix so that several results can be added
    input:
    adata_to_add - adata object where results will be added
    path - path at wich PCA, UMAP and leiden results are saved
    suffix - suffix is added to names, i.e. UMPA results is saved at X_umap_suffix
    id_study - if id column is just the cell id or a merged cell id + study
    pcs - no. of pcs

    output:
    adata object with adde results
    """
    harm_results = pd.read_csv(path, 
                     index_col = 0, 
                     dtype = {'id':'str','leiden_0.2':str, 
                              'leiden_0.3':str,'leiden_0.4':str,
                              'leiden_0.5':str,'leiden_0.7':str})

    if id_study == False:
        assert(np.all(adata_to_add.obs_names == harm_results['id']))
    elif id_study == True:
        harm_results = harm_results.set_index('Unnamed: 1')
        harm_results = harm_results.loc[(np.array(adata_to_add.obs.index) + '_' + np.array(adata_to_add.obs['study']))]
        print('success')
    else:
        assert(np.all((np.array(adata_to_add.obs.index) + '_' + np.array(adata_to_add.obs['study'])) == harm_results['id']))
        harm_results = harm_results.set_index('id',drop = False)
        print('success')

    adata_to_add.obsm[f'X_harmony_{suffix}'] = np.array(harm_results.iloc[:,:pcs])
    adata_to_add.obsm[f'X_umap_{suffix}'] = np.array(harm_results.iloc[:,pcs + 1:pcs + 3])

    adata_to_add.obs[[f'leiden_0.2_{suffix}',f'leiden_0.3_{suffix}',f'leiden_0.4_{suffix}',
                      f'leiden_0.5_{suffix}',f'leiden_0.7_{suffix}']] = harm_results[['leiden_0.2',
                                                                                        'leiden_0.3',
                                                                                        'leiden_0.4',
                                                                                        'leiden_0.5',
                                                                                        'leiden_0.7']].to_numpy()

    return adata_to_add


org_adata = read_and_add_harmony_results(org_adata, scvi_results_path, f'{organ}_scvi', id_study = True, pcs = 30)

# make new column with cluster resolutions to use per organ
res = resolutions[organ]
org_adata.obs.loc[:, 'cluster_organ_to_use'] = org_adata.obs.loc[:, f'leiden_{str(res)}_{organ}_scvi']



marker_dict = {}
for resolution in [0.2,0.3,0.4,0.5,0.7]:

    sc.tl.rank_genes_groups(org_adata, f'leiden_{str(resolution)}_{organ}_scvi', method='wilcoxon',use_raw=False)

    if (resolution == res):
        # also save all marker info
        groups = org_adata.obs[f'leiden_{str(res)}_{organ}_scvi'].unique()

        full = sc.get.rank_genes_groups_df(org_adata, groups)
        full_pivot = pd.pivot(full, index = 'names', columns = ['group'], values = ['logfoldchanges','pvals_adj', 'scores'])

        full_pivot.to_csv(full_marker_output_path_scvi)



    marker_genes = pd.DataFrame(org_adata.uns['rank_genes_groups']['names']).iloc[:1000, :]
    new_colnames = [f'{col}_leiden_{resolution}' for col in marker_genes.columns]
    marker_genes.columns = new_colnames
    marker_dict[resolution] = marker_genes.T

# save marker genes per cluster
all_markers = pd.concat(marker_dict.values()).T


all_markers.to_csv(marker_output_path_scvi)

# OCEAN data was SoupX corrected, following steps cannot all handle floats
org_adata.layers["counts"] = np.round(org_adata.layers["counts"])
org_adata.obs['sample_study'] = org_adata.obs['sample'].astype(str) + '_' + org_adata.obs['study'].astype(str)

# also make pseudobulks
padata = dc.get_pseudobulk(org_adata, sample_col='sample_study', groups_col='cluster_organ_to_use', layer='counts', min_prop=0, min_smpls=0)
pb_dat = pd.DataFrame(padata.X)
pb_dat.columns = padata.var.index.values
pb_dat.index = padata.obs.index.values
# save pseudobulks
pb_dat.to_csv(pb_data)

