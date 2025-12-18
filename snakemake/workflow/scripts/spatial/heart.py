import os
from pathlib import Path

import pandas as pd
import scanpy as sc
import decoupler as dc  # noqa: F401  (kept for environment consistency)
import matplotlib.pyplot as plt
import seaborn as sns  # noqa: F401  (kept for environment consistency)
import pickle
import liana as li
from liana.method.sp import RandomForestModel, LinearModel, RobustLinearModel  # noqa: F401
import plotnine as p9  # noqa: F401

ccc_path = snakemake.input["ccc"]
geneset_paths = snakemake.input["genesets"]
metadata_path = snakemake.input["metadata"]

cosine_out = snakemake.output["cosine"]
cosine_extra_out = snakemake.output["cosine_extra"]
visium_plot_out = snakemake.output["visium_plot"]

organs = snakemake.params["organs"]
visium_dir = snakemake.params["visium_dir"]
plot_genes = snakemake.params["plot_genes"]


def load_gmt_to_dataframe(gmt_file_path):
    data = []  # To store each row of data
    with open(gmt_file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')  # Splitting each line by tab
            gene_set_name = parts[0]
            description = parts[1]
            genes = list(parts[2:])  # The rest of the parts are genes
        # Convert the list of dictionaries to a DataFrame
        df = pd.DataFrame(genes, columns = ['genesymbol'])
        df['geneset'] = gene_set_name
        df['description'] = description 
    return df

def concatenate_gmt_files(gmt_files_list):
    dfs = []  # To store DataFrames loaded from each file
    for file_path in gmt_files_list:
        df = load_gmt_to_dataframe(file_path)
        dfs.append(df)
    # Concatenate all DataFrames into one
    concatenated_df = pd.concat(dfs, ignore_index=True)
    return concatenated_df

def clean_list(interactions_to_plot):
    interactions_to_plot_clean = []
    for interactor, col in interactions_to_plot:
        if '_' in interactor:
            print(interactor)
            
            interactors = interactor.split('_')
            for i in interactors:
                pair = (i, col)
                interactions_to_plot_clean.append(pair)
        else:
            interactions_to_plot_clean.append((interactor, col))

    return(interactions_to_plot_clean)



ccc = pd.read_pickle(ccc_path)
results_list = [ccc[i][(ccc[i]['interaction_eff'] > 0.5) & (ccc[i]['interaction_sd_eff'] > 0)]
                 for i in organs]
all_results = pd.concat(results_list)
sorted_top_pairs_tuple_all = list(set(list(zip(all_results.ligand, all_results.receptor))))




resource = li.resource.select_resource(resource_name='consensus')
ligands = resource['ligand'].unique()
receptors =  resource['receptor'].unique()
ecm = concatenate_gmt_files(geneset_paths)
ligands_filtered = set(ligands) - set(ecm['genesymbol'])
receptors_filtered = set(receptors) - set(ecm['genesymbol'])

all_filtered = ligands_filtered | receptors_filtered
extra_pairs = list(set(list(zip(all_filtered, ['COL1A1' for i in range(len(all_filtered))]))))
extra_pairs = clean_list(extra_pairs)


metadata = pd.read_csv(metadata_path)

for pair_list in [sorted_top_pairs_tuple_all, extra_pairs]:
    lrdict = {}
    for count, slide in enumerate(metadata['slide_name']):
        adata = sc.read(Path(visium_dir) / f"{slide}.h5ad")
        sc.pp.filter_cells(adata, min_genes=100)
        sc.pp.filter_genes(adata, min_cells=5)
        
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        region = metadata['major_labl'][metadata['slide_name'] == slide].item()

        plot, _ = li.ut.query_bandwidth(coordinates=adata.obsm['spatial'], start=0, end=600, interval_n=20)
        bandwidth = _[_['neighbours'] > 5]['bandwith'].iloc[0]
        print(bandwidth)

        li.ut.spatial_neighbors(adata, bandwidth=bandwidth, set_diag=True, cutoff=0.1)
        lr = li.mt.bivariate(adata,
                    local_name='cosine', # Name of the function
                    global_name="morans",
                    nz_prop=0.05,
                    n_perms=None,
                    use_raw=False,
                    add_categories=False,
                    interactions = pair_list
                    ).var
        lr['cond_test'] = region

        lrdict[slide] = lr

    all_lr_results = pd.concat(lrdict).reset_index().rename(columns = {'level_0':'sample'})

    if pair_list == sorted_top_pairs_tuple_all:
        Path(cosine_out).parent.mkdir(parents=True, exist_ok=True)
        all_lr_results.to_csv(cosine_out)
    else:
        Path(cosine_extra_out).parent.mkdir(parents=True, exist_ok=True)
        all_lr_results.to_csv(cosine_extra_out)





# plot example visium slide

fibr_sample = 'Visium_12_CK290'
ref_sample =  'Visium_1_CK279'

plot_sin = plot_genes
joined = '_'.join(plot_sin)

fig, ax = plt.subplots(2,2, figsize = (9,8), tight_layout = True)

for count, slide in enumerate([ref_sample, fibr_sample]):
     adata = sc.read(Path(visium_dir) / f"{slide}.h5ad")
     sc.pp.filter_cells(adata, min_genes=100)
     sc.pp.filter_genes(adata, min_cells=5)
    
     sc.pp.normalize_total(adata, target_sum=1e4)
     sc.pp.log1p(adata)

     region = metadata['major_labl'][metadata['slide_name'] == slide].item()
     if region == 'FZ':
          cond_test = 'fibrosis'
     else:
          cond_test = 'reference'
     gene_list_present = list(set(adata.var_names) & set(plot_sin))
     print(gene_list_present)

     for genecount, gene in enumerate(gene_list_present):
          sc.pl.spatial(adata, color=gene, library_id=list(adata.uns['spatial'].keys())[0], show = False,
               size=1.3, img_key=None, title = f'{cond_test} {gene}', use_raw = False, ax = ax[count, genecount])

Path(visium_plot_out).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(visium_plot_out)
