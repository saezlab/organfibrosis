from pathlib import Path
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import liana as li
import anndata as ad
import decoupler as dc
from scipy import sparse
from itertools import product

# Snakemake I/O
ccc_path = snakemake.input["ccc"]
geneset_paths = snakemake.input["genesets"]
# per-sample spatial h5ad files written by spatial_pp/dot.py: each already carries the DOT
# deconvolution weights and the (organ-wide) niche labels in .obs, raw counts in
# .X, and the deconv cell types / niche names in .uns.
# Snakemake waits for the whole DOT directory before discovering its samples.
deconv_h5ads = sorted(Path(snakemake.input["deconv_dir"]).glob("*.h5ad"))
if not deconv_h5ads:
    raise ValueError(f"No DOT H5AD files found in {snakemake.input['deconv_dir']}")
extra_genesets = snakemake.input['extra_genesets']
progeny_genesets = snakemake.input['progeny']

cosine_out = snakemake.output["cosine"]
cosine_extra_out = snakemake.output["cosine_extra"]
visium_plot_out = snakemake.output["visium_plot"]
cosine_genesets_out = snakemake.output['cosine_genesets_out']
ulm_estimate_out = snakemake.output['ulm_estimate']
full_deconv_output = snakemake.output['deconv_output']
cosine_out_ctypes = snakemake.output['cosine_out_ctypes']

# Parameters
organs = snakemake.params["organs"]
plot_genes = snakemake.params["plot_genes"]
plot_enrichment = snakemake.params.get("plot_enrichment", "NABA_CORE_MATRISOME")
organ_name = snakemake.params.get("organ_name", "liver")
cond_map = {"primary sclerosing cholangitis": "fibrosis"}
plt.rcParams.update({"font.size": 20})


def load_gmt_to_dataframe(gmt_file_path):
    data = []
    with open(gmt_file_path, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            gene_set_name = parts[0]
            description = parts[1]
            genes = list(parts[2:])
    df = pd.DataFrame(genes, columns=["genesymbol"])
    df["geneset"] = gene_set_name
    df["description"] = description
    return df


def concatenate_gmt_files(gmt_files_list):
    dfs = []
    for file_path in gmt_files_list:
        dfs.append(load_gmt_to_dataframe(file_path))
    return pd.concat(dfs, ignore_index=True)


def clean_list(interactions_to_plot):
    interactions_to_plot_clean = []
    for interactor, col in interactions_to_plot:
        if "_" in interactor:
            interactors = interactor.split("_")
            for i in interactors:
                interactions_to_plot_clean.append((i, col))
        else:
            interactions_to_plot_clean.append((interactor, col))
    return interactions_to_plot_clean


# Interactions
ccc = pd.read_pickle(ccc_path)
gsets = pd.read_csv(extra_genesets, index_col = 0)
# add the progeny pathway footprints as extra genesets
progeny = pd.read_csv(progeny_genesets, index_col = 0)
progeny = progeny.reset_index().rename(columns = {'source':'geneset',
                                    'target':'genesymbol'}).drop(columns='padj')
progeny['collection'] = 'progeny'
gsets = pd.concat([progeny, gsets])
results_list = [
    ccc[i][(ccc[i]["interaction_eff"] > 0.5) & (ccc[i]["interaction_sd_eff"] > 0)]
    for i in organs
]
all_results = pd.concat(results_list)
sorted_top_pairs_tuple_all = list(set(zip(all_results.ligand, all_results.receptor)))

resource = li.resource.select_resource(resource_name="consensus")
ligands = resource["ligand"].unique()
receptors = resource["receptor"].unique()
ecm = concatenate_gmt_files(geneset_paths)
ligands_filtered = set(ligands) - set(ecm["genesymbol"])
receptors_filtered = set(receptors) - set(ecm["genesymbol"])
all_filtered = ligands_filtered | receptors_filtered
extra_pairs = list(set(zip(all_filtered, ["COL1A1" for _ in range(len(all_filtered))])))
extra_pairs = clean_list(extra_pairs)
interactions_enrichment = list(product(all_filtered, gsets['geneset'].unique()))

combined_pairs = extra_pairs + interactions_enrichment
combined_pairs = clean_list(combined_pairs)

# Deconvolution cell types and niche one-hot columns are identical across all
# samples of the organ (niches are clustered once, organ-wide, in spatial_pp/dot.py), so
# read them from the first file and colocalize every one of them to each geneset.
_probe = sc.read_h5ad(deconv_h5ads[0])
deconv_celltypes = list(_probe.uns['deconv_celltypes'])
niche_names = [str(n) for n in _probe.uns['niche_names']]
del _probe
# the niche one-hot columns are built per sample from the categorical obs['niche']
niche_cols = [f'niche_{n}' for n in niche_names]
feature_cols = deconv_celltypes + niche_cols
cell_enr_pairs = list(product(feature_cols, gsets['geneset'].unique()))

ulm_estimate_dfs = []
ref_sample = snakemake.params.get('plot_samples', {}).get('ref', 'Visium_C73_D1')
fibr_sample = snakemake.params.get('plot_samples', {}).get('fib', 'Visium_PSC011_A1')
plot_slides = [ref_sample, fibr_sample]
slide_to_row = {ref_sample: 0, fibr_sample: 1}
enrichment_name = plot_enrichment
plot_gene = plot_genes[0] if isinstance(plot_genes, (list, tuple)) and plot_genes else plot_genes
fig, ax = plt.subplots(2, 2, figsize=(9, 8), tight_layout=True)

for pair_list in [sorted_top_pairs_tuple_all, combined_pairs, cell_enr_pairs]:
    lrdict = {}
    full_deconv_df_list = []
    for h5ad_path in deconv_h5ads:
        slide = Path(h5ad_path).stem
        adata = sc.read_h5ad(h5ad_path)

        # raw counts are already in .X; one-hot encode the niche labels, then
        # normalize (deconv weights + niche columns are already in adata.obs and
        # stay row-aligned through cell/gene filtering)
        for _n, _c in zip(niche_names, niche_cols):
            adata.obs[_c] = (adata.obs['niche'].astype(str) == _n).astype(float)
        sc.pp.filter_cells(adata, min_genes=400)
        sc.pp.filter_genes(adata, min_cells=5)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        region = adata.obs['cond_test'][0]

        plot, _ = li.ut.query_bandwidth(coordinates=adata.obsm['spatial'], start=0, end=600, interval_n=20)
        bandwidth = _[_['neighbours'] > 5]['bandwith'].iloc[0]
        print(bandwidth)

        dc.run_ulm(
            mat=adata,
            net=gsets,
            source='geneset',
            target='genesymbol',
            weight='weight',
            verbose=True,
            use_raw=False
        )
        acts_enrichment = li.ut.obsm_to_adata(adata, 'ulm_estimate')

        if pair_list == sorted_top_pairs_tuple_all:
            ulm_df = pd.DataFrame(
                adata.obsm['ulm_estimate'],
                index=adata.obs_names,
            )
            ulm_df.insert(0, 'condition', region)
            ulm_df.insert(0, 'cell', adata.obs_names.to_numpy())
            ulm_df.insert(0, 'sample', slide)
            ulm_df.insert(0, 'organ', organ_name)
            # squidpy.ipynb prefixes every barcode with its sample name when it
            # concatenates the per-sample h5ads, so key the enrichment scores the
            # same way to make them joinable onto combined.obs_names there.
            ulm_df.index = pd.Index(
                [f'{slide}_{bc}' for bc in adata.obs_names], name='spot_id'
            )
            ulm_estimate_dfs.append(ulm_df)

            if slide in plot_slides:
                adata.obs[plot_enrichment] = adata.obsm['ulm_estimate'][plot_enrichment].values
                plot_order = [plot_enrichment, plot_gene]
                row = slide_to_row[slide]
                for col, color_key in enumerate(plot_order):
                    title = f"{region} {color_key}"
                    sc.pl.spatial(
                        adata,
                        color=color_key,
                        library_id=list(adata.uns['spatial'].keys())[0],
                        show=False,
                        size=1.3,
                        img_key=None,
                        title=title,
                        use_raw=False,
                        ax=ax[row, col],
                    )

        # Convert to sparse matrix
        new_X = sparse.csr_matrix(adata.obsm['ulm_estimate'].values)
        # cell-type proportions AND niche one-hot columns become extra "features"
        feature_matrix = sparse.csr_matrix(adata.obs.loc[:, feature_cols].astype(float).values)

        # Append to existing matrix
        X_combined = sparse.hstack([adata.X, new_X, feature_matrix]).tocsr()
        new_var = pd.DataFrame(index=adata.obsm['ulm_estimate'].columns)
        feature_var = pd.DataFrame(index=feature_cols)
        var_combined = pd.concat([adata.var, new_var, feature_var], axis=0)

        adata_sub_extended = ad.AnnData(
            X=X_combined,
            obs=adata.obs.copy(),
            var=var_combined,
            obsm=adata.obsm.copy()
        )

        full_deconv_df = adata.obs
        full_deconv_df[adata.obsm['ulm_estimate'].columns] = adata.obsm['ulm_estimate'].values
        full_deconv_df['spatial_x'] = adata.obsm['spatial'][:, 0]
        full_deconv_df['spatial_y'] = adata.obsm['spatial'][:, 1]
        full_deconv_df['sample'] = slide
        full_deconv_df_list.append(full_deconv_df)


        li.ut.spatial_neighbors(adata_sub_extended, bandwidth=bandwidth, set_diag=True, cutoff=0.1)
        lr = li.mt.bivariate(adata_sub_extended,
                    local_name='cosine',
                    global_name="morans",
                    n_perms=None,
                    use_raw=False,
                    add_categories=False,
                    interactions=pair_list,
                    ).var
        lr['cond_test'] = region

        lrdict[slide] = lr
    full_deconv_df_full = pd.concat(full_deconv_df_list)
    full_deconv_df_full.to_csv(full_deconv_output)

    all_lr_results = pd.concat(lrdict).reset_index().rename(columns = {'level_0':'sample'})
    all_lr_results.loc[all_lr_results['cond_test'] == 'primary sclerosing cholangitis', 'cond_test'] = 'PSC'

    if pair_list == sorted_top_pairs_tuple_all:
        Path(cosine_out).parent.mkdir(parents=True, exist_ok=True)
        all_lr_results.to_csv(cosine_out)
    elif pair_list == combined_pairs:
        Path(cosine_extra_out).parent.mkdir(parents=True, exist_ok=True)
        all_lr_results.to_csv(cosine_genesets_out)
        all_lr_results_col = all_lr_results[all_lr_results['receptor'] == 'COL1A1']
        all_lr_results_col.to_csv(cosine_extra_out)
    elif pair_list == cell_enr_pairs:
        Path(cosine_out_ctypes).parent.mkdir(parents=True, exist_ok=True)
        all_lr_results.to_csv(cosine_out_ctypes)

ulm_estimate_all = pd.concat(ulm_estimate_dfs)
Path(ulm_estimate_out).parent.mkdir(parents=True, exist_ok=True)
ulm_estimate_all.to_csv(ulm_estimate_out)


Path(visium_plot_out).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(visium_plot_out)
