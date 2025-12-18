from pathlib import Path

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import liana as li

# Snakemake I/O
ccc_path = snakemake.input["ccc"]
geneset_paths = snakemake.input["genesets"]
visium_files = snakemake.input["visium_files"]

cosine_out = snakemake.output["cosine"]
cosine_extra_out = snakemake.output["cosine_extra"]
visium_plot_out = snakemake.output["visium_plot"]

# Parameters
organs = snakemake.params["organs"]
slides = snakemake.params["slides"]
visium_dir = snakemake.params["visium_dir"]
plot_genes = snakemake.params["plot_genes"]
cond_map = {"primary sclerosing cholangitis": "fibrosis"}



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

for pair_list in [sorted_top_pairs_tuple_all, extra_pairs]:
    lrdict = {}
    for count, slide in enumerate(slides):
        adata = sc.read(Path(visium_dir) / f"{slide}.h5ad")
        adata.var_names = adata.var["feature_name"].str.rsplit("_", n=1).str[0]
        adata.var.index.names = ["index"]
        sc.pp.filter_cells(adata, min_genes=400)
        sc.pp.filter_genes(adata, min_cells=5)

        region = adata.obs['disease'][0]

        plot, _ = li.ut.query_bandwidth(coordinates=adata.obsm['spatial'], start=0, end=600, interval_n=20)
        bandwidth = _[_['neighbours'] > 5]['bandwith'].iloc[0]
        print(bandwidth)
        
        li.ut.spatial_neighbors(adata, bandwidth=bandwidth, set_diag=True, cutoff=0.1)
        lr = li.mt.bivariate(adata,
                    local_name='cosine',
                    global_name="morans",
                    n_perms=None,
                    use_raw=False,
                    add_categories=False,
                    interactions=pair_list,
                    ).var
        lr['cond_test'] = region

        lrdict[slide] = lr

    all_lr_results = pd.concat(lrdict).reset_index().rename(columns = {'level_0':'sample'})
    all_lr_results.loc[all_lr_results['cond_test'] == 'primary sclerosing cholangitis', 'cond_test'] = 'PSC'

    if pair_list == sorted_top_pairs_tuple_all:
        Path(cosine_out).parent.mkdir(parents=True, exist_ok=True)
        all_lr_results.to_csv(cosine_out)
    else:
        Path(cosine_extra_out).parent.mkdir(parents=True, exist_ok=True)
        all_lr_results.to_csv(cosine_extra_out)


# Plot example visium slides
fibr_sample = snakemake.params.get("plot_samples", {}).get("fib", "Visium_PSC011_A1")
ref_sample = snakemake.params.get("plot_samples", {}).get("ref", "Visium_C73_D1")

plot_sin = plot_genes
joined = "_".join(plot_sin)

fig, ax = plt.subplots(2, 2, figsize=(9, 8), tight_layout=True)

for count, slide in enumerate([ref_sample, fibr_sample]):
    adata = sc.read(Path(visium_dir) / f"{slide}.h5ad")
    sc.pp.filter_cells(adata, min_genes=400)
    sc.pp.filter_genes(adata, min_cells=5)
    adata.var_names = adata.var["feature_name"].str.rsplit("_", n=1).str[0]
    adata.var.index.names = ["index"]

    condition = adata.obs["disease"][0]
    cond_test = "fibrosis" if cond_map.get(condition, condition) == "fibrosis" else "reference"
    gene_list_present = list(set(adata.var_names) & set(plot_sin))

    for genecount, gene in enumerate(gene_list_present):
        sc.pl.spatial(
            adata,
            color=gene,
            library_id=list(adata.uns["spatial"].keys())[0],
            show=False,
            size=1.3,
            img_key=None,
            title=f"{cond_test} {gene}",
            use_raw=False,
            ax=ax[count, genecount],
        )

Path(visium_plot_out).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(visium_plot_out)
