from pathlib import Path

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import liana as li

# Snakemake I/O
ccc_path = snakemake.input["ccc"]
geneset_paths = snakemake.input["genesets"]
adata_path = snakemake.input["adata"]

cosine_out = snakemake.output["cosine"]
cosine_extra_out = snakemake.output["cosine_extra"]
visium_plot_out = snakemake.output["visium_plot"]

# Parameters
organs = snakemake.params["organs"]
plot_genes = snakemake.params["plot_genes"]
sample_var = snakemake.params["sample_var"]
condition_var = snakemake.params["condition_var"]
exclude_samples = snakemake.params["exclude_samples"]

print(exclude_samples)
condition_map = {
        "Living Donor Healthy Reference Tissue (KPMP)": "control",
        "Reference Tissue (Unknown Clinical Status)": "control",
        "Healthy Reference Tissue": "control",
        "Diabetes Mellitus - Resilient": "control",
        "Hypertensive Chronic Kidney Disease": "CKD",
        "Diabetic Kidney Disease": "CKD",
        "Acute Kidney Injury": "AKI",
    }
plot_samples = snakemake.params["plot_samples"]


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


# Build interactions
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

# Load data and preprocess
adata = sc.read_h5ad(adata_path)
adata.layers["counts"] = adata.X.copy()
sc.pp.filter_cells(adata, min_genes=400)
sc.pp.filter_genes(adata, min_cells=5)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

dataset_names = list(adata.obs[sample_var].unique())


for pair_list in [sorted_top_pairs_tuple_all, extra_pairs]:
    lrdict = {}
    for dataset_name in dataset_names:
        adata_sub = adata[adata.obs[sample_var] == dataset_name]
        plot, bandwidth_df = li.ut.query_bandwidth(
            coordinates=adata_sub.obsm["spatial"], start=0, end=600, interval_n=20
        )
        bandwidth = bandwidth_df[bandwidth_df["neighbours"] > 5]["bandwith"].iloc[0]
        condition = adata_sub.obs[condition_var][0]

        li.ut.spatial_neighbors(adata_sub, bandwidth=bandwidth, set_diag=True, cutoff=0.1)
        lr = li.mt.bivariate(
            adata_sub,
            local_name="cosine",
            global_name="morans",
            nz_prop=0,
            n_perms=None,
            use_raw=False,
            add_categories=False,
            interactions=pair_list,
        ).var
        lr["cond_test"] = condition
        lrdict[dataset_name] = lr

    all_lr_results = pd.concat(lrdict).reset_index().rename(columns = {'level_0':'sample'})
    all_lr_results["cond_test"] = all_lr_results["cond_test"].replace(condition_map)
    print(all_lr_results.shape)
    all_lr_results = all_lr_results[~all_lr_results["sample"].isin(exclude_samples)]
    print(all_lr_results.shape)
    if pair_list == sorted_top_pairs_tuple_all:
        Path(cosine_out).parent.mkdir(parents=True, exist_ok=True)
        all_lr_results.to_csv(cosine_out)
    else:
        Path(cosine_extra_out).parent.mkdir(parents=True, exist_ok=True)
        all_lr_results.to_csv(cosine_extra_out)


# Plot example visium slide
fibr_sample = plot_samples["fib"]
ref_sample = plot_samples["ref"]
plot_sin = plot_genes
joined = "_".join(plot_sin)

fig, ax = plt.subplots(2, 2, figsize=(9, 8), tight_layout=True)

for count, slide in enumerate([ref_sample, fibr_sample]):
    adata_sub = adata[adata.obs[sample_var] == slide]
    condition = adata_sub.obs[condition_var][0]
    cond_test = "reference" if condition_map[condition] == "control" else "fibrosis"
    gene_list_present = list(set(adata.var_names) & set(plot_sin))

    for genecount, gene in enumerate(gene_list_present):
        sc.pl.spatial(
            adata_sub,
            color=gene,
            library_id=slide,
            show=False,
            size=1.3,
            img_key=None,
            title=f"{cond_test} {gene}",
            use_raw=False,
            ax=ax[count, genecount],
        )

Path(visium_plot_out).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(visium_plot_out)
