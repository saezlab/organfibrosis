# Author: Leonie Küchenhoff

# In this script, each patient is assigned to a fibrosis score,
# based on the enrichment of core matrisome genes (https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/NABA_CORE_MATRISOME.html)

import scanpy as sc
import decoupler as dc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


pb_path = snakemake.input[0]
meta_path = snakemake.input[1]


filter_by_expr = snakemake.params[0].get("filter_by_expr")
filter_by_prop = snakemake.params[0].get("filter_by_prop")
naba_path = snakemake.params[1]

out_path = snakemake.output[0]


def load_gmt_to_dataframe(gmt_file_path):
    """
    Load a GMT file and convert it to a DataFrame.
    """
    data = []  
    with open(gmt_file_path, "r") as file:
        for line in file:
            parts = line.strip().split("\t")  # Splitting each line by tab
            gene_set_name = parts[0]
            description = parts[1]
            genes = list(parts[2:])  # The rest of the parts are genes
        # Convert the list of dictionaries to a DataFrame
        df = pd.DataFrame(genes, columns=["genesymbol"])
        df["geneset"] = gene_set_name
        df["description"] = description
    return df


naba = load_gmt_to_dataframe(naba_path)


pdata = sc.read_csv(pb_path)
metadata = pd.read_csv(meta_path)

index_names = pdata.obs.index.str.split("_")[:].tolist()

pdata.obs["ctype"] = [item[-1] for item in index_names]
pdata.obs["sample"] = [
    "_".join(item[:-1]) if len(item) > 1 else item[1] for item in index_names
]
pdata.obs = pdata.obs.merge(metadata, on="sample", how="left")

pdata.obs.index = pdata.obs["sample"]
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
ctype_pdata = pdata[pdata.obs["ctype"] == "mesenchymal"].copy()


# Run ULM (Univariate Linear Model) enrichment for the mesenchymal cell type
dc.run_ulm(
    ctype_pdata,
    net=naba,
    source="geneset",
    target="genesymbol",
    use_raw=False,
    weight=None,
    verbose=True,
)

acts = dc.get_acts(ctype_pdata, obsm_key="ulm_estimate")

# Replace infinite values in the activity matrix with the maximum finite value
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

# Prepare output column names (including the NABA score)
col_names = list(acts.obs.columns)
col_names.append("NABA_CORE_MATRISOME")

# Extract scores for all samples and save to CSV
scores = sc.get.obs_df(acts, keys=col_names)

scores.to_csv(out_path)
