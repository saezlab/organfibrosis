# Author: Leonie Küchenhoff

# In this script, the study by Muto et al. is preprocessed to that the format is comparable to other study.
# Anoter preprocessing step follows after this script.
import scanpy as sc
import pandas as pd
import anndata
import numpy as np


# define snakemake inputs and outputs
paths = snakemake.input["adata"]
input_metadata = snakemake.input["meta"]

output = snakemake.output[0]


patient_names = [path.split("_")[2] for path in paths]

print(paths)
print(patient_names)

# read in all h5 files and merge
adata_list = []
for path, patient in zip(paths, patient_names):
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    adata.obs["patient"] = patient
    adata.obs["id"] = adata.obs.index
    adata_list.append(adata)
# combine together
adata_merge = anndata.concat(adata_list, join="outer")


# merge with metadata
meta = pd.read_csv(input_metadata)
meta["patient"] = meta["patient"].replace(
    {
        "control1": "cont1",
        "control2": "cont2",
        "control3": "cont3",
        "control4": "cont4",
        "control5": "cont5",
    }
)

adata_merge.obs = adata_merge.obs.merge(
    meta, how="left", left_on=["patient", "id"], right_on=["patient", "barcode"]
)
adata_merge.obs.index = adata_merge.obs["name"]


# some cells did not include any metadata info
#  (most likely sorted out in original study. Will also be excluded here)
adata_merge = adata_merge[~adata_merge.obs["disease"].isna()]

adata_merge.write(output)
