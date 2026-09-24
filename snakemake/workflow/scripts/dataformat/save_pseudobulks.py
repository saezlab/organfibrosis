import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import scanpy as sc
import decoupler as dc

# Get input and output paths from snakemake
pb_path = snakemake.input[0]
meta_path = snakemake.input[1]
out_path = snakemake.output[0]

metadata = pd.read_csv(meta_path)
pdata = sc.read_csv(pb_path)

# Map organ codes to display names from config
organ_name_map = snakemake.config["general_plotting"]["organ_names"]


index_names = pdata.obs.index.str.split("_")[:].tolist()

pdata.obs["ctype"] = [item[-1] for item in index_names]
pdata.obs["sample"] = [
    "_".join(item[:-1]) if len(item) > 1 else item[1] for item in index_names
]
pdata.obs = pdata.obs.merge(metadata, left_on=["sample","ctype"], right_on=["sample","annotation_MOFA"], how="left")
pdata.obs.index = pdata.obs["sample"]
pdata.obs.index.name = None

# Also ensure organ column in obs uses the real names after merge
pdata.obs["organ"] = organ_name_map[snakemake.wildcards['organ']]
pdata.obs = pdata.obs.drop(columns=["Unnamed: 0", "counts"])
# Save AnnData to h5ad
pdata.write(out_path)
