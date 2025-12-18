import pandas as pd

# snakmake inputs
spatial_path = snakemake.input["data_path"]

# snakemake params
organ_mapping = snakemake.params["organ_mapping"]

# snakmake outputs
out_path = snakemake.output["out_path"]

spatial_df = pd.read_csv(spatial_path, index_col=0, header=0)

df = spatial_df[
    [
        "sample",
        # "interaction",
        "ligand",
        "morans",
        "mean",
        "std",
        "cond_test",
        "organ",
    ]
]
df = df.rename(
    columns={
        "mean": "mean_cosine_similarity",
        "std": "std_cosine_similarity",
        "ligand": "gene",
    }
)
df["organ"] = df["organ"].map(organ_mapping)

df.to_pickle(out_path)
