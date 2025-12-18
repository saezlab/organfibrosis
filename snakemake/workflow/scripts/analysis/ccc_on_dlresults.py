# Author: Leonie Küchenhoff
# This script processes the results of a mixed effects model
# to extract cell-cell communication (CCC) data

import pickle
from itertools import product

import liana as li
import matplotlib.pyplot as plt
import pandas as pd
from liana.method._pipe_utils import filter_resource

# Parameters and thresholds
eff_cutoff_pos = 0.5
eff_cutoff_neg = -0.5

# Snakemake inputs/outputs/params
input_file = snakemake.input["organ_spec_dl"]
views = snakemake.params["views"]
organs = snakemake.params["organs"]
output_file = snakemake.output["ccc_result"]

# Plot parameters
plt.rcParams.update({"font.size": 18})

# Utility functions
def _join_stats(source, target, dedict, resource):
    """
    Joins and renames source-ligand and target-receptor stats to the ligand-receptor resource.

    Parameters
    ----------
    source
        Source/Sender cell type
    target
        Target/Receiver cell type
    dedict
        dictionary
    resource
        Ligand-receptor Resource

    Returns
    -------
    Ligand-Receptor stats

    """
    source_stats = dedict[source].copy()
    source_stats.columns = source_stats.columns.map(lambda x: "ligand_" + str(x))
    source_stats = source_stats.rename(
        columns={"ligand_gene": "ligand", "ligand_label": "source"}
    )

    target_stats = dedict[target].copy()
    target_stats.columns = target_stats.columns.map(lambda x: "receptor_" + str(x))
    target_stats = target_stats.rename(
        columns={"receptor_gene": "receptor", "receptor_label": "target"}
    )

    bound = resource.merge(source_stats).merge(target_stats)

    return bound


def filter_reassemble_complexes(
    lr_res,
    _key_cols,
    complex_cols,
    expr_prop,
    return_all_lrs=False,
    complex_policy="min",
):
    """
    Reassemble complexes from exploded long-format pandas DataFrame.

    Parameters
    ----------
    lr_res
        long-format pandas dataframe with exploded complex subunits
    _key_cols
        primary key for lr_res, typically a list with the following elements -
        ['source', 'target', 'ligand_complex', 'receptor_complex']
    complex_cols
        method/complex-relevant columns
    %(expr_prop)s
    %(return_all_lrs)s
    complex_policy
        approach by which the complexes are reassembled

    Return
    -----------
    lr_res: a reduced long-format pandas dataframe
    """

    aggs = {complex_policy, "min"}
    for col in complex_cols:
        lr_res = _reduce_complexes(
            col=col, lr_res=lr_res, key_cols=_key_cols, aggs=aggs
        )

    # Check for duplicated subunits
    duplicate_mask = lr_res.duplicated(subset=_key_cols, keep=False)
    if duplicate_mask.any():
        # Check for non-equal subunit values
        if (
            not lr_res[duplicate_mask]
            .groupby(_key_cols)[complex_cols]
            .transform(lambda x: x.duplicated(keep=False))
            .all()
            .all()
        ):
            _logg(
                "There were duplicated subunits in the complexes. "
                + "The subunits were reduced to only the minimum expression subunit. "
                + "However, there were subunits that were not the same within a complex. ",
                level="warn",
            )
        lr_res = lr_res.drop_duplicates(subset=_key_cols, keep="first")

    return lr_res


def _reduce_complexes(
    col: str, lr_res: pd.DataFrame, key_cols: list, aggs: (dict | str)
):
    """
    Reduce complex subunits by aggregation policy (e.g., min).

    Parameters
    ----------
    col
        Column to reduce
    lr_res
        DataFrame containing the data
    key_cols
        Key columns for grouping
    aggs
        Aggregation functions

    Returns
    -------
    Reduced DataFrame

    """
    lr_res = lr_res.groupby(key_cols)

    # Get min cols by which we will join
    # then rename from agg name to column name (e.g. 'min' to 'ligand_min')
    lr_min = (
        lr_res[col]
        .agg(aggs)
        .reset_index()
        .copy()
        .rename(columns={agg: col.split("_")[0] + "_" + agg for agg in aggs})
    )

    # right is the min subunit for that column
    join_key = col.split("_")[0] + "_min"  # ligand_min or receptor_min

    # Here, I join the min value and keep only those rows that match
    lr_res = lr_res.obj.merge(lr_min, on=key_cols, how="inner")
    lr_res = lr_res[lr_res[col] == lr_res[join_key]].drop(join_key, axis=1)

    return lr_res


def ccc_from_model(model_results, resource):
    """
    Assemble cell-cell communication results from model output and resource.

    Parameters
    ----------
    model_results
        Results from the mixed effects model
    resource
        Ligand-receptor resource

    Returns
    -------
    DataFrame containing cell-cell communication results

    """
    # Filter resource to only genes present in model results
    resource = filter_resource(
        resource, pd.concat(model_results)["gene"].unique().tolist()
    )

    stat_names = ["sd_eff", "eff"]
    complex_col = "eff"
    lr_sep = "^"

    # Join stats to ligand-receptor pairs
    lr_res = pd.concat(
        [
            _join_stats(source, target, model_results, resource)
            for source, target in zip(pairs["source"], pairs["target"])
        ]
    )

    # ligand_ or receptor + stat_keys
    complex_cols = list(product(["ligand", "receptor"], [complex_col]))
    complex_cols = [f"{x}_{y}" for x, y in complex_cols]

    # Assign receptor and ligand absolutes (handles missing values)
    _placeholders = ["ligand_absolute", "receptor_absolute"]
    lr_res[_placeholders] = lr_res[complex_cols].apply(lambda x: x.abs())

    lr_res = filter_reassemble_complexes(
        lr_res=lr_res,
        _key_cols=["source", "target", "ligand_complex", "receptor_complex"],
        expr_prop=0,
        return_all_lrs=False,
        complex_cols=_placeholders,
    )
    lr_res = lr_res.drop(["interaction", *_placeholders], axis=1)

    # Summarise stats for each ligand-receptor pair
    for key in stat_names:
        stat_columns = ["ligand_" + key, "receptor_" + key]
        lr_res.loc[:, f"interaction_{key}"] = lr_res.loc[:, stat_columns].mean(axis=1)

    lr_res["interaction"] = (
        lr_res["ligand_complex"] + lr_sep + lr_res["receptor_complex"]
    )

    lr_res["full_interaction"] = (
        lr_res["interaction"] + lr_res["source"] + "_" + lr_res["target"]
    )

    return lr_res


# --- Read in model results ---
with open(input_file, "rb") as fp:
    organ_spec = pickle.load(fp)

# Initialize a dictionary to hold the results for each organ
# Loop through organs and cell types to filter and prepare data
dedict = {}
for organ in organs:
    dedict[organ] = {}
    for ctype in views:
        # Prepare dataframe for each cell type and organ
        df = organ_spec[ctype][organ].reset_index().set_index("gene")
        # Keep genes where all w_re values are positive (after dropping NaNs)
        genes_to_keep = df.groupby(df.index)["w_re"].apply(
            lambda x: (x.dropna() > 0).all()
        )
        print(f"filtered out {len(df.index.unique()) - len(genes_to_keep)} genes.")
        working_model = df.loc[genes_to_keep]
        # Only keep rows corresponding to random effect summary
        working_model = working_model[working_model["summary_row"] == "random effect"]
        working_model["organ"] = organ
        working_model = working_model.loc[:, ["organ", "eff", "sd_eff"]]
        for_analysis = working_model.reset_index()
        for_analysis["label"] = ctype
        array_for_analysis = for_analysis[["gene", "label", "eff", "sd_eff"]]
        dedict[organ][ctype] = array_for_analysis

# Prepare resource and pairs
resource = li.resource.select_resource(resource_name="consensus")
pairs = pd.DataFrame(list(product(views, views))).rename(
    columns={0: "source", 1: "target"}
)
resource["interaction"] = resource["ligand"] + "&" + resource["receptor"]
resource = (
    resource.set_index("interaction")
    .apply(lambda x: x.str.split("_"))
    .explode(["ligand"])
    .explode("receptor")
    .reset_index()
)
resource[["ligand_complex", "receptor_complex"]] = resource["interaction"].str.split(
    "&", expand=True
)

# Run CCC analysis for each organ
ccc_results = {}
for organ in organs:
    ccc_organ = ccc_from_model(dedict[organ], resource)
    ccc_results[organ] = ccc_organ

# Save results 
with open(output_file, "wb") as fp:
    pickle.dump(ccc_results, fp)
