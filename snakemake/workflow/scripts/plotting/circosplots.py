import os

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from nichenetpy.wrappers import (
    create_ligand_receptor_links_prioritization_circos_plot,
)


ccc_results_path = snakemake.input["ccc_results"]
output_paths = snakemake.output["plots"]
ccc_table_path = snakemake.output["ccc_table"]

organs = snakemake.params["organs"]
real_names = snakemake.params["real_names"]
ctype_colors = snakemake.params["ctype_colors"]
eff_cutoff_pos = snakemake.params["eff_cutoff_pos"]
top_n = snakemake.params["top_n"]

ccc_results = pd.read_pickle(ccc_results_path)

plt.rcParams.update({"font.size": 20})

ccc_results_all = []

for organ, output_path in zip(organs, output_paths):
    ccc_result_organ = ccc_results[organ]
    ccc_results_all.append(ccc_result_organ.assign(organ=real_names[organ]))
    ctype_colors_org = {
        cell_type: ctype_colors[cell_type]
        for cell_type in ccc_result_organ["source"].unique()
    }
    significant = ccc_result_organ[
        (ccc_result_organ["interaction_eff"] > eff_cutoff_pos)
        & (
            ccc_result_organ["interaction_eff"]
            - ccc_result_organ["interaction_sd_eff"]
            > 0
        )
    ]
    significant = significant.sort_values(by="interaction_eff").iloc[-top_n:]

    ax, fig = create_ligand_receptor_links_prioritization_circos_plot(
        senders=significant["source"],
        receivers=significant["target"],
        ligands=significant["ligand"],
        receptors=significant["receptor"],
        colors=ctype_colors_org,
    )
    for text in fig.findobj(match=matplotlib.text.Text):
        text.set_fontsize(11)
    ax.set_title(f"{real_names[organ]} \n", fontsize=20)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)

# all organs' interactions in one long table, organ as an extra column
os.makedirs(os.path.dirname(ccc_table_path), exist_ok=True)
pd.concat(ccc_results_all, ignore_index=True).to_csv(ccc_table_path, index=False)
