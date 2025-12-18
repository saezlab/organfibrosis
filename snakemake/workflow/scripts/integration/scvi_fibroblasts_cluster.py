"""Run UMAP and Leiden clustering on pretrained scVI embeddings."""

import argparse
import os
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
import torch
from rich import print
from typing import Dict

SCVI_LATENT_KEY = "X_scVI"
NEIGHBORS_KEY = "neighbors_scvi"
LEIDEN_RESOLUTIONS = (0.2, 0.3, 0.4, 0.5, 0.7)


def configure_deterministic_seed(seed: int) -> None:
    """Set seeds across all libraries to make runs deterministic."""
    # Ensure cuBLAS runs deterministically on CUDA >= 10.2.
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    scvi.settings.seed = seed
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    try:
        torch.use_deterministic_algorithms(True)
    except (AttributeError, RuntimeError):
        pass


def load_trained_organs(intermediate_dir: Path) -> Dict[str, sc.AnnData]:
    """Load per-organ AnnData objects produced during the training stage."""
    if not intermediate_dir.exists():
        raise FileNotFoundError(f"Intermediate directory {intermediate_dir} not found")

    organ_adatas: Dict[str, sc.AnnData] = {}
    for path in sorted(intermediate_dir.glob("*_scvi.h5ad")):
        organ = path.stem.removesuffix("_scvi")
        organ_adatas[organ] = sc.read_h5ad(str(path))
        print(f"[bold]Loaded trained organ AnnData:[/bold] {organ} from {path}")

    if not organ_adatas:
        raise ValueError(
            f"No trained AnnData files found in {intermediate_dir}. "
            "Expected files named <organ>_scvi.h5ad."
        )
    return organ_adatas


def run_embedding_and_clustering(
    organ_adatas: Dict[str, sc.AnnData],
    output_csv: Path,
    output_plot_path: Path,
    seed: int,
) -> None:
    """Compute neighbors, UMAP, Leiden clusters, plots, and CSV summaries."""
    output_plot_path.mkdir(parents=True, exist_ok=True)
    harm_dict: Dict[str, pd.DataFrame] = {}

    for organ, organ_adata in organ_adatas.items():
        print(f"[bold]Computing embeddings for organ:[/bold] {organ}")
        if SCVI_LATENT_KEY not in organ_adata.obsm:
            raise KeyError(
                f"{SCVI_LATENT_KEY} missing in obsm for organ {organ}. "
                "Run the training stage first."
            )

        sc.pp.neighbors(
            organ_adata,
            use_rep=SCVI_LATENT_KEY,
            key_added=NEIGHBORS_KEY,
            random_state=seed,
        )

        sc.tl.umap(organ_adata, neighbors_key=NEIGHBORS_KEY, random_state=seed)
        for resolution in LEIDEN_RESOLUTIONS:
            sc.tl.leiden(
                organ_adata,
                neighbors_key=NEIGHBORS_KEY,
                resolution=resolution,
                key_added=f"leiden_{resolution}",
                random_state=seed,
            )

        fig = sc.pl.umap(
            organ_adata,
            color=[
                "study",
                "leiden_0.2",
                "leiden_0.3",
                "leiden_0.4",
                "leiden_0.5",
                "leiden_0.7",
                "cond_test",
                "POSTN",
                "annotation_MOFA",
            ],
            frameon=False,
            ncols=1,
            show=False,
            return_fig=True,
        )
        plot_path = output_plot_path / f"{organ}.pdf"
        fig.savefig(plot_path, bbox_inches="tight", dpi=200)
        plt.close(fig)

        harm = pd.DataFrame(
            organ_adata.obsm[SCVI_LATENT_KEY],
            index=organ_adata.obs.index,
        )
        harm["id"] = organ_adata.obs.index
        harm[["UMAP1", "UMAP2"]] = organ_adata.obsm["X_umap"]
        harm[["leiden_0.2","leiden_0.3","leiden_0.4", "leiden_0.5", "leiden_0.7"]] = organ_adata.obs[
            ["leiden_0.2","leiden_0.3","leiden_0.4", "leiden_0.5", "leiden_0.7"]
        ].to_numpy()
        harm_dict[organ] = harm

    full_harm = pd.concat(harm_dict, names=["organ"])
    full_harm.to_csv(output_csv)
    print(f"Saved embeddings and clustering summary to {output_csv}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--intermediate",
        required=True,
        help="Directory containing trained per-organ AnnData (<organ>_scvi.h5ad).",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="CSV output path for embeddings and clustering summary.",
    )
    parser.add_argument(
        "-p",
        "--plot",
        required=True,
        help="Directory for per-organ UMAP plots.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=7,
        help="Global random seed for deterministic execution.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    sns.set_theme()
    sc.set_figure_params(figsize=(6, 6), frameon=False)
    torch.set_float32_matmul_precision("high")

    configure_deterministic_seed(args.seed)

    intermediate_dir = Path(args.intermediate)
    organ_adatas = load_trained_organs(intermediate_dir)
    output_csv = Path(args.output)
    output_plot_path = Path(args.plot)

    run_embedding_and_clustering(
        organ_adatas,
        output_csv=output_csv,
        output_plot_path=output_plot_path,
        seed=args.seed,
    )
    print("Clustering stage complete.")


if __name__ == "__main__":
    main()
