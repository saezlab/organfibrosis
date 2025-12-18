"""Train scVI models for mesenchymal cells."""

import argparse
import os
import random
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
import torch
from rich import print
from typing import Dict, Optional

assert torch.cuda.is_available(), "No GPU found!"

SCVI_LATENT_KEY = "X_scVI"


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


def prepare_organs(input_path: str, hvg_dir: str) -> Dict[str, sc.AnnData]:
    """Load merged AnnData and prepare one AnnData per organ."""
    merged_adata = sc.read_h5ad(input_path)
    merged_adata.obs.index = (
        np.array(merged_adata.obs.index) + "_" + np.array(merged_adata.obs["study"])
    )
    merged_adata.obsm = pd.DataFrame()
    merged_adata.layers["counts"] = merged_adata.X.copy()
    merged_adata.obs["sample_study"] = (
        np.array(merged_adata.obs["sample"]) + "_" + np.array(merged_adata.obs["study"])
    )

    organ_adatas: Dict[str, sc.AnnData] = {}
    for organ in merged_adata.obs["organ"].unique():
        print(f"[bold]Preparing organ:[/bold] {organ}")
        organ_adata = merged_adata[merged_adata.obs["organ"] == organ].copy()
        sc.pp.filter_cells(organ_adata, min_genes=200)
        sc.pp.filter_genes(organ_adata, min_cells=10)

        organ_adata.var["rb"] = organ_adata.var_names.str.startswith(("RPS", "RPL"))
        organ_adata.var["mt"] = organ_adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(organ_adata, qc_vars=["rb", "mt"], inplace=True)

        hvg_path = Path(hvg_dir) / f"{organ}.csv"
        genes = list(
            pd.read_csv(hvg_path, index_col=0, skiprows=1, names=["genes"]).loc[
                :, "genes"
            ]
        )
        organ_adatas[organ] = organ_adata[:, genes].copy()
    return organ_adatas


def train_scvi_models(
    organ_adatas: Dict[str, sc.AnnData],
    intermediate_dir: Optional[Path],
    max_epochs: int = 50,
) -> Dict[str, sc.AnnData]:
    """Train scVI per organ and optionally write trained AnnData with latent space."""
    if intermediate_dir is not None:
        intermediate_dir.mkdir(parents=True, exist_ok=True)

    trained_adatas: Dict[str, sc.AnnData] = {}
    for organ, organ_adata in organ_adatas.items():
        print(f"[bold]Training scVI for organ:[/bold] {organ}")
        scvi.model.SCVI.setup_anndata(
            organ_adata,
            layer="counts",
            batch_key="sample_study",
            categorical_covariate_keys=["study", "modality", "tech"],
            continuous_covariate_keys=[
                "total_counts",
                "n_genes_by_counts",
                "pct_counts_mt",
                "pct_counts_rb",
            ],
        )
        model = scvi.model.SCVI(organ_adata, n_layers=2, n_latent=30)
        model.view_anndata_setup()
        model.train(max_epochs=max_epochs)

        organ_adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
        trained_adatas[organ] = organ_adata

        if intermediate_dir is not None:
            output_path = intermediate_dir / f"{organ}_scvi.h5ad"
            organ_adata.write_h5ad(str(output_path))
            print(f"Saved trained AnnData to {output_path}")
    return trained_adatas


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Path to merged AnnData", required=True)
    parser.add_argument("-v", "--hvg", help="Directory with organ HVGs", required=True)
    parser.add_argument(
        "--intermediate",
        required=True,
        help="Directory to store trained per-organ AnnData",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=7,
        help="Global random seed for deterministic execution.",
    )
    parser.add_argument(
        "--max-epochs",
        type=int,
        default=50,
        help="Maximum epochs for scVI training.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    sns.set_theme()
    sc.set_figure_params(figsize=(6, 6), frameon=False)
    torch.set_float32_matmul_precision("high")

    configure_deterministic_seed(args.seed)

    intermediate_dir = Path(args.intermediate)
    organ_adatas = prepare_organs(args.input, args.hvg)
    train_scvi_models(
        organ_adatas,
        intermediate_dir=intermediate_dir,
        max_epochs=args.max_epochs,
    )
    print("Training stage complete.")


if __name__ == "__main__":
    main()
