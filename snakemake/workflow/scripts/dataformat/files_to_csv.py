from pathlib import Path
import pandas as pd


def _load_table(path: Path) -> pd.DataFrame:
    """Load parquet or pickled dataframe-like content into a DataFrame."""
    suffix = path.suffix.lower()
    if suffix in {".parquet", ".parq"}:
        return pd.read_parquet(path, engine="pyarrow")

    if suffix in {".pckl", ".pkl", ".pickle"}:
        obj = pd.read_pickle(path)
    else:
        raise ValueError(f"Unsupported file type for {path}")

    # Direct DataFrame
    if isinstance(obj, pd.DataFrame):
        return obj

    # Dict cases
    if isinstance(obj, dict):
        values = list(obj.values())

        # Dict of DataFrames -> concatenate; ensure celltype column present
        if values and all(isinstance(v, pd.DataFrame) for v in values):
            frames = []
            for key, df in obj.items():
                if "celltype" not in df.columns:
                    df = df.copy()
                    df["celltype"] = key
                # keep gene index as column if present
                if df.index.name and df.index.name not in df.columns:
                    df = df.reset_index()
                frames.append(df)
            return pd.concat(frames, axis=0, ignore_index=False)

        # Dict of scalars -> single-row DataFrame
        if all(
            not isinstance(v, (list, tuple, set, pd.Series, pd.Index, range, dict))
            for v in values
        ):
            return pd.DataFrame([obj])

        # Otherwise let pandas align lists/arrays by keys
        return pd.DataFrame(obj)

    # Iterable of records (list/tuple of dicts) or array-like
    return pd.DataFrame(obj)


def main():
    input_paths = [Path(p) for p in snakemake.input]  # type: ignore[name-defined]
    output_paths = [Path(p) for p in snakemake.output]  # type: ignore[name-defined]

    if len(input_paths) != len(output_paths):
        raise ValueError("Number of outputs must match number of inputs")

    for in_path, out_path in zip(input_paths, output_paths):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        df = _load_table(in_path)
        if in_path.stem in {"enrich_collectri", "enrich_difib_collectri"} and "FC" in df.columns:
            df = df.rename(columns={"FC": "enrichment"})
        df.to_csv(out_path, index=False)


if __name__ == "__main__":
    main()
