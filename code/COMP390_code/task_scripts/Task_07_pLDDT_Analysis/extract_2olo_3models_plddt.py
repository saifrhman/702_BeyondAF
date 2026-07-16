from pathlib import Path
import argparse
import re
import pandas as pd


MODEL_NAMES = {
    "alphafold_params": "AlphaFold",
    "openfold_official": "OpenFold",
    "trained_model_epoch_32": "Retrained_epoch_32",
}

PATH_PATTERN = re.compile(
    r"(?P<model_dir>alphafold_params|openfold_official|trained_model_epoch_32)"
    r"/repeat_(?P<repeat>\d+)"
    r"/predictions/(?P<target>.+?)_model_1_unrelaxed\.pdb$"
)


def extract_ca_plddt(pdb_path: Path):
    residues = []

    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue

            if line[12:16].strip() != "CA":
                continue

            try:
                residues.append(
                    {
                        "chain_id": line[21].strip(),
                        "residue_number": int(line[22:26].strip()),
                        "insertion_code": line[26].strip(),
                        "residue_name": line[17:20].strip(),
                        "plddt": float(line[60:66].strip()),
                    }
                )
            except ValueError:
                continue

    return residues


def main():
    parser = argparse.ArgumentParser(
        description="Extract pLDDT from the 2OLO three-model repeated experiment."
    )
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    root = args.root.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb_files = sorted(root.rglob("*_model_1_unrelaxed.pdb"))

    matched_files = []
    for pdb_path in pdb_files:
        relative_path = pdb_path.relative_to(root).as_posix()
        match = PATH_PATTERN.fullmatch(relative_path)

        if match:
            matched_files.append((pdb_path, match.groupdict()))

    if len(matched_files) != 15:
        raise RuntimeError(
            f"Expected 15 prediction files but found {len(matched_files)} "
            f"under {root}"
        )

    prediction_rows = []
    residue_rows = []

    for pdb_path, metadata in matched_files:
        residues = extract_ca_plddt(pdb_path)

        if not residues:
            raise RuntimeError(f"No CA pLDDT values found in {pdb_path}")

        model = MODEL_NAMES[metadata["model_dir"]]
        repeat = int(metadata["repeat"])
        target = metadata["target"]

        values = [row["plddt"] for row in residues]

        prediction_rows.append(
            {
                "model": model,
                "repeat": repeat,
                "seed": 1000 + repeat,
                "target": target,
                "mean_plddt": sum(values) / len(values),
                "min_plddt": min(values),
                "max_plddt": max(values),
                "num_residues": len(values),
                "pdb_path": str(pdb_path),
            }
        )

        for residue_index, residue in enumerate(residues, start=1):
            residue_rows.append(
                {
                    "model": model,
                    "repeat": repeat,
                    "seed": 1000 + repeat,
                    "target": target,
                    "residue_index": residue_index,
                    **residue,
                    "pdb_path": str(pdb_path),
                }
            )

    prediction_df = pd.DataFrame(prediction_rows)
    residue_df = pd.DataFrame(residue_rows)

    model_order = {
        "AlphaFold": 0,
        "OpenFold": 1,
        "Retrained_epoch_32": 2,
    }

    prediction_df["_order"] = prediction_df["model"].map(model_order)
    prediction_df = (
        prediction_df.sort_values(["_order", "repeat"])
        .drop(columns="_order")
        .reset_index(drop=True)
    )

    residue_df["_order"] = residue_df["model"].map(model_order)
    residue_df = (
        residue_df.sort_values(["_order", "repeat", "residue_index"])
        .drop(columns="_order")
        .reset_index(drop=True)
    )

    model_summary = (
        prediction_df.groupby("model", sort=False)
        .agg(
            predictions=("mean_plddt", "count"),
            mean_plddt=("mean_plddt", "mean"),
            repeat_std=("mean_plddt", "std"),
            lowest_run_mean=("mean_plddt", "min"),
            highest_run_mean=("mean_plddt", "max"),
        )
        .reset_index()
    )

    residue_summary = (
        residue_df.groupby(
            [
                "model",
                "target",
                "residue_index",
                "chain_id",
                "residue_number",
                "insertion_code",
                "residue_name",
            ],
            dropna=False,
            sort=False,
        )
        .agg(
            mean_plddt=("plddt", "mean"),
            repeat_std=("plddt", "std"),
            min_plddt=("plddt", "min"),
            max_plddt=("plddt", "max"),
            predictions=("plddt", "count"),
        )
        .reset_index()
    )

    prediction_csv = out_dir / "2olo_plddt_per_prediction.csv"
    model_csv = out_dir / "2olo_plddt_model_summary.csv"
    residue_csv = out_dir / "2olo_plddt_per_residue.csv"
    residue_summary_csv = out_dir / "2olo_plddt_residue_summary.csv"

    prediction_df.to_csv(prediction_csv, index=False)
    model_summary.to_csv(model_csv, index=False)
    residue_df.to_csv(residue_csv, index=False)
    residue_summary.to_csv(residue_summary_csv, index=False)

    print(f"Processed predictions: {len(prediction_df)}")
    print(f"Processed residue records: {len(residue_df)}")
    print()
    print(model_summary.to_string(index=False))
    print()
    print(f"Saved outputs to: {out_dir}")


if __name__ == "__main__":
    main()
