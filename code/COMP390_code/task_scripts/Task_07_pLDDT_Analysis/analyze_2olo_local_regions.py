from pathlib import Path
import argparse

import numpy as np
import pandas as pd


CENTRES = [38, 39, 268, 386]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--accuracy-csv", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--half-window", type=int, default=5)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.accuracy_csv)
    rows = []

    for (model, repeat), prediction in df.groupby(["model", "repeat"]):
        for centre in CENTRES:
            window = prediction[
                prediction["reference_residue_number"].between(
                    centre - args.half_window,
                    centre + args.half_window,
                )
            ]

            centre_row = prediction[
                prediction["reference_residue_number"] == centre
            ]

            if window.empty or centre_row.empty:
                raise RuntimeError(
                    f"Missing residue {centre} for {model}, repeat {repeat}"
                )

            distances = window["ca_error_angstrom"].to_numpy(float)

            rows.append({
                "model": model,
                "repeat": int(repeat),
                "centre_residue": centre,
                "window_start": centre - args.half_window,
                "window_end": centre + args.half_window,
                "residues_in_window": len(window),
                "centre_ca_error_angstrom": float(
                    centre_row.iloc[0]["ca_error_angstrom"]
                ),
                "window_mean_ca_error_angstrom": float(
                    distances.mean()
                ),
                "window_rmsd_angstrom": float(
                    np.sqrt(np.mean(distances ** 2))
                ),
                "window_max_ca_error_angstrom": float(
                    distances.max()
                ),
                "centre_plddt": float(
                    centre_row.iloc[0]["plddt"]
                ),
            })

    result = pd.DataFrame(rows)

    summary = (
        result.groupby(["model", "centre_residue"], sort=False)
        .agg(
            predictions=("repeat", "count"),
            mean_centre_error=(
                "centre_ca_error_angstrom",
                "mean",
            ),
            std_centre_error=(
                "centre_ca_error_angstrom",
                "std",
            ),
            mean_window_rmsd=(
                "window_rmsd_angstrom",
                "mean",
            ),
            std_window_rmsd=(
                "window_rmsd_angstrom",
                "std",
            ),
            mean_centre_plddt=("centre_plddt", "mean"),
        )
        .reset_index()
    )

    result.to_csv(
        args.out_dir / "2olo_local_regions_per_prediction.csv",
        index=False,
    )
    summary.to_csv(
        args.out_dir / "2olo_local_regions_summary.csv",
        index=False,
    )

    print(f"Local-region records: {len(result)}")
    print()
    print(summary.to_string(index=False))
    print()
    print(f"Saved outputs to: {args.out_dir}")


if __name__ == "__main__":
    main()
