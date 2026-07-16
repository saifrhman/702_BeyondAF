from pathlib import Path
import argparse
import re
import subprocess

import pandas as pd


MODEL_NAMES = {
    "alphafold_params": "AlphaFold",
    "openfold_official": "OpenFold",
    "trained_model_epoch_32": "Retrained_epoch_32",
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--usalign", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = args.out_dir / "raw_outputs"
    raw_dir.mkdir(exist_ok=True)

    prediction_files = sorted(
        args.root.rglob("2olo_A_model_1_unrelaxed.pdb")
    )

    if len(prediction_files) != 15:
        raise RuntimeError(
            f"Expected 15 predictions, found {len(prediction_files)}"
        )

    rows = []

    for prediction_path in prediction_files:
        relative = prediction_path.relative_to(args.root)
        model_directory = relative.parts[0]
        model = MODEL_NAMES[model_directory]

        repeat_match = re.search(
            r"repeat_(\d+)",
            relative.as_posix(),
        )
        if not repeat_match:
            raise RuntimeError(
                f"Cannot identify repeat: {prediction_path}"
            )

        repeat = int(repeat_match.group(1))

        result = subprocess.run(
            [
                str(args.usalign),
                str(prediction_path),
                str(args.reference),
                "-mol",
                "prot",
            ],
            text=True,
            capture_output=True,
            check=True,
        )

        output = result.stdout

        raw_path = raw_dir / (
            f"{model_directory}_repeat_{repeat}.txt"
        )
        raw_path.write_text(output)

        length_matches = re.findall(
            r"Length of Structure_[12]:\s+(\d+)",
            output,
        )
        alignment_match = re.search(
            r"Aligned length=\s*(\d+),\s*"
            r"RMSD=\s*([\d.]+),\s*"
            r"Seq_ID=n_identical/n_aligned=\s*([\d.]+)",
            output,
        )
        tm_matches = re.findall(
            r"TM-score=\s*([\d.]+)",
            output,
        )

        if (
            len(length_matches) != 2
            or alignment_match is None
            or len(tm_matches) < 2
        ):
            raise RuntimeError(
                f"Could not parse US-align output for {prediction_path}"
            )

        rows.append(
            {
                "model": model,
                "repeat": repeat,
                "seed": 1000 + repeat,
                "prediction_length": int(length_matches[0]),
                "reference_length": int(length_matches[1]),
                "aligned_length": int(alignment_match.group(1)),
                "usalign_rmsd_angstrom": float(
                    alignment_match.group(2)
                ),
                "sequence_identity": float(
                    alignment_match.group(3)
                ),
                "tm_score_prediction_normalized": float(
                    tm_matches[0]
                ),
                "tm_score_reference_normalized": float(
                    tm_matches[1]
                ),
                "prediction_path": str(prediction_path),
                "raw_output_path": str(raw_path),
            }
        )

    result_df = pd.DataFrame(rows)

    order = {
        "AlphaFold": 0,
        "OpenFold": 1,
        "Retrained_epoch_32": 2,
    }

    result_df["_order"] = result_df["model"].map(order)
    result_df = (
        result_df.sort_values(["_order", "repeat"])
        .drop(columns="_order")
        .reset_index(drop=True)
    )

    summary_df = (
        result_df.groupby("model", sort=False)
        .agg(
            predictions=("tm_score_reference_normalized", "count"),
            mean_tm_score=("tm_score_reference_normalized", "mean"),
            std_tm_score=("tm_score_reference_normalized", "std"),
            min_tm_score=("tm_score_reference_normalized", "min"),
            max_tm_score=("tm_score_reference_normalized", "max"),
            mean_usalign_rmsd=(
                "usalign_rmsd_angstrom",
                "mean",
            ),
            mean_aligned_length=("aligned_length", "mean"),
        )
        .reset_index()
    )

    result_df.to_csv(
        args.out_dir / "2olo_tmscore_per_prediction.csv",
        index=False,
    )
    summary_df.to_csv(
        args.out_dir / "2olo_tmscore_model_summary.csv",
        index=False,
    )

    print(f"Predictions processed: {len(result_df)}")
    print()
    print(summary_df.to_string(index=False))
    print()
    print(f"Saved outputs to: {args.out_dir}")


if __name__ == "__main__":
    main()
