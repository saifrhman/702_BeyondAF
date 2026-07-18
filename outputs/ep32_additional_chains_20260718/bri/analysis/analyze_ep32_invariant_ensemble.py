#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FEATURE_GROUPS = {
    "BRI": [
        "x(AN)",
        "x(AC)",
        "y(AC)",
        "x(N)",
        "y(N)",
        "z(N)",
        "x(A)",
        "y(A)",
        "z(A)",
        "x(C)",
        "y(C)",
        "z(C)",
    ],
    "LAI": [
        "length(N)",
        "length(A)",
        "length(C)",
        "angle(N)",
        "angle(A)",
        "angle(C)",
    ],
    "BTI": [
        "tau(NA)",
        "tau(AC)",
        "tau(CN)",
    ],
}


def repeat_number(path: Path) -> int:
    match = re.search(r"repeat_(\d+)", path.name)

    if not match:
        raise ValueError(f"Cannot extract repeat number from {path.name}")

    return int(match.group(1))


def wrap_degrees(values):
    values = np.asarray(values, dtype=float)
    return ((values + 180.0) % 360.0) - 180.0


def circular_mean_degrees(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    radians = np.deg2rad(values)

    mean_sin = np.nanmean(np.sin(radians), axis=0)
    mean_cos = np.nanmean(np.cos(radians), axis=0)

    mean_angle = np.rad2deg(np.arctan2(mean_sin, mean_cos))
    resultant_length = np.sqrt(mean_sin ** 2 + mean_cos ** 2)

    return wrap_degrees(mean_angle), resultant_length


def sequence_match_mask(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(
        {"true", "1", "yes"}
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compare an ensemble of BRI/LAI/BTI invariant CSV files "
            "with an experimental reference."
        )
    )
    parser.add_argument("--prediction-dir", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--target-label", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--include-mismatches",
        action="store_true",
        help="Include sequence mismatches from the alignment mapping.",
    )
    args = parser.parse_args()

    prediction_dir = Path(args.prediction_dir).expanduser().resolve()
    reference_path = Path(args.reference).expanduser().resolve()
    mapping_path = Path(args.mapping).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()

    output_dir.mkdir(parents=True, exist_ok=True)

    prediction_files = sorted(
        prediction_dir.glob("repeat_*.csv"),
        key=repeat_number,
    )

    if not prediction_files:
        raise FileNotFoundError(
            f"No repeat_*.csv files found in {prediction_dir}"
        )

    reference = pd.read_csv(reference_path)
    mapping = pd.read_csv(mapping_path)

    if not args.include_mismatches:
        mapping = mapping[
            sequence_match_mask(mapping["sequence_match"])
        ].copy()

    if mapping.empty:
        raise RuntimeError("No aligned residue pairs remain after filtering.")

    prediction_indices = (
        mapping["prediction_position"].astype(int).to_numpy() - 1
    )
    reference_indices = (
        mapping["reference_position"].astype(int).to_numpy() - 1
    )

    if prediction_indices.min() < 0:
        raise ValueError("Prediction mapping contains an invalid position.")

    if reference_indices.min() < 0:
        raise ValueError("Reference mapping contains an invalid position.")

    frames = []
    expected_columns = None
    expected_labels = None

    for path in prediction_files:
        frame = pd.read_csv(path)

        if expected_columns is None:
            expected_columns = frame.columns.tolist()
            expected_labels = (
                frame["residue_label"].astype(str).tolist()
            )
        else:
            if frame.columns.tolist() != expected_columns:
                raise RuntimeError(
                    f"Column mismatch in prediction file {path}"
                )

            labels = frame["residue_label"].astype(str).tolist()

            if labels != expected_labels:
                raise RuntimeError(
                    f"Sequence mismatch between prediction repeats: {path}"
                )

        if prediction_indices.max() >= len(frame):
            raise IndexError(
                f"Mapping exceeds prediction length for {path}"
            )

        frames.append(frame)

    if reference_indices.max() >= len(reference):
        raise IndexError("Mapping exceeds experimental reference length.")

    required_features = [
        feature
        for features in FEATURE_GROUPS.values()
        for feature in features
    ]

    missing_prediction = [
        feature
        for feature in required_features
        if feature not in expected_columns
    ]
    missing_reference = [
        feature
        for feature in required_features
        if feature not in reference.columns
    ]

    if missing_prediction:
        raise KeyError(
            f"Prediction CSVs lack columns: {missing_prediction}"
        )

    if missing_reference:
        raise KeyError(
            f"Reference CSV lacks columns: {missing_reference}"
        )

    residue_ids = pd.to_numeric(
        mapping["reference_residue_id"],
        errors="coerce",
    ).to_numpy()

    if np.isnan(residue_ids).any():
        residue_ids = (
            mapping["reference_position"].astype(int).to_numpy()
        )

    per_residue_rows = []
    metric_rows = []

    warnings.filterwarnings(
        "ignore",
        message="All-NaN slice encountered",
        category=RuntimeWarning,
    )

    for group_name, features in FEATURE_GROUPS.items():
        if group_name == "BRI":
            rows, columns = 4, 3
            figure_size = (16, 12)
        elif group_name == "LAI":
            rows, columns = 3, 2
            figure_size = (15, 10)
        else:
            rows, columns = 3, 1
            figure_size = (15, 10)

        figure, axes = plt.subplots(
            rows,
            columns,
            figsize=figure_size,
            squeeze=False,
        )

        axes_flat = axes.ravel()

        for feature_index, feature in enumerate(features):
            axis = axes_flat[feature_index]

            prediction_values = np.vstack(
                [
                    frame.iloc[prediction_indices][feature]
                    .to_numpy(dtype=float)
                    for frame in frames
                ]
            )

            reference_values = (
                reference.iloc[reference_indices][feature]
                .to_numpy(dtype=float)
            )

            if group_name != "BTI":
                centre = np.nanmedian(
                    prediction_values,
                    axis=0,
                )
                lower = np.nanpercentile(
                    prediction_values,
                    5,
                    axis=0,
                )
                upper = np.nanpercentile(
                    prediction_values,
                    95,
                    axis=0,
                )
                difference = centre - reference_values

                axis.plot(
                    residue_ids,
                    reference_values,
                    label="Experimental",
                    linewidth=1.2,
                )
                axis.plot(
                    residue_ids,
                    centre,
                    label="Prediction median",
                    linewidth=1.0,
                )
                axis.fill_between(
                    residue_ids,
                    lower,
                    upper,
                    alpha=0.2,
                    label="Prediction 5–95%",
                )

                axis.set_ylabel(feature)

                for index in range(len(mapping)):
                    per_residue_rows.append(
                        {
                            "target": args.target_label,
                            "group": group_name,
                            "feature": feature,
                            "reference_residue_id": residue_ids[index],
                            "prediction_position": (
                                prediction_indices[index] + 1
                            ),
                            "reference_position": (
                                reference_indices[index] + 1
                            ),
                            "reference_value": reference_values[index],
                            "prediction_centre": centre[index],
                            "prediction_p05": lower[index],
                            "prediction_p95": upper[index],
                            "difference": difference[index],
                            "resultant_length": np.nan,
                        }
                    )

            else:
                circular_differences = wrap_degrees(
                    prediction_values
                    - reference_values[np.newaxis, :]
                )

                centre, resultant_length = circular_mean_degrees(
                    circular_differences
                )

                centred_differences = wrap_degrees(
                    circular_differences
                    - centre[np.newaxis, :]
                )

                lower = centre + np.nanpercentile(
                    centred_differences,
                    5,
                    axis=0,
                )
                upper = centre + np.nanpercentile(
                    centred_differences,
                    95,
                    axis=0,
                )

                difference = centre

                axis.axhline(
                    0.0,
                    linewidth=1.2,
                    label="Experimental reference",
                )
                axis.plot(
                    residue_ids,
                    centre,
                    linewidth=1.0,
                    label="Prediction circular mean difference",
                )
                axis.fill_between(
                    residue_ids,
                    lower,
                    upper,
                    alpha=0.2,
                    label="Prediction 5–95%",
                )
                axis.set_ylim(-180, 180)
                axis.set_ylabel("Circular difference (degrees)")

                predicted_angle = wrap_degrees(
                    reference_values + centre
                )

                for index in range(len(mapping)):
                    per_residue_rows.append(
                        {
                            "target": args.target_label,
                            "group": group_name,
                            "feature": feature,
                            "reference_residue_id": residue_ids[index],
                            "prediction_position": (
                                prediction_indices[index] + 1
                            ),
                            "reference_position": (
                                reference_indices[index] + 1
                            ),
                            "reference_value": reference_values[index],
                            "prediction_centre": predicted_angle[index],
                            "prediction_p05": lower[index],
                            "prediction_p95": upper[index],
                            "difference": difference[index],
                            "resultant_length": resultant_length[index],
                        }
                    )

            valid = np.isfinite(difference)

            if valid.any():
                errors = difference[valid]

                metric_rows.append(
                    {
                        "target": args.target_label,
                        "group": group_name,
                        "feature": feature,
                        "mapped_residues": int(valid.sum()),
                        "mae": float(np.mean(np.abs(errors))),
                        "rmse": float(
                            np.sqrt(np.mean(errors ** 2))
                        ),
                        "max_absolute_error": float(
                            np.max(np.abs(errors))
                        ),
                    }
                )

            axis.set_title(feature)
            axis.set_xlabel("Experimental residue ID")
            axis.grid(alpha=0.2)

            if feature_index == 0:
                axis.legend(fontsize=8)

        for unused_axis in axes_flat[len(features):]:
            unused_axis.set_visible(False)

        if group_name == "BTI":
            figure.suptitle(
                f"{args.target_label}: BTI circular differences "
                "from experimental structure"
            )
        else:
            figure.suptitle(
                f"{args.target_label}: {group_name} comparison"
            )

        figure.tight_layout()
        figure.savefig(
            output_dir / f"{group_name}_comparison.png",
            dpi=200,
            bbox_inches="tight",
        )
        plt.close(figure)

    per_residue = pd.DataFrame(per_residue_rows)
    metrics = pd.DataFrame(metric_rows)

    per_residue.to_csv(
        output_dir / "per_residue_invariant_summary.csv",
        index=False,
    )
    metrics.to_csv(
        output_dir / "feature_error_summary.csv",
        index=False,
    )

    manifest = output_dir / "analysis_manifest.txt"

    manifest.write_text(
        "\n".join(
            [
                f"Target: {args.target_label}",
                f"Prediction directory: {prediction_dir}",
                f"Experimental reference: {reference_path}",
                f"Alignment mapping: {mapping_path}",
                f"Prediction CSV count: {len(prediction_files)}",
                f"Mapped residue count used: {len(mapping)}",
                (
                    "Sequence mismatches included: "
                    f"{args.include_mismatches}"
                ),
                (
                    "BTI treatment: signed circular differences "
                    "relative to experimental angles"
                ),
            ]
        )
        + "\n"
    )

    print("Target:", args.target_label)
    print("Prediction CSVs:", len(prediction_files))
    print("Mapped residues used:", len(mapping))
    print("Output directory:", output_dir)
    print()
    print(metrics.to_string(index=False))
    print()
    print("Saved outputs:")

    for path in sorted(output_dir.iterdir()):
        print(path)


if __name__ == "__main__":
    main()
