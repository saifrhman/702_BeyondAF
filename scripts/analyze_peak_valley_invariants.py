#!/usr/bin/env python3
"""Residue-wise peak-valley analysis across repeated BRI invariant CSV files."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


METADATA_COLUMNS = {
    "model_id",
    "chain_id",
    "residue_id",
    "residue_label",
    "chain_length",
}

TAU_COLUMNS = {"tau(NA)", "tau(AC)", "tau(CN)"}
KNOWN_CANDIDATES = {38, 39, 43, 268, 386}
REPEAT_PATTERN = re.compile(r"repeat[_-]?(\d+)", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate residue-wise minima, maxima, ranges, standard "
            "deviations and peak/valley repeat files from invariant CSVs."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing repeated invariant CSV files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for analysis outputs.",
    )
    parser.add_argument(
        "--expected-files",
        type=int,
        default=None,
        help="Fail unless this number of CSV files is present.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top candidates retained per descriptor.",
    )
    return parser.parse_args()


def extract_repeat(filename: str) -> int | None:
    match = REPEAT_PATTERN.search(filename)
    return int(match.group(1)) if match else None


def structure_stem(filename: str) -> str:
    if filename.endswith("_inv.csv"):
        return filename[: -len("_inv.csv")]
    return Path(filename).stem


def circular_difference_degrees(value_a: float, value_b: float) -> float:
    """Smallest absolute angular difference in degrees."""
    return float(abs(((value_b - value_a + 180.0) % 360.0) - 180.0))


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not input_dir.is_dir():
        print(f"ERROR: input directory does not exist: {input_dir}", file=sys.stderr)
        return 1

    csv_files = sorted(input_dir.glob("*.csv"))

    if not csv_files:
        print(f"ERROR: no CSV files found in {input_dir}", file=sys.stderr)
        return 1

    if args.expected_files is not None and len(csv_files) != args.expected_files:
        print(
            f"ERROR: expected {args.expected_files} CSV files, "
            f"but found {len(csv_files)}.",
            file=sys.stderr,
        )
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)

    long_tables: list[pd.DataFrame] = []
    manifest_records: list[dict[str, object]] = []

    for csv_file in csv_files:
        dataframe = pd.read_csv(csv_file)

        required = {"residue_id", "residue_label"}
        missing_required = required - set(dataframe.columns)
        if missing_required:
            print(
                f"ERROR: {csv_file.name} is missing columns: "
                f"{sorted(missing_required)}",
                file=sys.stderr,
            )
            return 1

        if dataframe["residue_id"].duplicated().any():
            print(
                f"ERROR: duplicate residue IDs in {csv_file.name}",
                file=sys.stderr,
            )
            return 1

        descriptor_columns: list[str] = []

        for column in dataframe.columns:
            if column in METADATA_COLUMNS:
                continue

            converted = pd.to_numeric(dataframe[column], errors="coerce")
            if converted.notna().any():
                dataframe[column] = converted
                descriptor_columns.append(column)

        if not descriptor_columns:
            print(
                f"ERROR: no numeric descriptor columns in {csv_file.name}",
                file=sys.stderr,
            )
            return 1

        repeat_number = extract_repeat(csv_file.name)

        long_table = dataframe[
            ["residue_id", "residue_label"] + descriptor_columns
        ].melt(
            id_vars=["residue_id", "residue_label"],
            value_vars=descriptor_columns,
            var_name="descriptor",
            value_name="value",
        )

        long_table["value"] = pd.to_numeric(
            long_table["value"], errors="coerce"
        )
        long_table = long_table.dropna(subset=["value"])
        long_table["invariant_file"] = csv_file.name
        long_table["repeat"] = repeat_number

        long_tables.append(long_table)

        manifest_records.append(
            {
                "invariant_file": csv_file.name,
                "repeat": repeat_number,
                "rows": len(dataframe),
                "minimum_residue_id": int(dataframe["residue_id"].min()),
                "maximum_residue_id": int(dataframe["residue_id"].max()),
                "descriptor_count": len(descriptor_columns),
            }
        )

    combined = pd.concat(long_tables, ignore_index=True)

    # Ensure each residue number has one consistent amino-acid label.
    label_counts = combined.groupby("residue_id")["residue_label"].nunique()
    inconsistent = label_counts[label_counts > 1]

    if not inconsistent.empty:
        print(
            "ERROR: inconsistent residue labels across repeats for residues: "
            + ", ".join(map(str, inconsistent.index.tolist())),
            file=sys.stderr,
        )
        return 1

    result_records: list[dict[str, object]] = []

    grouped = combined.groupby(
        ["descriptor", "residue_id", "residue_label"],
        sort=True,
        dropna=False,
    )

    for (descriptor, residue_id, residue_label), group in grouped:
        group = group.sort_values(
            ["value", "invariant_file"],
            kind="stable",
        ).reset_index(drop=True)

        values = group["value"].to_numpy(dtype=float)
        valley_row = group.iloc[0]
        peak_row = group.iloc[-1]

        valley_value = float(valley_row["value"])
        peak_value = float(peak_row["value"])
        raw_range = peak_value - valley_value

        is_tau = descriptor in TAU_COLUMNS
        circular_difference = (
            circular_difference_degrees(valley_value, peak_value)
            if is_tau
            else np.nan
        )

        ranking_difference = (
            circular_difference if is_tau else raw_range
        )

        result_records.append(
            {
                "descriptor": descriptor,
                "residue_id": int(residue_id),
                "residue_label": residue_label,
                "number_of_values": len(values),
                "valley_value": valley_value,
                "peak_value": peak_value,
                "raw_range": raw_range,
                "circular_difference_deg": circular_difference,
                "ranking_difference": ranking_difference,
                "standard_deviation_sample": (
                    float(np.std(values, ddof=1))
                    if len(values) > 1
                    else np.nan
                ),
                "standard_deviation_population": float(
                    np.std(values, ddof=0)
                ),
                "mean": float(np.mean(values)),
                "median": float(np.median(values)),
                "valley_repeat": valley_row["repeat"],
                "peak_repeat": peak_row["repeat"],
                "valley_invariant_file": valley_row["invariant_file"],
                "peak_invariant_file": peak_row["invariant_file"],
                "valley_structure_stem": structure_stem(
                    valley_row["invariant_file"]
                ),
                "peak_structure_stem": structure_stem(
                    peak_row["invariant_file"]
                ),
                "is_tau_descriptor": is_tau,
                "possible_wraparound": bool(is_tau and raw_range > 180.0),
                "wraparound_reduction_deg": (
                    raw_range - circular_difference
                    if is_tau
                    else np.nan
                ),
            }
        )

    results = pd.DataFrame(result_records)

    results["rank_within_descriptor"] = (
        results.groupby("descriptor")["ranking_difference"]
        .rank(method="min", ascending=False)
        .astype("Int64")
    )

    results = results.sort_values(
        ["descriptor", "rank_within_descriptor", "residue_id"]
    ).reset_index(drop=True)

    tau_results = results[
        results["is_tau_descriptor"]
    ].sort_values(
        ["ranking_difference", "raw_range"],
        ascending=[False, False],
    )

    top_candidates = results[
        results["rank_within_descriptor"] <= args.top_n
    ].copy()

    known_candidates = results[
        results["residue_id"].isin(KNOWN_CANDIDATES)
    ].sort_values(
        ["descriptor", "residue_id"]
    )

    manifest = pd.DataFrame(manifest_records).sort_values(
        ["repeat", "invariant_file"],
        na_position="last",
    )

    manifest.to_csv(output_dir / "input_manifest.csv", index=False)
    results.to_csv(
        output_dir / "peak_valley_all_descriptors.csv",
        index=False,
    )
    tau_results.to_csv(
        output_dir / "tau_peak_valley_ranked.csv",
        index=False,
    )
    top_candidates.to_csv(
        output_dir / "top_candidates_per_descriptor.csv",
        index=False,
    )
    known_candidates.to_csv(
        output_dir / "known_candidate_residues.csv",
        index=False,
    )

    summary_lines = [
        "2OLO peak-valley analysis",
        "==========================",
        f"Input directory: {input_dir}",
        f"CSV files analysed: {len(csv_files)}",
        f"Residues represented: {combined['residue_id'].nunique()}",
        f"Descriptors represented: {combined['descriptor'].nunique()}",
        "",
        "Tau descriptors are ranked using circular angular difference.",
        "Other descriptors are ranked using their raw maximum-minus-minimum range.",
        "",
    ]

    for descriptor in sorted(TAU_COLUMNS):
        subset = tau_results[
            tau_results["descriptor"] == descriptor
        ].head(10)

        summary_lines.append(f"Top 10 for {descriptor}")
        summary_lines.append("-" * (11 + len(descriptor)))

        for row in subset.itertuples(index=False):
            summary_lines.append(
                f"Residue {row.residue_id} {row.residue_label}: "
                f"valley={row.valley_value:.2f} "
                f"(repeat {row.valley_repeat}), "
                f"peak={row.peak_value:.2f} "
                f"(repeat {row.peak_repeat}), "
                f"raw={row.raw_range:.2f}, "
                f"circular={row.circular_difference_deg:.2f}, "
                f"std={row.standard_deviation_sample:.2f}, "
                f"wraparound={row.possible_wraparound}"
            )

        summary_lines.append("")

    summary_path = output_dir / "summary.txt"
    summary_path.write_text("\n".join(summary_lines), encoding="utf-8")

    print("\n".join(summary_lines))
    print(f"Full outputs written to: {output_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
