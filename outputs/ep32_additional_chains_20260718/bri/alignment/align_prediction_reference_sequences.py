#!/usr/bin/env python3

from pathlib import Path
import argparse

import pandas as pd
from Bio.Align import PairwiseAligner


def sequence_from_frame(frame: pd.DataFrame) -> str:
    labels = frame["residue_label"].astype(str).str.strip().tolist()

    invalid = [
        (index, label)
        for index, label in enumerate(labels)
        if len(label) != 1
    ]

    if invalid:
        raise ValueError(
            f"Invalid one-letter residue labels encountered: {invalid[:10]}"
        )

    return "".join(labels)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prediction", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    prediction_path = Path(args.prediction)
    reference_path = Path(args.reference)
    output_path = Path(args.output)

    prediction = pd.read_csv(prediction_path)
    reference = pd.read_csv(reference_path)

    prediction_sequence = sequence_from_frame(prediction)
    reference_sequence = sequence_from_frame(reference)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(
        prediction_sequence,
        reference_sequence,
    )[0]

    prediction_blocks, reference_blocks = alignment.aligned

    mapping_rows = []

    for prediction_block, reference_block in zip(
        prediction_blocks,
        reference_blocks,
    ):
        prediction_start, prediction_end = map(int, prediction_block)
        reference_start, reference_end = map(int, reference_block)

        prediction_length = prediction_end - prediction_start
        reference_length = reference_end - reference_start

        if prediction_length != reference_length:
            raise RuntimeError(
                "Aligned blocks have inconsistent lengths: "
                f"{prediction_block} versus {reference_block}"
            )

        for offset in range(prediction_length):
            prediction_position = prediction_start + offset
            reference_position = reference_start + offset

            prediction_row = prediction.iloc[prediction_position]
            reference_row = reference.iloc[reference_position]

            prediction_label = str(prediction_row["residue_label"])
            reference_label = str(reference_row["residue_label"])

            mapping_rows.append(
                {
                    "prediction_position": prediction_position + 1,
                    "reference_position": reference_position + 1,
                    "prediction_residue_id": prediction_row["residue_id"],
                    "reference_residue_id": reference_row["residue_id"],
                    "prediction_residue_label": prediction_label,
                    "reference_residue_label": reference_label,
                    "sequence_match": prediction_label == reference_label,
                }
            )

    mapping = pd.DataFrame(mapping_rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(output_path, index=False)

    mapped = len(mapping)
    matches = int(mapping["sequence_match"].sum())

    prediction_coverage = mapped / len(prediction)
    reference_coverage = mapped / len(reference)
    identity = matches / mapped if mapped else 0.0

    print(f"Prediction CSV: {prediction_path}")
    print(f"Reference CSV:  {reference_path}")
    print(f"Prediction length: {len(prediction)}")
    print(f"Reference length:  {len(reference)}")
    print(f"Aligned residue pairs: {mapped}")
    print(f"Exact residue matches: {matches}")
    print(f"Sequence identity over aligned positions: {identity:.4%}")
    print(f"Prediction coverage: {prediction_coverage:.4%}")
    print(f"Reference coverage:  {reference_coverage:.4%}")
    print(f"Alignment score: {alignment.score}")
    print(f"Saved mapping: {output_path}")
    print()
    print(alignment)


if __name__ == "__main__":
    main()
