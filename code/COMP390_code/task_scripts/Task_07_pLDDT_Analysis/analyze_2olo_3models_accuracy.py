from pathlib import Path
import argparse
import re

import numpy as np
import pandas as pd
from Bio import Align
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1
from Bio.SVDSuperimposer import SVDSuperimposer


MODEL_NAMES = {
    "alphafold_params": "AlphaFold",
    "openfold_official": "OpenFold",
    "trained_model_epoch_32": "Retrained_epoch_32",
}


def get_ca_residues(structure, chain_id=None):
    model = next(structure.get_models())

    if chain_id is not None:
        chain = model[chain_id]
    else:
        chains = list(model.get_chains())
        if len(chains) != 1:
            raise RuntimeError(
                f"Expected one prediction chain, found {len(chains)}"
            )
        chain = chains[0]

    return [
        residue
        for residue in chain
        if residue.id[0] == " " and "CA" in residue
    ]


def residue_sequence(residues):
    return "".join(
        seq1(
            residue.resname,
            custom_map={"MSE": "M"},
            undef_code="X",
        )
        for residue in residues
    )


def matched_indices(prediction_sequence, reference_sequence):
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(
        prediction_sequence,
        reference_sequence,
    )[0]

    prediction_indices, reference_indices = alignment.indices

    matches = []
    prediction_only = []
    reference_only = []
    substitutions = []

    for pred_i, ref_i in zip(prediction_indices, reference_indices):
        pred_i = int(pred_i)
        ref_i = int(ref_i)

        if pred_i >= 0 and ref_i >= 0:
            matches.append((pred_i, ref_i))

            if prediction_sequence[pred_i] != reference_sequence[ref_i]:
                substitutions.append(
                    (
                        pred_i + 1,
                        prediction_sequence[pred_i],
                        ref_i + 1,
                        reference_sequence[ref_i],
                    )
                )
        elif pred_i >= 0:
            prediction_only.append(
                (pred_i + 1, prediction_sequence[pred_i])
            )
        elif ref_i >= 0:
            reference_only.append(
                (ref_i + 1, reference_sequence[ref_i])
            )

    return matches, prediction_only, reference_only, substitutions


def align_coordinates(reference_coordinates, prediction_coordinates):
    superimposer = SVDSuperimposer()
    superimposer.set(reference_coordinates, prediction_coordinates)
    superimposer.run()

    rotation, translation = superimposer.get_rotran()
    aligned_prediction = prediction_coordinates @ rotation + translation

    distances = np.linalg.norm(
        aligned_prediction - reference_coordinates,
        axis=1,
    )

    return float(superimposer.get_rms()), distances


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    root = args.root.resolve()
    reference_path = args.reference.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    reference_structure = MMCIFParser(QUIET=True).get_structure(
        "2olo_reference",
        reference_path,
    )
    reference_residues = get_ca_residues(reference_structure, "A")
    reference_sequence = residue_sequence(reference_residues)

    prediction_files = sorted(
        root.rglob("2olo_A_model_1_unrelaxed.pdb")
    )

    if len(prediction_files) != 15:
        raise RuntimeError(
            f"Expected 15 predictions, found {len(prediction_files)}"
        )

    prediction_rows = []
    residue_rows = []
    alignment_lines = []

    for prediction_path in prediction_files:
        relative_path = prediction_path.relative_to(root)
        model_directory = relative_path.parts[0]

        repeat_match = re.search(
            r"repeat_(\d+)",
            relative_path.as_posix(),
        )
        if not repeat_match:
            raise RuntimeError(
                f"Could not identify repeat from {prediction_path}"
            )

        repeat = int(repeat_match.group(1))
        model_name = MODEL_NAMES[model_directory]

        prediction_structure = PDBParser(QUIET=True).get_structure(
            f"{model_name}_{repeat}",
            prediction_path,
        )
        prediction_residues = get_ca_residues(prediction_structure)
        prediction_sequence = residue_sequence(prediction_residues)

        (
            matches,
            prediction_only,
            reference_only,
            substitutions,
        ) = matched_indices(prediction_sequence, reference_sequence)

        if len(matches) != 388:
            raise RuntimeError(
                f"{prediction_path}: expected 388 matched residues, "
                f"found {len(matches)}"
            )

        reference_coordinates = np.array(
            [
                reference_residues[ref_i]["CA"].coord
                for _, ref_i in matches
            ],
            dtype=float,
        )
        prediction_coordinates = np.array(
            [
                prediction_residues[pred_i]["CA"].coord
                for pred_i, _ in matches
            ],
            dtype=float,
        )

        rmsd, distances = align_coordinates(
            reference_coordinates,
            prediction_coordinates,
        )

        plddt_values = np.array(
            [
                prediction_residues[pred_i]["CA"].bfactor
                for pred_i, _ in matches
            ],
            dtype=float,
        )

        pearson = pd.Series(plddt_values).corr(
            pd.Series(distances),
            method="pearson",
        )
        spearman = pd.Series(plddt_values).corr(
            pd.Series(distances),
            method="spearman",
        )

        prediction_rows.append(
            {
                "model": model_name,
                "repeat": repeat,
                "seed": 1000 + repeat,
                "matched_residues": len(matches),
                "ca_rmsd_angstrom": rmsd,
                "mean_ca_error_angstrom": float(distances.mean()),
                "median_ca_error_angstrom": float(
                    np.median(distances)
                ),
                "max_ca_error_angstrom": float(distances.max()),
                "mean_plddt": float(plddt_values.mean()),
                "pearson_plddt_vs_error": pearson,
                "spearman_plddt_vs_error": spearman,
                "pdb_path": str(prediction_path),
            }
        )

        for pair_index, ((pred_i, ref_i), distance) in enumerate(
            zip(matches, distances),
            start=1,
        ):
            pred_residue = prediction_residues[pred_i]
            ref_residue = reference_residues[ref_i]

            residue_rows.append(
                {
                    "model": model_name,
                    "repeat": repeat,
                    "matched_pair_index": pair_index,
                    "prediction_position": pred_i + 1,
                    "prediction_residue": prediction_sequence[pred_i],
                    "reference_sequence_position": ref_i + 1,
                    "reference_residue_number": ref_residue.id[1],
                    "reference_residue": reference_sequence[ref_i],
                    "plddt": pred_residue["CA"].bfactor,
                    "ca_error_angstrom": float(distance),
                    "pdb_path": str(prediction_path),
                }
            )

        alignment_lines.append(
            f"{model_name}, repeat {repeat}: "
            f"matched={len(matches)}, "
            f"prediction_only={prediction_only}, "
            f"reference_only={reference_only}, "
            f"substitutions={substitutions}"
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
        residue_df.sort_values(
            ["_order", "repeat", "matched_pair_index"]
        )
        .drop(columns="_order")
        .reset_index(drop=True)
    )

    model_summary = (
        prediction_df.groupby("model", sort=False)
        .agg(
            predictions=("ca_rmsd_angstrom", "count"),
            mean_ca_rmsd=("ca_rmsd_angstrom", "mean"),
            std_ca_rmsd=("ca_rmsd_angstrom", "std"),
            min_ca_rmsd=("ca_rmsd_angstrom", "min"),
            max_ca_rmsd=("ca_rmsd_angstrom", "max"),
            mean_ca_error=("mean_ca_error_angstrom", "mean"),
            mean_plddt=("mean_plddt", "mean"),
            mean_pearson_plddt_error=(
                "pearson_plddt_vs_error",
                "mean",
            ),
            mean_spearman_plddt_error=(
                "spearman_plddt_vs_error",
                "mean",
            ),
        )
        .reset_index()
    )

    prediction_df.to_csv(
        out_dir / "2olo_accuracy_per_prediction.csv",
        index=False,
    )
    residue_df.to_csv(
        out_dir / "2olo_accuracy_per_residue.csv",
        index=False,
    )
    model_summary.to_csv(
        out_dir / "2olo_accuracy_model_summary.csv",
        index=False,
    )

    with (out_dir / "2olo_sequence_alignment_summary.txt").open(
        "w"
    ) as handle:
        handle.write("\n".join(alignment_lines) + "\n")

    print(f"Predictions processed: {len(prediction_df)}")
    print(f"Matched residue records: {len(residue_df)}")
    print()
    print(model_summary.to_string(index=False))
    print()
    print(f"Saved outputs to: {out_dir}")


if __name__ == "__main__":
    main()
