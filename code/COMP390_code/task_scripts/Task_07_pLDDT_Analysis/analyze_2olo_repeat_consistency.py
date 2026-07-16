from pathlib import Path
import argparse
import itertools
import re

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer


MODEL_NAMES = {
    "alphafold_params": "AlphaFold",
    "openfold_official": "OpenFold",
    "trained_model_epoch_32": "Retrained_epoch_32",
}


def get_ca_coordinates(pdb_path):
    structure = PDBParser(QUIET=True).get_structure("prediction", pdb_path)
    model = next(structure.get_models())
    chain = next(model.get_chains())

    residues = [
        residue for residue in chain
        if residue.id[0] == " " and "CA" in residue
    ]

    return np.array(
        [residue["CA"].coord for residue in residues],
        dtype=float,
    )


def aligned_rmsd(reference, mobile):
    sup = SVDSuperimposer()
    sup.set(reference, mobile)
    sup.run()
    return float(sup.get_rms())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    rows = []

    for model_dir, model_name in MODEL_NAMES.items():
        files = sorted(
            (args.root / model_dir).rglob(
                "2olo_A_model_1_unrelaxed.pdb"
            )
        )

        if len(files) != 5:
            raise RuntimeError(
                f"{model_name}: expected 5 files, found {len(files)}"
            )

        structures = {}

        for path in files:
            match = re.search(r"repeat_(\d+)", path.as_posix())
            repeat = int(match.group(1))
            coordinates = get_ca_coordinates(path)

            if len(coordinates) != 389:
                raise RuntimeError(
                    f"{path}: expected 389 residues, "
                    f"found {len(coordinates)}"
                )

            structures[repeat] = coordinates

        for repeat_1, repeat_2 in itertools.combinations(
            sorted(structures), 2
        ):
            rmsd = aligned_rmsd(
                structures[repeat_1],
                structures[repeat_2],
            )

            rows.append(
                {
                    "model": model_name,
                    "repeat_1": repeat_1,
                    "repeat_2": repeat_2,
                    "matched_residues": 389,
                    "pairwise_ca_rmsd_angstrom": rmsd,
                }
            )

    pairwise_df = pd.DataFrame(rows)

    summary_df = (
        pairwise_df.groupby("model", sort=False)
        .agg(
            comparisons=("pairwise_ca_rmsd_angstrom", "count"),
            mean_pairwise_rmsd=(
                "pairwise_ca_rmsd_angstrom",
                "mean",
            ),
            std_pairwise_rmsd=(
                "pairwise_ca_rmsd_angstrom",
                "std",
            ),
            min_pairwise_rmsd=(
                "pairwise_ca_rmsd_angstrom",
                "min",
            ),
            max_pairwise_rmsd=(
                "pairwise_ca_rmsd_angstrom",
                "max",
            ),
        )
        .reset_index()
    )

    pairwise_df.to_csv(
        args.out_dir / "2olo_repeat_pairwise_rmsd.csv",
        index=False,
    )
    summary_df.to_csv(
        args.out_dir / "2olo_repeat_consistency_summary.csv",
        index=False,
    )

    print(f"Pairwise comparisons: {len(pairwise_df)}")
    print()
    print(summary_df.to_string(index=False))
    print()
    print(f"Saved outputs to: {args.out_dir}")


if __name__ == "__main__":
    main()
