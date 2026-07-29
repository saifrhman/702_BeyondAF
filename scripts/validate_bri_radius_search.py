from __future__ import annotations

import gzip
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from bri.invariant_compare import (
    coordinate_value_reshape,
    extract_chain_info,
    group_invariant_compare,
)
from bri.pdbx2df import MiniChain, StructureBase

LENGTH = 162
RADIUS = 1.0

SCRATCH = Path(os.environ["COMP702_SCRATCH_ROOT"]) / "bri_geometric"
INPUT_CSV = SCRATCH / "pilot_100" / "length162_chains.csv"
MMCIF_ROOT = SCRATCH / "mmcif_raw"
OUTPUT_DIR = SCRATCH / "pilot_100" / "radius_validation"


def load_invariants() -> tuple[pd.DataFrame, pd.DataFrame]:
    dataset = pd.read_csv(INPUT_CSV, keep_default_na=False)
    invariants: list[pd.DataFrame] = []
    failures: list[dict[str, object]] = []

    with tempfile.TemporaryDirectory(prefix="bri_radius_", dir="/tmp") as tmp:
        tmp_dir = Path(tmp)

        for pdb_id, rows in dataset.groupby("pdb_id", sort=False):
            pdb_lower = str(pdb_id).lower()
            source = MMCIF_ROOT / pdb_lower[1:3] / f"{pdb_lower}.cif.gz"
            local_cif = tmp_dir / f"{pdb_lower}.cif"

            try:
                with gzip.open(source, "rb") as src, local_cif.open("wb") as dst:
                    shutil.copyfileobj(src, dst)

                structure = StructureBase(str(local_cif))
                coordinates = structure.coordinates

                for _, row in rows.iterrows():
                    model_id = int(row["model_id"])
                    chain_id = str(row["chain_id"])
                    start = int(row["start_residue"])
                    chain_length = int(row["chain_length"])

                    segment = coordinates[
                        (coordinates["model_id"] == model_id)
                        & (coordinates["chain_id"] == chain_id)
                        & (coordinates["residue_id"] >= start)
                        & (coordinates["residue_id"] < start + chain_length)
                    ]

                    chain = MiniChain(
                        str(local_cif),
                        model_id=model_id,
                        chain_id=chain_id,
                        start_residue=start,
                        chain_length=chain_length,
                        data=segment,
                    )

                    invariant = chain.invariant.copy()
                    invariant["pdb_id"] = str(pdb_id).upper()
                    invariants.append(invariant)

            except Exception as exc:
                failures.append({"pdb_id": pdb_id, "error": repr(exc)})
            finally:
                local_cif.unlink(missing_ok=True)

    return pd.concat(invariants, ignore_index=True), pd.DataFrame(failures)


def compare_methods(invariants: pd.DataFrame):
    original = group_invariant_compare(
        invariants.copy(),
        seq_compare=True,
    )
    original_lt1 = original[original["distance"] < RADIUS].copy()

    chains = extract_chain_info(invariants.copy())
    vectors = coordinate_value_reshape(invariants, LENGTH)

    pairs = cKDTree(vectors).query_pairs(
        r=RADIUS,
        p=np.inf,
        output_type="ndarray",
    )

    rows = []
    for i, j in pairs:
        distance = float(np.max(np.abs(vectors[i] - vectors[j])))
        if distance < RADIUS:
            a = chains.iloc[i]
            b = chains.iloc[j]
            rows.append(
                {
                    "distance": distance,
                    "pdb_id1": a["pdb_id"],
                    "model_id1": a["model_id"],
                    "chain_id1": a["chain_id"],
                    "pdb_id2": b["pdb_id"],
                    "model_id2": b["model_id"],
                    "chain_id2": b["chain_id"],
                }
            )

    radius_lt1 = pd.DataFrame(rows)
    return original, original_lt1, radius_lt1


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    invariants, failures = load_invariants()
    failures.to_csv(OUTPUT_DIR / "failures.csv", index=False)

    original, original_lt1, radius_lt1 = compare_methods(invariants)

    pair_cols = [
        "pdb_id1", "model_id1", "chain_id1",
        "pdb_id2", "model_id2", "chain_id2",
    ]

    if radius_lt1.empty:
        radius_lt1 = pd.DataFrame(columns=["distance", *pair_cols])

    original_check = original_lt1[pair_cols + ["distance"]].copy()
    radius_check = radius_lt1[pair_cols + ["distance"]].copy()

    merged = original_check.merge(
        radius_check,
        on=pair_cols,
        how="outer",
        suffixes=("_bri", "_radius"),
        indicator=True,
    )

    common = merged[merged["_merge"] == "both"].copy()
    max_diff = (
        float(
            np.max(
                np.abs(
                    common["distance_bri"].to_numpy(float)
                    - common["distance_radius"].to_numpy(float)
                )
            )
        )
        if len(common)
        else 0.0
    )

    original.to_csv(OUTPUT_DIR / "bri_original_all_pairs.csv", index=False)
    original_lt1.to_csv(OUTPUT_DIR / "bri_original_lt1.csv", index=False)
    radius_lt1.to_csv(OUTPUT_DIR / "radius_lt1.csv", index=False)
    merged.to_csv(OUTPUT_DIR / "method_comparison.csv", index=False)

    print("Invariant chains:", invariants["pdb_id"].astype(str).str.cat(
        invariants["chain_id"].astype(str), sep="_"
    ).nunique())
    print("Failures:", len(failures))
    print("BRI original all pairs:", len(original))
    print("BRI original distance < 1:", len(original_lt1))
    print("Radius search distance < 1:", len(radius_lt1))
    print("Pairs only in BRI:", int((merged["_merge"] == "left_only").sum()))
    print("Pairs only in radius search:", int((merged["_merge"] == "right_only").sum()))
    print("Maximum distance difference:", max_diff)


if __name__ == "__main__":
    main()
