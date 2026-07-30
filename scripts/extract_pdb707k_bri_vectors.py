#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from bri.invariant_compare import coordinate_value_reshape, extract_chain_info
from bri.pdbx2df import MiniChain, StructureBase


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract compact BRI invariant vectors from a PDB707K manifest."
    )
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--mmcif-root", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--task-index", type=int, required=True)
    parser.add_argument("--batch-size", type=int, required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, keep_default_na=False)
    manifest["pdb_id"] = manifest["pdb_id"].astype(str).str.upper()

    pdb_ids = manifest["pdb_id"].drop_duplicates().tolist()

    start_idx = args.task_index * args.batch_size
    end_idx = min(start_idx + args.batch_size, len(pdb_ids))

    if start_idx >= len(pdb_ids):
        raise ValueError(
            f"Task index {args.task_index} starts beyond "
            f"{len(pdb_ids)} available PDB entries."
        )

    selected_ids = pdb_ids[start_idx:end_idx]
    selected = manifest[manifest["pdb_id"].isin(selected_ids)].copy()

    mmcif_root = Path(args.mmcif_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stem = f"shard_{args.task_index:05d}"
    vectors_path = output_dir / f"{stem}_vectors.npy"
    metadata_path = output_dir / f"{stem}_metadata.csv"
    failures_path = output_dir / f"{stem}_failures.csv"
    summary_path = output_dir / f"{stem}_summary.json"

    vectors: list[np.ndarray] = []
    metadata: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    vector_offset = 0

    def record_failure(row: pd.Series, stage: str, error: object) -> None:
        failures.append(
            {
                "pdb_id": str(row["pdb_id"]),
                "model_id": int(row["model_id"]),
                "chain_id": str(row["chain_id"]),
                "start_residue": int(row["start_residue"]),
                "chain_length": int(row["chain_length"]),
                "stage": stage,
                "error": repr(error),
            }
        )

    with tempfile.TemporaryDirectory(
        prefix=f"bri_extract_{args.task_index:05d}_",
        dir="/tmp",
    ) as tmp:
        tmp_dir = Path(tmp)

        for pdb_id, rows in selected.groupby("pdb_id", sort=False):
            pdb_lower = pdb_id.lower()
            source = mmcif_root / pdb_lower[1:3] / f"{pdb_lower}.cif.gz"
            local_cif = tmp_dir / f"{pdb_lower}.cif"

            if not source.is_file():
                for _, row in rows.iterrows():
                    record_failure(row, "source", FileNotFoundError(source))
                continue

            try:
                with gzip.open(source, "rb") as src, local_cif.open("wb") as dst:
                    shutil.copyfileobj(src, dst)

                structure = StructureBase(str(local_cif))
                coordinates = structure.coordinates

            except Exception as exc:
                for _, row in rows.iterrows():
                    record_failure(row, "parse", exc)
                local_cif.unlink(missing_ok=True)
                continue

            for _, row in rows.iterrows():
                try:
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
                    invariant["pdb_id"] = pdb_id

                    if len(invariant) != chain_length:
                        raise ValueError(
                            f"Expected {chain_length} invariant residues, "
                            f"got {len(invariant)}"
                        )

                    vector = coordinate_value_reshape(
                        invariant,
                        chain_length,
                    )[0].astype(np.float64, copy=False)

                    expected_size = 9 * chain_length
                    if vector.size != expected_size:
                        raise ValueError(
                            f"Expected vector size {expected_size}, "
                            f"got {vector.size}"
                        )

                    if not np.isfinite(vector).all():
                        raise ValueError("Invariant vector contains non-finite values")

                    chain_info = extract_chain_info(invariant.copy())
                    if len(chain_info) != 1:
                        raise ValueError(
                            f"Expected one chain-info row, got {len(chain_info)}"
                        )

                    bri_seq = str(chain_info.iloc[0]["seq"])

                    vectors.append(vector)

                    metadata.append(
                        {
                            "pdb_id": pdb_id,
                            "model_id": model_id,
                            "chain_id": chain_id,
                            "start_residue": start,
                            "chain_length": chain_length,
                            "dataset_seq": str(row["seq"]),
                            "bri_seq": bri_seq,
                            "vector_offset": vector_offset,
                            "vector_size": vector.size,
                        }
                    )

                    vector_offset += vector.size

                except Exception as exc:
                    record_failure(row, "chain", exc)

            local_cif.unlink(missing_ok=True)

    flat_vectors = (
        np.concatenate(vectors).astype(np.float64, copy=False)
        if vectors
        else np.empty(0, dtype=np.float64)
    )

    np.save(vectors_path, flat_vectors)

    metadata_df = pd.DataFrame(metadata)
    failures_df = pd.DataFrame(
        failures,
        columns=[
            "pdb_id",
            "model_id",
            "chain_id",
            "start_residue",
            "chain_length",
            "stage",
            "error",
        ],
    )

    metadata_df.to_csv(metadata_path, index=False)
    failures_df.to_csv(failures_path, index=False)

    accounted = len(metadata_df) + len(failures_df)
    expected = len(selected)

    if accounted != expected:
        raise RuntimeError(
            f"Accounting failure: expected {expected} chains, "
            f"but successes + failures = {accounted}"
        )

    summary = {
        "task_index": args.task_index,
        "batch_size": args.batch_size,
        "pdb_start_index": start_idx,
        "pdb_end_index_exclusive": end_idx,
        "pdb_entries": len(selected_ids),
        "expected_chains": expected,
        "successful_chains": len(metadata_df),
        "failed_chains": len(failures_df),
        "vector_values": int(flat_vectors.size),
    }

    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print(f"Vectors:  {vectors_path}")
    print(f"Metadata: {metadata_path}")
    print(f"Failures: {failures_path}")
    print(f"Summary:  {summary_path}")


if __name__ == "__main__":
    main()
