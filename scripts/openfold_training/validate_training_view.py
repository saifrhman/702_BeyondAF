#!/usr/bin/env python3
"""Prove OpenFold's view of every Gold chain is the retained chain itself.

The invariant this establishes, for every one of the 499,770 training chains:

    OpenFold projected sequence == Gold retained_sequence == MSA query row

Character for character.  The first equality is what this script measures at
full population; the third is guaranteed by construction, because the MSA store
is content-addressed by ``sha256(retained_sequence)`` and
``verify_alignment_index.py`` reads MSAs back through the index to confirm it.

Why it can fail is worth stating, because a length check alone would miss it.
OpenFold indexes the deposited polymer, ``seq_idx = label_seq_id - min(
_entity_poly_seq.num)``.  Gold records ``label_seq_id`` values.  If the
numbering origin were assumed rather than read, every residue would shift
against its own coordinates while all lengths stayed correct, and no shape
check anywhere in the pipeline would notice.

A subsample additionally builds real feature dicts and asserts the structure
and MSA feature widths agree, since that is the failure this whole exercise
exists to remove.

Sharding is by mmCIF entry so each file is parsed exactly once.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import traceback
from collections import defaultdict
from pathlib import Path

import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from pdbclean.openfold_training_view import (
    ProjectionError,
    expand_ranges,
    project_mmcif_object,
    seqres_start_number,
)


def load_chain_rows(path: Path) -> dict[str, list[dict]]:
    """Group Gold chains by mmCIF entry."""

    columns = [
        "openfold_chain_id",
        "pdb_id",
        "auth_chain_id",
        "label_chain_id",
        "entity_id",
        "retained_sequence",
        "retained_residue_count",
        "terminal_trimmed",
    ]

    by_entry: dict[str, list[dict]] = defaultdict(list)

    parquet = pq.ParquetFile(path)

    for batch in parquet.iter_batches(batch_size=20000, columns=columns):
        data = batch.to_pydict()

        for i, chain in enumerate(data["openfold_chain_id"]):
            by_entry[data["pdb_id"][i]].append(
                {
                    "chain": chain,
                    "auth": data["auth_chain_id"][i],
                    "label": data["label_chain_id"][i],
                    "entity": str(data["entity_id"][i]),
                    "sequence": data["retained_sequence"][i],
                    "count": data["retained_residue_count"][i],
                    "trimmed": bool(data["terminal_trimmed"][i]),
                }
            )

    return by_entry


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-chains", required=True, type=Path)
    parser.add_argument("--projection-index", required=True, type=Path)
    parser.add_argument("--mmcif-dir", required=True, type=Path)
    parser.add_argument("--alignment-index", default=None, type=Path)
    parser.add_argument("--msa-store", default=None, type=Path)
    parser.add_argument("--shard-id", required=True, type=int)
    parser.add_argument("--shard-count", required=True, type=int)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument(
        "--feature-checks",
        default=25,
        type=int,
        help="How many chains in this shard to build full feature dicts for.",
    )

    args = parser.parse_args()

    from openfold.data import mmcif_parsing
    from openfold.data.data_pipeline import DataPipeline

    by_entry = load_chain_rows(args.training_chains)

    with args.projection_index.open() as handle:
        projection = json.load(handle)

    alignment_index = None
    if args.alignment_index is not None:
        with args.alignment_index.open() as handle:
            alignment_index = json.load(handle)

    entries = sorted(by_entry)
    mine = [e for i, e in enumerate(entries) if i % args.shard_count == args.shard_id]

    pipeline = DataPipeline(template_featurizer=None)

    report: dict[str, object] = {
        "stage": "openfold_training_view_validation",
        "shard_id": args.shard_id,
        "shard_count": args.shard_count,
        "entries_in_shard": len(mine),
    }

    chains_checked = 0
    sequence_ok = 0
    feature_checks = 0
    feature_ok = 0

    failures: list[dict] = []
    trimmed_checked = 0
    auth_differs_checked = 0

    def fail(chain: str, kind: str, detail: str) -> None:
        failures.append({"chain": chain, "kind": kind, "detail": detail})

    for pdb_id in mine:
        rows = by_entry[pdb_id]
        path = args.mmcif_dir / f"{pdb_id}.cif"

        try:
            with path.open() as handle:
                parsed = mmcif_parsing.parse(
                    file_id=pdb_id, mmcif_string=handle.read()
                )
        except OSError as error:
            for row in rows:
                fail(row["chain"], "mmcif_unreadable", str(error))
            continue

        if parsed.mmcif_object is None:
            detail = str(list(parsed.errors.values())[:1])

            for row in rows:
                fail(row["chain"], "mmcif_unparsable", detail)
            continue

        obj = parsed.mmcif_object

        for row in rows:
            chains_checked += 1
            chain = row["chain"]

            if row["trimmed"]:
                trimmed_checked += 1

            if row["auth"] != row["label"]:
                auth_differs_checked += 1

            entry = projection.get(chain)

            if entry is None:
                fail(chain, "no_projection_entry", "absent from projection index")
                continue

            residues = expand_ranges(entry["r"])

            if len(residues) != row["count"]:
                fail(
                    chain,
                    "lineage_count_mismatch",
                    f"index {len(residues)} vs manifest {row['count']}",
                )
                continue

            try:
                start = seqres_start_number(obj, entry["e"])
                projected = project_mmcif_object(
                    obj,
                    row["auth"],
                    residues,
                    seq_start_num=start,
                )
            except ProjectionError as error:
                fail(chain, "projection_error", str(error))
                continue
            except Exception as error:  # noqa: BLE001 - reported, not swallowed
                fail(chain, "projection_crash", f"{type(error).__name__}: {error}")
                continue

            projected_sequence = projected.chain_to_seqres[row["auth"]]

            if projected_sequence != row["sequence"]:
                fail(
                    chain,
                    "sequence_mismatch",
                    f"projected len {len(projected_sequence)} vs retained "
                    f"len {len(row['sequence'])}; "
                    f"projected[:40]={projected_sequence[:40]!r} "
                    f"retained[:40]={row['sequence'][:40]!r}",
                )
                continue

            # Residue mapping must be dense and in order: get_atom_coords walks
            # range(num_res) and would KeyError or silently shift otherwise.
            keys = sorted(projected.seqres_to_structure[row["auth"]])

            if keys != list(range(len(projected_sequence))):
                fail(chain, "index_not_dense", f"keys[:5]={keys[:5]}")
                continue

            sequence_ok += 1

            if feature_checks < args.feature_checks and alignment_index:
                align = alignment_index.get(chain)

                if align is None:
                    fail(chain, "no_alignment_entry", "absent from alignment index")
                    continue

                feature_checks += 1

                try:
                    feats = pipeline.process_mmcif(
                        mmcif=projected,
                        alignment_dir=str(args.msa_store),
                        chain_id=row["auth"],
                        alignment_index=align,
                    )
                except Exception as error:  # noqa: BLE001
                    fail(
                        chain,
                        "feature_build_failed",
                        f"{type(error).__name__}: {error}",
                    )
                    continue

                width = feats["aatype"].shape[0]
                msa_width = feats["msa"].shape[1]
                positions = feats["all_atom_positions"].shape[0]
                mask = feats["all_atom_mask"].shape[0]

                if not (width == msa_width == positions == mask == row["count"]):
                    fail(
                        chain,
                        "feature_width_disagreement",
                        f"aatype={width} msa={msa_width} pos={positions} "
                        f"mask={mask} retained={row['count']}",
                    )
                    continue

                feature_ok += 1

    report.update(
        {
            "chains_checked": chains_checked,
            "sequence_matches": sequence_ok,
            "terminal_trimmed_checked": trimmed_checked,
            "auth_differs_from_label_checked": auth_differs_checked,
            "feature_checks": feature_checks,
            "feature_matches": feature_ok,
            "failure_count": len(failures),
            "failures": failures[:50],
            "status": "FAIL" if failures else "PASS",
        }
    )

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print(f"shard {args.shard_id}/{args.shard_count}")
    print(f"  entries          : {len(mine):,}")
    print(f"  chains checked   : {chains_checked:,}")
    print(f"  sequence matches : {sequence_ok:,}")
    print(f"  trimmed covered  : {trimmed_checked:,}")
    print(f"  auth!=label      : {auth_differs_checked:,}")
    print(f"  feature checks   : {feature_ok:,}/{feature_checks:,}")
    print(f"  failures         : {len(failures):,}")

    for failure in failures[:10]:
        print(f"    {failure['chain']}: {failure['kind']}: {failure['detail'][:160]}",
              file=sys.stderr)

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
