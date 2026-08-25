#!/usr/bin/env python3
"""Record each Gold chain's retained lineage compactly enough to load per worker.

The training view needs, for every chain, the mmCIF ``label_seq_id`` values
Protocol 3.2 retained.  Stored literally that is ~125 million integers, which
as JSON is hundreds of megabytes and is paid again in every dataloader worker.

Terminal trimming, though, removes residues from the ends, so a chain's lineage
is usually one contiguous run.  Storing runs instead of residues makes the
common case two integers.  Whether that assumption actually holds across the
corpus is not assumed here -- it is measured, and the run form is exact either
way, so a chain with internal excisions simply costs more runs.

Output is one JSON object per chain::

    {"101m_A": {"e": "1", "r": [[1, 154]], "s": "3f2a...16hex"}}

where ``e`` is the mmCIF entity id, needed to establish the seqres numbering
origin, ``r`` is the inclusive run list, and ``s`` is a truncated SHA256 of the
retained sequence.

``s`` is what makes the projection self-checking at training time.  A length
check is not sufficient: entries with point microheterogeneity carry two
``_entity_poly_seq`` rows for one residue number, OpenFold appends both, and
every index past that point shifts by one -- producing a projected sequence of
exactly the right *length* and the wrong *content*.  Comparing the digest
catches that, so a chain can never be trained on against an MSA built for a
different sequence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from pdbclean.openfold_training_view import contiguous_ranges


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)

    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-chains", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--manifest", default=None, type=Path)
    parser.add_argument("--expected-chains", default=None, type=int)
    parser.add_argument("--batch-rows", default=20000, type=int)

    args = parser.parse_args()

    columns = [
        "openfold_chain_id",
        "entity_id",
        "auth_chain_id",
        "retained_label_seq_ids",
        "retained_residue_count",
        "retained_sequence",
        "retained_start_label_seq_id",
        "retained_end_label_seq_id",
        "terminal_trimmed",
    ]

    index: dict[str, dict] = {}

    single_run = 0
    multi_run = 0
    max_runs = 0
    total_runs = 0
    count_disagreements: list[str] = []
    unsorted_lineage: list[str] = []

    parquet = pq.ParquetFile(args.training_chains)

    for batch in parquet.iter_batches(
        batch_size=args.batch_rows, columns=columns
    ):
        data = batch.to_pydict()

        for i, chain in enumerate(data["openfold_chain_id"]):
            residues = data["retained_label_seq_ids"][i]

            if list(residues) != sorted(residues):
                unsorted_lineage.append(chain)

            runs = contiguous_ranges(residues)

            total_runs += len(runs)
            max_runs = max(max_runs, len(runs))

            if len(runs) == 1:
                single_run += 1
            else:
                multi_run += 1

            # The manifest carries a residue count; the lineage must agree.
            if len(residues) != data["retained_residue_count"][i]:
                count_disagreements.append(chain)

            digest = hashlib.sha256(
                data["retained_sequence"][i].encode()
            ).hexdigest()[:16]

            index[chain] = {
                "e": str(data["entity_id"][i]),
                "r": [[int(a), int(b)] for a, b in runs],
                "s": digest,
            }

    print(f"chains                    : {len(index):,}")
    print(f"single-run lineages       : {single_run:,}")
    print(f"multi-run lineages        : {multi_run:,}")
    print(f"max runs for one chain    : {max_runs}")
    print(f"total runs                : {total_runs:,}")
    print(f"count disagreements       : {len(count_disagreements):,}")
    print(f"unsorted lineages         : {len(unsorted_lineage):,}")

    if count_disagreements:
        print(
            f"ERROR: {len(count_disagreements)} chains whose lineage length "
            "disagrees with retained_residue_count",
            file=sys.stderr,
        )

        for chain in count_disagreements[:10]:
            print(f"  {chain}", file=sys.stderr)

        return 1

    if unsorted_lineage:
        print(
            f"ERROR: {len(unsorted_lineage)} chains have an unsorted lineage; "
            "run collapsing assumes ascending residue ids",
            file=sys.stderr,
        )

        for chain in unsorted_lineage[:10]:
            print(f"  {chain}", file=sys.stderr)

        return 1

    if args.expected_chains is not None and len(index) != args.expected_chains:
        print(
            f"ERROR: indexed {len(index)} chains, expected "
            f"{args.expected_chains}",
            file=sys.stderr,
        )
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w") as handle:
        json.dump(index, handle, separators=(",", ":"), sort_keys=True)

    size = args.output.stat().st_size

    print(f"index written             : {args.output}")
    print(f"index size                : {size:,} bytes")

    if args.manifest is not None:
        manifest = {
            "stage": "openfold_projection_index",
            "stage_version": "1.0",
            "training_chains": str(args.training_chains.resolve()),
            "training_chains_sha256": sha256_of(args.training_chains),
            "projection_index": str(args.output.resolve()),
            "projection_index_sha256": sha256_of(args.output),
            "projection_index_bytes": size,
            "chains": len(index),
            "single_run_lineages": single_run,
            "multi_run_lineages": multi_run,
            "max_runs_for_one_chain": max_runs,
            "total_runs": total_runs,
            "lineage_source": "Gold retained_label_seq_ids "
                              "(openfold_training_chains.parquet)",
        }

        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )

        print(f"manifest written          : {args.manifest}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
