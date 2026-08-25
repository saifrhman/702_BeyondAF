#!/usr/bin/env python3
"""Map every training chain onto the MSA built for its sequence.

OpenFold's training dataset addresses alignments by chain (``101m_A``), but the
MSA corpus is content-addressed by sequence SHA256, because 499,770 chains
collapse to 142,056 distinct sequences and searching the duplicates would have
been 3.5x the compute for identical answers.  Something has to bridge the two,
and that something is OpenFold's *alignment index*.

The index format is fixed by ``DataPipeline._parse_msa_data``::

    {chain_name: {"db": <path relative to alignment_dir>,
                  "files": [[name, start, size], ...]}}

read as ``open(join(alignment_dir, db))``, ``seek(start)``, ``read(size)``.

Upstream builds that index over a *packed* database: every alignment file is
concatenated into a handful of huge ``.db`` files and the index records byte
offsets into them.  Nothing in the format requires packing, though.  Pointing
``db`` at an individual ``<sha256>.a3m`` with ``start=0`` and ``size`` equal to
the file length satisfies the same contract exactly, and costs nothing: the
corpus is 161 GiB across 142,056 files that already exist, so packing would
duplicate 161 GiB of a 500 GiB quota and buy nothing this pipeline needs.  The
equivalence is not assumed: ``tests/pdbclean/test_alignment_index.py`` asserts
it at the byte level, which is the level OpenFold reads at, and it was
confirmed once against OpenFold's own ``DataPipeline._parse_msa_data`` by
parsing the same MSA both ways and comparing sequences, descriptions and
deletion matrices.

Sequence reuse costs nothing either.  Chains sharing a sequence share one index
entry, so the 357,714 duplicate chains add index keys but not one byte of MSA.

The alignment filename is free-form: the training read path dispatches on the
``.a3m`` extension alone (specific names like ``uniref90_hits.a3m`` are used
only by OpenFold's own alignment *generation* code, which this corpus replaces).
It is named for what it actually is, so a reader of the index can tell.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pyarrow.parquet as pq


#: What the alignment is, recorded in the index so provenance survives.
ALIGNMENT_FILENAME = "uniref30_mmseqs_hits.a3m"


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk)

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


def build_index(
    chains: list[str],
    sequence_shas: list[str],
    store: Path,
) -> tuple[dict[str, dict], list[str], list[str]]:
    """Return (index, missing, empty).

    One entry object is created per *sequence* and shared by every chain that
    uses it.  JSON has no references so the written file repeats them, but the
    build itself stays proportional to the 142,056 sequences rather than the
    499,770 chains.
    """

    entries: dict[str, dict] = {}
    missing: list[str] = []
    empty: list[str] = []

    index: dict[str, dict] = {}

    for chain, sha in zip(chains, sequence_shas):
        entry = entries.get(sha)

        if entry is None:
            path = store / f"{sha}.a3m"

            try:
                size = path.stat().st_size
            except FileNotFoundError:
                missing.append(sha)
                entries[sha] = {}
                continue

            if size == 0:
                empty.append(sha)

            entry = {
                "db": path.name,
                "files": [[ALIGNMENT_FILENAME, 0, size]],
            }
            entries[sha] = entry

        if entry:
            index[chain] = entry

    return index, missing, empty


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--chain-map",
        required=True,
        type=Path,
        help="chain_to_sequence.parquet: openfold_chain_id -> sequence_sha256",
    )
    parser.add_argument(
        "--msa-store",
        required=True,
        type=Path,
        help="Directory of <sequence_sha256>.a3m; becomes OpenFold's "
             "--train_alignment_dir",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--manifest", default=None, type=Path)
    parser.add_argument(
        "--expected-chains",
        default=None,
        type=int,
        help="Fail unless exactly this many chains are indexed.",
    )
    parser.add_argument(
        "--expected-sequences",
        default=None,
        type=int,
        help="Fail unless exactly this many distinct sequences are referenced.",
    )

    args = parser.parse_args()

    # -- inputs ------------------------------------------------------------
    table = pq.read_table(
        args.chain_map, columns=["openfold_chain_id", "sequence_sha256"]
    )

    chains = table.column("openfold_chain_id").to_pylist()
    sequence_shas = table.column("sequence_sha256").to_pylist()

    print(f"chains in mapping         : {len(chains):,}")
    print(f"distinct sequences        : {len(set(sequence_shas)):,}")

    duplicate_chains = len(chains) - len(set(chains))

    if duplicate_chains:
        print(
            f"ERROR: {duplicate_chains} duplicate chain ids in the mapping",
            file=sys.stderr,
        )
        return 1

    # -- build -------------------------------------------------------------
    index, missing, empty = build_index(chains, sequence_shas, args.msa_store)

    print(f"chains indexed            : {len(index):,}")

    if missing:
        print(f"ERROR: {len(missing):,} sequences have no MSA", file=sys.stderr)

        for sha in missing[:10]:
            print(f"  missing: {sha}.a3m", file=sys.stderr)

        return 1

    if empty:
        print(f"ERROR: {len(empty):,} MSAs are empty", file=sys.stderr)

        for sha in empty[:10]:
            print(f"  empty: {sha}.a3m", file=sys.stderr)

        return 1

    referenced = {entry["db"] for entry in index.values()}
    indexed_bytes = sum(
        entry["files"][0][2]
        for entry in {id(e): e for e in index.values()}.values()
    )

    print(f"distinct MSAs referenced  : {len(referenced):,}")
    print(f"MSA bytes addressed       : {indexed_bytes:,}")

    # -- expectation gates --------------------------------------------------
    if args.expected_chains is not None and len(index) != args.expected_chains:
        print(
            f"ERROR: indexed {len(index)} chains, expected "
            f"{args.expected_chains}",
            file=sys.stderr,
        )
        return 1

    if (
        args.expected_sequences is not None
        and len(referenced) != args.expected_sequences
    ):
        print(
            f"ERROR: referenced {len(referenced)} sequences, expected "
            f"{args.expected_sequences}",
            file=sys.stderr,
        )
        return 1

    # -- write --------------------------------------------------------------
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Compact: OpenFold json.loads it, and indentation would roughly triple a
    # file this size for no reader's benefit.
    with args.output.open("w") as handle:
        json.dump(index, handle, separators=(",", ":"), sort_keys=True)

    print(f"index written             : {args.output}")
    print(f"index size                : {args.output.stat().st_size:,} bytes")

    if args.manifest is not None:
        manifest = {
            "stage": "openfold_alignment_index",
            "stage_version": "1.0",
            "alignment_dir": str(args.msa_store.resolve()),
            "alignment_filename": ALIGNMENT_FILENAME,
            "alignment_index": str(args.output.resolve()),
            "alignment_index_sha256": sha256_of(args.output),
            "chain_map": str(args.chain_map.resolve()),
            "chain_map_sha256": sha256_of(args.chain_map),
            "chains_indexed": len(index),
            "distinct_msas_referenced": len(referenced),
            "chains_reusing_a_shared_msa": len(index) - len(referenced),
            "msa_bytes_addressed": indexed_bytes,
            "storage_model": "zero-copy: each entry addresses its own .a3m "
                             "in full (start=0), so no bytes are duplicated",
        }

        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )

        print(f"manifest written          : {args.manifest}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
