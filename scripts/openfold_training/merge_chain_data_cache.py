#!/usr/bin/env python3
"""Merge sharded chain-cache fragments and prove the population is exact.

The cache decides which chains OpenFold trains on: a chain absent from it is
removed by ``OpenFoldSingleDataset.__init__`` with a log line and nothing else,
so an incomplete merge silently shrinks the training set.  The expected
population is therefore asserted rather than reported.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pyarrow.parquet as pq


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)

    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards", required=True, type=Path)
    parser.add_argument("--shard-count", required=True, type=int)
    parser.add_argument("--training-chains", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--manifest", default=None, type=Path)

    args = parser.parse_args()

    failures: list[str] = []

    paths = sorted(args.shards.glob("shard_*.json"))

    if len(paths) != args.shard_count:
        failures.append(
            f"found {len(paths)} shard files, expected {args.shard_count}"
        )

    cache: dict[str, dict] = {}
    collisions: list[str] = []

    for path in paths:
        with path.open() as handle:
            fragment = json.load(handle)

        for chain, entry in fragment.items():
            if chain in cache:
                collisions.append(chain)

            cache[chain] = entry

    if collisions:
        failures.append(
            f"{len(collisions)} chains appeared in more than one shard "
            f"(e.g. {collisions[:5]}); sharding is not disjoint"
        )

    expected = set(
        pq.read_table(args.training_chains, columns=["openfold_chain_id"])
        .column("openfold_chain_id")
        .to_pylist()
    )

    missing = sorted(expected - set(cache))
    extra = sorted(set(cache) - expected)

    if missing:
        failures.append(
            f"{len(missing)} retained chains are absent from the cache and "
            f"would be dropped from training (e.g. {missing[:5]})"
        )

    if extra:
        failures.append(
            f"{len(extra)} cache entries are not retained Gold chains "
            f"(e.g. {extra[:5]})"
        )

    # Every entry must carry what the sampler reads, or training dies mid-run.
    malformed = [
        chain
        for chain, entry in cache.items()
        if "seq" not in entry or entry.get("resolution") is None
    ]

    if malformed:
        failures.append(
            f"{len(malformed)} entries lack 'seq' or 'resolution' "
            f"(e.g. {malformed[:5]})"
        )

    print(f"shard files               : {len(paths)}")
    print(f"chains in cache           : {len(cache):,}")
    print(f"retained chains expected  : {len(expected):,}")
    print(f"missing                   : {len(missing):,}")
    print(f"unexpected                : {len(extra):,}")
    print(f"malformed                 : {len(malformed):,}")

    if failures:
        print()
        print("STATUS: FAIL")

        for failure in failures:
            print(f"  - {failure}")
            print(f"  - {failure}", file=sys.stderr)

        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w") as handle:
        json.dump(cache, handle, separators=(",", ":"), sort_keys=True)

    size = args.output.stat().st_size

    print(f"cache written             : {args.output}")
    print(f"cache size                : {size:,} bytes")

    if args.manifest is not None:
        manifest = {
            "stage": "openfold_train_chain_data_cache",
            "stage_version": "1.0",
            "chain_data_cache": str(args.output.resolve()),
            "chain_data_cache_sha256": sha256_of(args.output),
            "chain_data_cache_bytes": size,
            "chains": len(cache),
            "training_chains": str(args.training_chains.resolve()),
            "sequence_semantics": "Gold retained_sequence (the training view), "
                                  "not the deposited mmCIF seqres",
            "fields": sorted({k for entry in cache.values() for k in entry}),
        }

        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )

        print(f"manifest written          : {args.manifest}")

    print()
    print("STATUS: PASS")
    print("  the cache describes exactly the retained Gold population,")
    print("  with the retained sequence as 'seq'.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
