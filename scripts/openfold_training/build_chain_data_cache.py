#!/usr/bin/env python3
"""Build OpenFold's train chain-data cache over the retained training view.

OpenFold's sampler is not optional machinery: ``OpenFoldDataset.looped_samples``
indexes ``chain_data_cache[chain_id]`` directly, so a missing cache is a
TypeError on the first batch and a missing *chain* is silently dropped from
training by the filter in ``OpenFoldSingleDataset.__init__``.

Upstream's ``scripts/generate_chain_data_cache.py`` walks a directory of mmCIFs
and records ``mmcif.chain_to_seqres`` for every chain it finds.  That is the
deposited polymer, and using it here would describe a different population than
the one being trained on in two separate ways:

*wrong residues*
    the sampler's length-based sampling probability and its
    ``max_single_aa_prop`` filter would be computed over residues Protocol 3.2
    removed;

*wrong chains*
    every chain of every entry would be described, including chains Gold did
    not retain.

Both are corrected here by driving the cache from the Gold manifest rather than
from the directory listing: the population is exactly the retained chains, and
``seq`` is the retained sequence -- the same string the MSA was searched on and
the same one the projected mmCIF yields.

Only ``resolution`` and ``release_date`` come from the deposited file, because
they are properties of the experiment rather than of the chain, and they are
read with OpenFold's own ``_get_header`` so the values match what OpenFold
would have computed.  The structure itself is never built, which makes this
much cheaper than a full parse.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import pyarrow.parquet as pq


def header_for(path: Path) -> dict:
    """Resolution and release date, exactly as OpenFold computes them."""

    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from openfold.data.mmcif_parsing import _get_header

    parsed_info = MMCIF2Dict(str(path))

    # OpenFold normalises singletons to lists before reading the header; the
    # header helpers index [0], so the same normalisation is required here.
    for key, value in list(parsed_info.items()):
        if not isinstance(value, list):
            parsed_info[key] = [value]

    return _get_header(parsed_info)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-chains", required=True, type=Path)
    parser.add_argument("--mmcif-dir", required=True, type=Path)
    parser.add_argument("--shard-id", required=True, type=int)
    parser.add_argument("--shard-count", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)

    args = parser.parse_args()

    columns = ["openfold_chain_id", "pdb_id", "retained_sequence"]

    by_entry: dict[str, list[tuple[str, str]]] = defaultdict(list)

    parquet = pq.ParquetFile(args.training_chains)

    for batch in parquet.iter_batches(batch_size=20000, columns=columns):
        data = batch.to_pydict()

        for i, chain in enumerate(data["openfold_chain_id"]):
            by_entry[data["pdb_id"][i]].append(
                (chain, data["retained_sequence"][i])
            )

    entries = sorted(by_entry)
    mine = [
        e for i, e in enumerate(entries)
        if i % args.shard_count == args.shard_id
    ]

    cache: dict[str, dict] = {}
    failures: list[str] = []

    for pdb_id in mine:
        path = args.mmcif_dir / f"{pdb_id}.cif"

        try:
            header = header_for(path)
        except Exception as error:  # noqa: BLE001 - recorded, not swallowed
            failures.append(f"{pdb_id}: {type(error).__name__}: {error}")
            continue

        release_date = header.get("release_date")
        resolution = header.get("resolution")

        for chain, sequence in by_entry[pdb_id]:
            entry = {
                "release_date": release_date,
                "seq": sequence,
                "resolution": resolution,
            }

            # A None resolution makes resolution_filter reject the chain and
            # it would vanish from training without comment. Surface it.
            if resolution is None:
                failures.append(f"{chain}: no resolution in header")

            cache[chain] = entry

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w") as handle:
        json.dump(cache, handle, separators=(",", ":"), sort_keys=True)

    print(f"shard {args.shard_id}/{args.shard_count}")
    print(f"  entries   : {len(mine):,}")
    print(f"  chains    : {len(cache):,}")
    print(f"  failures  : {len(failures):,}")

    for failure in failures[:10]:
        print(f"    {failure}", file=sys.stderr)

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
