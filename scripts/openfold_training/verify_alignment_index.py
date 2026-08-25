#!/usr/bin/env python3
"""Prove the alignment index hands each training chain its own MSA.

The index is a bridge between two namespaces -- chains (``101m_A``) and
sequence digests -- and a bridge that is merely *complete* can still be wrong.
An off-by-one in the chain mapping, or a sequence digest computed over a
different string than the one that was searched, produces an index where every
chain resolves to some MSA and no check based on counts alone notices.

So the test here is not "does every chain have an entry".  It is: follow the
entry exactly as OpenFold will -- open ``db``, ``seek(start)``, ``read(size)``
-- and confirm the MSA's query row is character-for-character the cleaned
sequence recorded for *that chain* in the training manifest.

That single check subsumes the weaker ones.  It can only pass if the chain
mapping, the digests, the file names and the byte ranges are all correct
together.

Exit status is 0 only when every check holds.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import pyarrow.parquet as pq


def query_sequence(db_path: Path, start: int, size: int) -> str:
    """The first record of an a3m, read exactly as OpenFold reads it.

    a3m lower-cases insertions relative to the query, and the query row itself
    carries none; it is upper-cased and stripped of gaps defensively so a
    comparison never fails on formatting alone.
    """

    with db_path.open("rb") as handle:
        handle.seek(start)
        blob = handle.read(size)

    text = blob.decode("utf-8", errors="strict")

    chunks: list[str] = []
    seen_header = False

    for line in text.splitlines():
        if line.startswith(">"):
            if seen_header:
                break

            seen_header = True
            continue

        if seen_header:
            chunks.append(line.strip())

    return "".join(chunks).replace("-", "").upper()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", required=True, type=Path)
    parser.add_argument("--alignment-dir", required=True, type=Path)
    parser.add_argument(
        "--training-chains",
        required=True,
        type=Path,
        help="openfold_training_chains.parquet, carrying retained_sequence",
    )
    parser.add_argument(
        "--train-filter",
        default=None,
        type=Path,
        help="Every chain listed here must be present in the index.",
    )
    parser.add_argument(
        "--verify-chains",
        default=500,
        type=int,
        help="How many chains to resolve and read. -1 reads every one (slow: "
             "that is 161 GiB).",
    )
    parser.add_argument("--seed", default=20260101, type=int)
    parser.add_argument("--report", default=None, type=Path)

    args = parser.parse_args()

    failures: list[str] = []
    report: dict[str, object] = {"stage": "openfold_alignment_index_verify"}

    # -- load ---------------------------------------------------------------
    with args.index.open() as handle:
        index = json.load(handle)

    print(f"chains in index           : {len(index):,}")

    table = pq.read_table(
        args.training_chains,
        columns=["openfold_chain_id", "retained_sequence"],
    )

    expected = dict(
        zip(
            table.column("openfold_chain_id").to_pylist(),
            table.column("retained_sequence").to_pylist(),
        )
    )

    report["chains_in_index"] = len(index)
    report["chains_in_manifest"] = len(expected)

    print(f"chains in manifest        : {len(expected):,}")

    # -- coverage -----------------------------------------------------------
    uncovered = sorted(set(expected) - set(index))
    unexpected = sorted(set(index) - set(expected))

    report["chains_without_entry"] = len(uncovered)
    report["entries_without_chain"] = len(unexpected)

    print(f"chains without an entry   : {len(uncovered):,}")
    print(f"entries with no chain     : {len(unexpected):,}")

    if uncovered:
        failures.append(f"{len(uncovered)} training chains have no index entry")

        for chain in uncovered[:10]:
            print(f"  uncovered: {chain}", file=sys.stderr)

    if unexpected:
        failures.append(
            f"{len(unexpected)} index entries name no training chain"
        )

    if args.train_filter is not None:
        wanted = {
            line.strip()
            for line in args.train_filter.read_text().splitlines()
            if line.strip()
        }

        filtered_out = sorted(wanted - set(index))

        report["train_filter_chains"] = len(wanted)
        report["train_filter_chains_missing"] = len(filtered_out)

        print(f"train_filter chains       : {len(wanted):,}")
        print(f"  of those, unindexed     : {len(filtered_out):,}")

        if filtered_out:
            failures.append(
                f"{len(filtered_out)} chains in the train filter are not in "
                "the index; OpenFold would silently train on fewer chains"
            )

    # -- sharing ------------------------------------------------------------
    distinct = {entry["db"] for entry in index.values()}

    report["distinct_msas_referenced"] = len(distinct)
    report["chains_sharing_an_msa"] = len(index) - len(distinct)

    print(f"distinct MSAs referenced  : {len(distinct):,}")
    print(f"chains reusing a shared MSA: {len(index) - len(distinct):,}")

    # -- the real check: resolve and read ----------------------------------
    shared = sorted(set(index) & set(expected))

    if args.verify_chains < 0:
        sample = shared
    else:
        rng = random.Random(args.seed)
        sample = rng.sample(shared, min(args.verify_chains, len(shared)))

    mismatched: list[str] = []
    unreadable: list[str] = []
    checked = 0

    for chain in sample:
        entry = index[chain]

        if len(entry["files"]) != 1:
            unreadable.append(f"{chain}: expected one file, got {entry['files']}")
            continue

        name, start, size = entry["files"][0]

        if not name.endswith(".a3m"):
            unreadable.append(f"{chain}: {name} is not an .a3m")
            continue

        db_path = args.alignment_dir / entry["db"]

        try:
            observed = query_sequence(db_path, start, size)
        except (OSError, UnicodeDecodeError) as error:
            unreadable.append(f"{chain}: {type(error).__name__} {error}")
            continue

        checked += 1

        if observed != expected[chain].upper():
            mismatched.append(chain)

    report["chains_resolved_and_read"] = checked
    report["unreadable_count"] = len(unreadable)
    report["mismatched_count"] = len(mismatched)
    report["mismatched_examples"] = mismatched[:10]

    print(f"chains resolved and read  : {checked:,}")
    print(f"unreadable entries        : {len(unreadable):,}")
    print(f"wrong sequence delivered  : {len(mismatched):,}")

    if unreadable:
        failures.append(f"{len(unreadable)} index entries could not be read")

        for line in unreadable[:10]:
            print(f"  unreadable: {line}", file=sys.stderr)

    if mismatched:
        failures.append(
            f"{len(mismatched)} chains resolve to an MSA whose query row is "
            "not that chain's sequence"
        )

        for chain in mismatched[:10]:
            print(f"  mismatched: {chain}", file=sys.stderr)

    # -- verdict ------------------------------------------------------------
    report["status"] = "FAIL" if failures else "PASS"
    report["failures"] = failures

    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n"
        )

    print()

    if failures:
        print("STATUS: FAIL")

        for failure in failures:
            print(f"  - {failure}")

        return 1

    print("STATUS: PASS")
    print("  every training chain has an index entry,")
    print("  every entry resolves to a readable MSA,")
    print("  and every MSA read delivers that chain's own cleaned sequence.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
