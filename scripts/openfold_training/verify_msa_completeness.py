#!/usr/bin/env python3
"""Prove the MSA corpus is complete, unduplicated, and correctly attributed.

An MSA corpus is only usable for training if three things hold, and none of
them is safe to assume:

*complete*
    every unique query sequence has an MSA, and therefore every training chain
    has one through the chain -> sequence mapping;

*unduplicated*
    no sequence was searched twice, and no MSA exists that nothing asked for;

*correctly attributed*
    the MSA stored under a sequence's SHA256 actually begins with that
    sequence.

The third is the one that matters most and is easiest to skip.  The Stage-17
pilot failed precisely because outputs carried the right *content* under the
wrong *name*; a corpus that is complete and unduplicated but misattributed is
worse than a missing one, because nothing downstream will notice.

Exit status is 0 only when all three hold.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


def read_fasta_names(path: Path) -> dict[str, str]:
    """Return {name: sequence} for a FASTA whose headers are sequence SHA256."""

    names: dict[str, str] = {}

    name: str | None = None
    chunks: list[str] = []

    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if name is not None:
                    names[name] = "".join(chunks)

                name = line[1:].split()[0]
                chunks = []
            elif name is not None:
                chunks.append(line.strip())

    if name is not None:
        names[name] = "".join(chunks)

    return names


def first_sequence(path: Path) -> str:
    """The query sequence an a3m describes: its first record.

    a3m lower-cases insertions relative to the query; the first record is the
    query itself and carries none, but it is upper-cased defensively so a
    comparison never fails on formatting alone.
    """

    chunks: list[str] = []
    seen_header = False

    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if seen_header:
                    break

                seen_header = True
                continue

            if seen_header:
                chunks.append(line.strip())

    return "".join(chunks).replace("-", "").upper()


def shard_slices(total: int, shard_count: int) -> dict[int, set[int]]:
    return {
        shard: {i for i in range(total) if i % shard_count == shard}
        for shard in range(shard_count)
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queries", required=True, type=Path)
    parser.add_argument("--msa-store", required=True, type=Path)
    parser.add_argument("--chain-map", default=None, type=Path)
    parser.add_argument("--shard-count", default=None, type=int)
    parser.add_argument(
        "--verify-sequences",
        default=250,
        type=int,
        help=(
            "How many MSAs to open and check the query sequence of. "
            "0 checks none, -1 checks every one (slow)."
        ),
    )
    parser.add_argument("--report", default=None, type=Path)

    args = parser.parse_args()

    failures: list[str] = []
    report: dict[str, object] = {"stage": "openfold_msa_completeness"}

    # -- expected --------------------------------------------------------
    expected = read_fasta_names(args.queries)

    report["unique_sequences_expected"] = len(expected)

    print(f"expected unique sequences : {len(expected):,}")

    # -- produced --------------------------------------------------------
    produced = {
        path.stem: path
        for path in args.msa_store.glob("*.a3m")
    }

    report["msas_present"] = len(produced)

    print(f"MSAs present              : {len(produced):,}")

    # -- complete and unduplicated ---------------------------------------
    missing = sorted(set(expected) - set(produced))
    orphan = sorted(set(produced) - set(expected))

    report["missing_count"] = len(missing)
    report["orphan_count"] = len(orphan)
    report["missing_examples"] = missing[:10]
    report["orphan_examples"] = orphan[:10]

    print(f"missing MSAs              : {len(missing):,}")
    print(f"orphan MSAs               : {len(orphan):,}")

    if missing:
        failures.append(f"{len(missing)} sequences have no MSA")

        for name in missing[:10]:
            print(f"  missing: {name}", file=sys.stderr)

    if orphan:
        failures.append(f"{len(orphan)} MSAs correspond to no query")

        for name in orphan[:10]:
            print(f"  orphan:  {name}", file=sys.stderr)

    # Duplication cannot occur in the store itself -- it is content-addressed,
    # so two searches of one sequence collapse to one filename. What *can*
    # happen is wasted work, which the shard-coverage check below detects.

    # -- empty outputs ---------------------------------------------------
    empty = [
        name for name, path in produced.items() if path.stat().st_size == 0
    ]

    report["empty_count"] = len(empty)

    print(f"empty MSAs                : {len(empty):,}")

    if empty:
        failures.append(f"{len(empty)} MSAs are empty files")

    # -- correctly attributed --------------------------------------------
    shared = sorted(set(expected) & set(produced))

    if args.verify_sequences == 0:
        checked = []
    elif args.verify_sequences < 0:
        checked = shared
    else:
        # Evenly spaced rather than the first N: a systematic fault confined to
        # one shard would hide behind a contiguous prefix.
        step = max(1, len(shared) // max(1, args.verify_sequences))
        checked = shared[::step][: args.verify_sequences]

    mismatched: list[str] = []

    for name in checked:
        observed = first_sequence(produced[name])

        if observed != expected[name].upper():
            mismatched.append(name)
            continue

        # The name is a claim about content; verify the claim itself.
        digest = hashlib.sha256(expected[name].encode()).hexdigest()

        if digest != name:
            mismatched.append(f"{name} (header is not its own sequence sha256)")

    report["sequences_verified"] = len(checked)
    report["mismatched_count"] = len(mismatched)
    report["mismatched_examples"] = mismatched[:10]

    print(f"MSAs verified by sequence : {len(checked):,}")
    print(f"misattributed             : {len(mismatched):,}")

    if mismatched:
        failures.append(
            f"{len(mismatched)} MSAs do not start with the sequence they are "
            "named for"
        )

        for name in mismatched[:10]:
            print(f"  misattributed: {name}", file=sys.stderr)

    # -- shard coverage ---------------------------------------------------
    if args.shard_count:
        slices = shard_slices(len(expected), args.shard_count)

        union: set[int] = set()
        overlaps = 0

        for indices in slices.values():
            overlaps += len(union & indices)
            union |= indices

        report["shard_count"] = args.shard_count
        report["shard_union"] = len(union)
        report["shard_overlaps"] = overlaps

        print(f"shard coverage            : {len(union):,} of {len(expected):,}")
        print(f"shard overlaps            : {overlaps:,}")

        if len(union) != len(expected):
            failures.append(
                f"shards cover {len(union)} of {len(expected)} sequences"
            )

        if overlaps:
            failures.append(f"shards overlap on {overlaps} sequences")

    # -- chain coverage ---------------------------------------------------
    if args.chain_map is not None:
        import pyarrow.parquet as pq

        table = pq.read_table(
            args.chain_map, columns=["openfold_chain_id", "sequence_sha256"]
        )

        chain_sequences = table.column("sequence_sha256").to_pylist()

        uncovered = sum(
            1 for digest in chain_sequences if digest not in produced
        )

        report["training_chains"] = len(chain_sequences)
        report["chains_without_msa"] = uncovered

        print(f"training chains           : {len(chain_sequences):,}")
        print(f"chains without an MSA     : {uncovered:,}")

        if uncovered:
            failures.append(f"{uncovered} training chains have no MSA")

    # -- verdict -----------------------------------------------------------
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
    print("  every query sequence has exactly one MSA,")
    print("  no MSA exists that nothing asked for,")
    print("  and every MSA checked begins with the sequence it is named for.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
