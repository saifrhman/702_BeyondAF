#!/usr/bin/env python3
"""Combine per-shard training-view reports into one verdict.

A sharded validation is only evidence if every shard ran.  A missing shard and
a passing shard look identical in a directory listing, so the expected shard
count is required rather than inferred.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports", required=True, type=Path)
    parser.add_argument("--shard-count", required=True, type=int)
    parser.add_argument("--expected-chains", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--exclusions",
        default=None,
        type=Path,
        help=(
            "Exclusion record. Failures are tolerated only when every failing "
            "chain is named there, so the verdict covers the population that "
            "will actually be trained and nothing is quietly forgiven."
        ),
    )

    args = parser.parse_args()

    failures: list[str] = []

    excluded: dict[str, str] = {}

    if args.exclusions is not None:
        with args.exclusions.open() as handle:
            excluded = json.load(handle).get("excluded", {})

    present = sorted(args.reports.glob("shard_*.json"))
    seen: dict[int, dict] = {}

    for path in present:
        with path.open() as handle:
            report = json.load(handle)

        seen[report["shard_id"]] = report

    missing = [i for i in range(args.shard_count) if i not in seen]

    if missing:
        failures.append(f"{len(missing)} shards produced no report: {missing[:10]}")

    totals = {
        "chains_checked": 0,
        "sequence_matches": 0,
        "terminal_trimmed_checked": 0,
        "auth_differs_from_label_checked": 0,
        "feature_checks": 0,
        "feature_matches": 0,
        "failure_count": 0,
    }

    examples: list[dict] = []

    for report in seen.values():
        for key in totals:
            totals[key] += report.get(key, 0)

        examples.extend(report.get("failures", []))

    if totals["chains_checked"] != args.expected_chains:
        failures.append(
            f"checked {totals['chains_checked']} chains, expected "
            f"{args.expected_chains}"
        )

    unaccounted = [
        failure["chain"]
        for failure in examples
        if failure["chain"] not in excluded
    ]

    mismatched = totals["chains_checked"] - totals["sequence_matches"]

    if unaccounted:
        failures.append(
            f"{len(unaccounted)} chains failed and are not recorded as "
            f"exclusions: {unaccounted[:10]}"
        )

    if totals["feature_matches"] != totals["feature_checks"]:
        failures.append(
            f"{totals['feature_checks'] - totals['feature_matches']} feature "
            "builds disagreed in width"
        )

    trainable = totals["chains_checked"] - mismatched

    if excluded and mismatched != len(excluded):
        failures.append(
            f"{mismatched} chains mismatched but {len(excluded)} are recorded "
            "as excluded; the two must agree"
        )

    summary = {
        "stage": "openfold_training_view_validation_summary",
        "shard_count": args.shard_count,
        "shards_reported": len(seen),
        "expected_chains": args.expected_chains,
        **totals,
        "mismatched_chains": mismatched,
        "excluded_chains": len(excluded),
        "trainable_population": trainable,
        "failure_examples": examples[:25],
        "failures": failures,
        "status": "FAIL" if failures else "PASS",
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    print(f"shards reported           : {len(seen)}/{args.shard_count}")
    print(f"chains checked            : {totals['chains_checked']:,}")
    print(f"projected == retained     : {totals['sequence_matches']:,}")
    print(f"terminal-trimmed covered  : {totals['terminal_trimmed_checked']:,}")
    print(f"auth != label covered     : {totals['auth_differs_from_label_checked']:,}")
    print(f"feature-width checks      : {totals['feature_matches']:,}/"
          f"{totals['feature_checks']:,}")
    print(f"mismatched chains         : {mismatched:,}")
    print(f"recorded exclusions       : {len(excluded):,}")
    print(f"trainable population      : {trainable:,}")
    print()

    if failures:
        print("STATUS: FAIL")

        for failure in failures:
            print(f"  - {failure}")

        return 1

    print("STATUS: PASS")
    print(f"  {trainable:,} chains project to exactly their Gold sequence,")
    print("  every feature build agrees in width with its MSA,")

    if excluded:
        print(f"  and the {len(excluded)} chains that do not are recorded as")
        print("  explicit exclusions rather than trained on.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
