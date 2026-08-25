#!/usr/bin/env python3
"""Derive the trainable filter, recording every chain held back and why.

Validation found a small number of chains whose projected sequence differs from
the Gold retained sequence at identical length.  The cause is upstream: entries
with point microheterogeneity carry two ``_entity_poly_seq`` rows for one
residue number (MET/MSE, THR/AEI, PYR/SER), OpenFold's ``_get_protein_chains``
appends both, and every seqres index past that point shifts by one.  OpenFold
therefore cannot build a correct feature set for those chains at all -- with or
without this pipeline's projection.

They are excluded rather than repaired.  Repairing would mean changing
OpenFold's core polymer parsing for every structure to fix a handful, and the
alternative to exclusion is training a chain against an MSA built for a
different sequence.

The exclusion is written down, with the chain list, so the trained population is
auditable and never silently smaller than the Gold population.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-filter", required=True, type=Path)
    parser.add_argument("--validation-reports", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--exclusions", required=True, type=Path)

    args = parser.parse_args()

    chains = [
        line.strip()
        for line in args.train_filter.read_text().splitlines()
        if line.strip()
    ]

    excluded: dict[str, str] = {}

    for path in sorted(args.validation_reports.glob("shard_*.json")):
        with path.open() as handle:
            report = json.load(handle)

        for failure in report.get("failures", []):
            excluded[failure["chain"]] = failure["kind"]

    kept = [c for c in chains if c not in excluded]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(kept) + "\n")

    record = {
        "stage": "openfold_train_filter_view",
        "stage_version": "1.0",
        "gold_population": len(chains),
        "trainable_population": len(kept),
        "excluded_count": len(excluded),
        "excluded": dict(sorted(excluded.items())),
        "exclusion_reason": (
            "OpenFold's _get_protein_chains appends every _entity_poly_seq row, "
            "so entries with point microheterogeneity (two mon_id values at one "
            "residue number) yield a seqres one residue too long and every "
            "index past that point shifts. The projected sequence then has the "
            "correct length and the wrong residues, which no shape check "
            "detects. These chains are held back rather than trained against "
            "an MSA built for a different sequence."
        ),
        "train_filter_view": str(args.output.resolve()),
        "train_filter_view_sha256": hashlib.sha256(
            args.output.read_bytes()
        ).hexdigest(),
    }

    args.exclusions.parent.mkdir(parents=True, exist_ok=True)
    args.exclusions.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")

    print(f"gold population       : {len(chains):,}")
    print(f"excluded              : {len(excluded):,}")
    print(f"trainable population  : {len(kept):,}")
    print(f"filter written        : {args.output}")
    print(f"exclusions recorded   : {args.exclusions}")

    for chain, kind in sorted(excluded.items()):
        print(f"  excluded {chain} ({kind})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
