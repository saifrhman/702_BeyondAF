"""Walk every validation chain through the real eval pipeline and report all failures.

The validation path has produced five distinct failures in this project, each
discovered by losing two hours of H100 time to it: Lightning tying checkpoints
to the validation schedule, uncropped full-length inference exhausting HBM, the
monomer data module silently dropping the chain-data cache, the missing
projection index featuring the deposited polymer instead of the retained one,
and a single-residue chain collapsing seq_length to a 0-dim tensor.

Every one of those was a data- or wiring-level fault that a dataloader pass
would have caught in minutes. This script is that pass. It builds the eval
dataset exactly as OpenFoldDataModule.setup() builds it -- same class, same
arguments -- and pulls every chain through featurisation, collecting failures
rather than stopping at the first.

It touches no GPU and trains nothing; it needs a CUDA node only because
openfold.config imports DeepSpeed, which shells out to nvcc on import.

Read-only with respect to the release, the mmCIF tree and the MSA store.
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
import traceback
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--val-data-dir", required=True, type=Path)
    ap.add_argument("--val-alignment-dir", required=True, type=Path)
    ap.add_argument("--val-cache", required=True, type=Path)
    ap.add_argument("--val-projection-index", required=True, type=Path)
    ap.add_argument("--template-mmcif-dir", required=True, type=Path)
    ap.add_argument("--max-template-date", default="2026-01-01")
    ap.add_argument("--config-preset", default="initial_training")
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after N chains (0 = every chain)")
    ap.add_argument("--report", type=Path, default=None)
    args = ap.parse_args()

    from openfold.config import model_config
    from openfold.data.data_modules import OpenFoldSingleDataset

    config = model_config(args.config_preset, train=True)

    with open(args.val_projection_index) as handle:
        projection_index = json.load(handle)

    # Exactly the construction in OpenFoldDataModule.setup(); if this drifts,
    # the preflight stops testing what actually runs.
    dataset = OpenFoldSingleDataset(
        data_dir=str(args.val_data_dir),
        alignment_dir=str(args.val_alignment_dir),
        template_mmcif_dir=str(args.template_mmcif_dir),
        max_template_date=args.max_template_date,
        config=config.data,
        chain_data_cache_path=str(args.val_cache),
        projection_index=projection_index,
        filter_path=None,
        max_template_hits=config.data.eval.max_template_hits,
        mode="eval",
    )

    total = len(dataset)
    print(f"eval dataset chains : {total:,}", flush=True)

    n = total if args.limit in (0, None) else min(args.limit, total)
    failures: list[tuple[str, str, str]] = []
    longest = (0, None)
    kinds: collections.Counter = collections.Counter()

    for i in range(n):
        name = dataset.idx_to_chain_id(i) if hasattr(dataset, "idx_to_chain_id") \
               else dataset._chain_ids[i]
        try:
            feats = dataset[i]
            L = int(feats["aatype"].shape[0])
            if L > longest[0]:
                longest = (L, name)
        except Exception as exc:                      # noqa: BLE001 - the point
            kind = type(exc).__name__
            kinds[kind] += 1
            failures.append((name, kind, str(exc)[:200]))

        if (i + 1) % 250 == 0:
            print(f"  {i+1:>6,} / {n:,}   failures so far: {len(failures)}",
                  flush=True)

    print(f"\nchains walked       : {n:,}")
    print(f"failures            : {len(failures):,}")
    print(f"longest featurised  : {longest[0]:,} residues ({longest[1]})")

    if kinds:
        print("\nfailures by kind")
        for kind, count in kinds.most_common():
            print(f"  {kind:<28} {count:>6,}")
        print("\nfirst 25 failing chains")
        for name, kind, msg in failures[:25]:
            print(f"  {name:<12} {kind:<24} {msg}")

    if args.report:
        args.report.write_text(json.dumps({
            "chains_walked": n,
            "failure_count": len(failures),
            "longest_featurised_residues": longest[0],
            "longest_featurised_chain": longest[1],
            "failures_by_kind": dict(kinds),
            "failures": [{"chain": c, "kind": k, "message": m}
                         for c, k, m in failures],
        }, indent=2, sort_keys=True))
        print(f"\nreport written      : {args.report}")

    if failures:
        print("\nVALIDATION SET IS NOT CLEAN -- fix these before spending GPU time.")
        return 1

    print("\nVALIDATION SET IS CLEAN -- every chain featurises.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception:
        traceback.print_exc()
        sys.exit(2)
