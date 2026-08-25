#!/usr/bin/env python3
"""Pull real training examples through OpenFold's dataset, on CPU.

The GPU smoke test is the expensive way to discover that the dataloader is
broken.  ``OpenFoldSingleDataset.__getitem__`` does not need a GPU and does not
need the chain-data cache -- that is only read by the sampler -- so the whole
input path can be exercised first: mmCIF parse, retained-chain projection, the
digest self-check, MSA lookup through the alignment index, template stubbing,
cropping and tensor construction.

What this asserts is the invariant the training run depends on: the residue
dimension of the structure features and of the MSA features are the same
number, on real chains, after the feature pipeline has cropped them.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import traceback
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mmcif-dir", required=True, type=Path)
    parser.add_argument("--msa-store", required=True, type=Path)
    parser.add_argument("--alignment-index", required=True, type=Path)
    parser.add_argument("--projection-index", required=True, type=Path)
    parser.add_argument("--train-filter", required=True, type=Path)
    parser.add_argument("--samples", default=12, type=int)
    parser.add_argument("--seed", default=20260101, type=int)
    parser.add_argument("--report", default=None, type=Path)

    args = parser.parse_args()

    from openfold.config import model_config
    from openfold.data.data_modules import OpenFoldSingleDataset

    config = model_config("initial_training", train=True)

    print("crop size        :", config.data.train.crop_size)
    print("max msa clusters :", config.data.train.max_msa_clusters)

    dataset = OpenFoldSingleDataset(
        data_dir=str(args.mmcif_dir),
        alignment_dir=str(args.msa_store),
        template_mmcif_dir=str(args.mmcif_dir),
        max_template_date="2026-01-01",
        config=config.data,
        chain_data_cache_path=None,
        filter_path=str(args.train_filter),
        mode="train",
        alignment_index=json.loads(args.alignment_index.read_text()),
        projection_index=json.loads(args.projection_index.read_text()),
        treat_pdb_as_distillation=False,
    )

    print(f"dataset size     : {len(dataset):,}")

    rng = random.Random(args.seed)
    picks = rng.sample(range(len(dataset)), min(args.samples, len(dataset)))

    failures: list[dict] = []
    ok = 0

    for idx in picks:
        name = dataset.idx_to_chain_id(idx)

        try:
            item = dataset[idx]
        except Exception as error:  # noqa: BLE001 - reported, not swallowed
            failures.append(
                {
                    "chain": name,
                    "error": f"{type(error).__name__}: {error}",
                    "traceback": traceback.format_exc()[-1200:],
                }
            )
            print(f"  {name}: FAILED {type(error).__name__}: {error}")
            continue

        aatype = tuple(item["aatype"].shape)
        msa = tuple(item["msa_feat"].shape)
        positions = tuple(item["all_atom_positions"].shape)
        mask = tuple(item["seq_mask"].shape)

        # aatype is [N_res, N_recycle]; msa_feat is [N_clust, N_res, C, N_recycle]
        n_res = aatype[0]
        msa_res = msa[1]

        agree = n_res == msa_res == positions[0] == mask[0]

        print(
            f"  {name}: aatype={aatype} msa_feat={msa} "
            f"positions={positions} agree={agree}"
        )

        if not agree:
            failures.append(
                {
                    "chain": name,
                    "error": f"residue dims disagree: aatype={n_res} "
                             f"msa={msa_res} positions={positions[0]} "
                             f"mask={mask[0]}",
                }
            )
            continue

        ok += 1

    report = {
        "stage": "openfold_dataloader_preflight",
        "dataset_size": len(dataset),
        "samples_attempted": len(picks),
        "samples_ok": ok,
        "failure_count": len(failures),
        "failures": failures[:10],
        "crop_size": int(config.data.train.crop_size),
        "status": "FAIL" if failures else "PASS",
    }

    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print()
    print(f"samples ok       : {ok}/{len(picks)}")

    if failures:
        print("STATUS: FAIL")

        for failure in failures[:5]:
            print(f"  {failure['chain']}: {failure['error']}", file=sys.stderr)

            if "traceback" in failure:
                print(failure["traceback"], file=sys.stderr)

        return 1

    print("STATUS: PASS")
    print("  real chains load, project, find their MSA, and produce")
    print("  structure and MSA features of the same residue length.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
