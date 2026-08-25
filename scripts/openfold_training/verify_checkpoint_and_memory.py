#!/usr/bin/env python3
"""Reload a smoke checkpoint and measure what one real training step costs.

Two things the smoke run itself does not prove.  That its checkpoint can be
loaded back -- a checkpoint that cannot be resumed is not a checkpoint, and a
multi-day run split across Slurm wall-time boundaries depends entirely on it.
And how much GPU memory a step actually needs, which decides the batch and
worker settings for the full run rather than being guessed from the model size.

The step here is real: a chain drawn from the training dataset, through the
retained-chain projection and its own MSA, forward, loss, backward.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import torch


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--mmcif-dir", required=True, type=Path)
    parser.add_argument("--msa-store", required=True, type=Path)
    parser.add_argument("--alignment-index", required=True, type=Path)
    parser.add_argument("--projection-index", required=True, type=Path)
    parser.add_argument("--train-filter", required=True, type=Path)
    parser.add_argument("--steps", default=2, type=int)
    parser.add_argument("--seed", default=20260101, type=int)
    parser.add_argument("--report", default=None, type=Path)

    args = parser.parse_args()

    from openfold.config import model_config
    from openfold.data.data_modules import OpenFoldBatchCollator, OpenFoldSingleDataset
    from openfold.model.model import AlphaFold
    from openfold.utils.loss import AlphaFoldLoss

    config = model_config("initial_training", train=True, low_prec=True)

    # -- reload -------------------------------------------------------------
    print(f"loading checkpoint: {args.checkpoint}")

    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)

    print(f"  checkpoint keys : {sorted(checkpoint)[:8]}")
    print(f"  global_step     : {checkpoint.get('global_step')}")
    print(f"  epoch           : {checkpoint.get('epoch')}")

    state = checkpoint["state_dict"]
    model_state = {
        k[len("model."):]: v for k, v in state.items() if k.startswith("model.")
    }

    model = AlphaFold(config)
    missing, unexpected = model.load_state_dict(model_state, strict=False)

    print(f"  missing keys    : {len(missing)}")
    print(f"  unexpected keys : {len(unexpected)}")

    if missing or unexpected:
        print("  FAILED: state dict does not match the model")
        return 1

    print("  checkpoint reloaded into AlphaFold cleanly")

    model = model.cuda()
    loss_fn = AlphaFoldLoss(config.loss)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, eps=1e-5)

    # -- one real step ------------------------------------------------------
    dataset = OpenFoldSingleDataset(
        data_dir=str(args.mmcif_dir),
        alignment_dir=str(args.msa_store),
        template_mmcif_dir=str(args.mmcif_dir),
        max_template_date="2026-01-01",
        config=config.data,
        filter_path=str(args.train_filter),
        mode="train",
        alignment_index=json.loads(args.alignment_index.read_text()),
        projection_index=json.loads(args.projection_index.read_text()),
        treat_pdb_as_distillation=False,
    )

    collator = OpenFoldBatchCollator()
    rng = random.Random(args.seed)

    torch.cuda.reset_peak_memory_stats()

    losses = []

    for step in range(args.steps):
        idx = rng.randrange(len(dataset))
        name = dataset.idx_to_chain_id(idx)

        batch = collator([dataset[idx]])
        batch = {k: v.cuda() for k, v in batch.items()}

        output = model(batch)

        # The loss reads the final recycling iteration of each feature.
        batch = {k: v[..., -1] for k, v in batch.items()}

        loss, _ = loss_fn(output, batch, _return_breakdown=True)

        optimizer.zero_grad()
        loss.backward()

        grads_finite = all(
            torch.isfinite(p.grad).all().item()
            for p in model.parameters()
            if p.grad is not None
        )

        optimizer.step()

        value = loss.item()
        losses.append(value)

        print(
            f"  step {step}: chain={name} loss={value:.4f} "
            f"finite={torch.isfinite(loss).item()} grads_finite={grads_finite}"
        )

        if not torch.isfinite(loss).item() or not grads_finite:
            print("  FAILED: non-finite loss or gradients")
            return 1

    peak_gib = torch.cuda.max_memory_allocated() / (1024 ** 3)
    reserved_gib = torch.cuda.max_memory_reserved() / (1024 ** 3)

    print()
    print(f"peak GPU allocated : {peak_gib:.2f} GiB")
    print(f"peak GPU reserved  : {reserved_gib:.2f} GiB")

    report = {
        "stage": "openfold_checkpoint_reload_and_memory",
        "checkpoint": str(args.checkpoint),
        "global_step": checkpoint.get("global_step"),
        "epoch": checkpoint.get("epoch"),
        "reload_ok": True,
        "steps_run": args.steps,
        "losses": losses,
        "all_losses_finite": True,
        "peak_gpu_allocated_gib": round(peak_gib, 3),
        "peak_gpu_reserved_gib": round(reserved_gib, 3),
        "device": torch.cuda.get_device_name(0),
        "status": "PASS",
    }

    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print()
    print("STATUS: PASS")
    print("  the smoke checkpoint reloads, and real forward/backward/optimizer")
    print("  steps run on it with finite loss and finite gradients.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
