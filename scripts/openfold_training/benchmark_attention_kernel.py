#!/usr/bin/env python3
"""Is DeepSpeed's fused Evoformer attention faster here, and does it agree?

This OpenFold version enables ``use_deepspeed_evo_attention`` only under
``long_sequence_inference``, which asserts ``not train`` -- so training always
runs the unfused attention path, even though the kernel supports training and
newer OpenFold releases default to it.

Switching kernels is not a scientific change: it is the same attention
arithmetic with a different reduction order, like changing a matmul backend.
But "not a scientific change" is a claim, so it is measured rather than
asserted: the same chains, in the same order, are run through both paths and
the losses compared. A kernel that is fast and wrong is worthless.

Reports the speedup and the loss agreement so the decision rests on numbers.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import torch


def build(config, checkpoint):
    from openfold.model.model import AlphaFold

    model = AlphaFold(config)

    if checkpoint is not None:
        state = torch.load(checkpoint, map_location="cpu", weights_only=False)["state_dict"]
        model.load_state_dict(
            {k[len("model."):]: v for k, v in state.items() if k.startswith("model.")},
            strict=False,
        )

    return model.cuda()


def run(model, loss_fn, batches, warmup=1):
    """Return (losses, seconds_per_step) over identical prepared batches."""

    losses = []
    times = []

    for i, batch in enumerate(batches):
        batch = {k: v.cuda() for k, v in batch.items()}

        torch.cuda.synchronize()
        start = time.time()

        out = model(batch)
        final = {k: v[..., -1] for k, v in batch.items()}
        loss, _ = loss_fn(out, final, _return_breakdown=True)
        loss.backward()

        torch.cuda.synchronize()
        elapsed = time.time() - start

        model.zero_grad(set_to_none=True)

        losses.append(loss.item())

        if i >= warmup:            # first step pays JIT and allocator costs
            times.append(elapsed)

    return losses, (sum(times) / len(times) if times else float("nan"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mmcif-dir", required=True, type=Path)
    parser.add_argument("--msa-store", required=True, type=Path)
    parser.add_argument("--alignment-index", required=True, type=Path)
    parser.add_argument("--projection-index", required=True, type=Path)
    parser.add_argument("--train-filter", required=True, type=Path)
    parser.add_argument("--checkpoint", default=None, type=Path)
    parser.add_argument("--steps", default=4, type=int)
    parser.add_argument("--seed", default=20260101, type=int)
    parser.add_argument("--report", default=None, type=Path)

    args = parser.parse_args()

    from openfold.config import model_config
    from openfold.data.data_modules import OpenFoldBatchCollator, OpenFoldSingleDataset
    from openfold.utils.loss import AlphaFoldLoss

    base = model_config("initial_training", train=True, low_prec=True)

    dataset = OpenFoldSingleDataset(
        data_dir=str(args.mmcif_dir),
        alignment_dir=str(args.msa_store),
        template_mmcif_dir=str(args.mmcif_dir),
        max_template_date="2026-01-01",
        config=base.data,
        filter_path=str(args.train_filter),
        mode="train",
        alignment_index=json.loads(args.alignment_index.read_text()),
        projection_index=json.loads(args.projection_index.read_text()),
        treat_pdb_as_distillation=False,
    )

    collator = OpenFoldBatchCollator()

    # Identical batches for both paths -- otherwise the loss comparison is noise.
    import random

    rng = random.Random(args.seed)
    picks = [rng.randrange(len(dataset)) for _ in range(args.steps)]
    names = [dataset.idx_to_chain_id(i) for i in picks]

    print("chains:", ", ".join(names))

    torch.manual_seed(args.seed)
    batches = [collator([dataset[i]]) for i in picks]

    results = {}

    for label, use_kernel in (("unfused (current)", False), ("deepspeed_evo", True)):
        config = model_config("initial_training", train=True, low_prec=True)
        config.globals.use_deepspeed_evo_attention = use_kernel

        print(f"\n=== {label}: use_deepspeed_evo_attention={use_kernel} ===")

        try:
            torch.manual_seed(args.seed)
            model = build(config, args.checkpoint)
            loss_fn = AlphaFoldLoss(config.loss)

            torch.cuda.reset_peak_memory_stats()
            losses, per_step = run(model, loss_fn, batches)
            peak = torch.cuda.max_memory_allocated() / (1024 ** 3)

            print(f"  losses     : {[round(l, 4) for l in losses]}")
            print(f"  s / step   : {per_step:.2f}")
            print(f"  peak GiB   : {peak:.2f}")

            results[label] = {
                "losses": losses,
                "seconds_per_step": per_step,
                "peak_gib": peak,
                "ok": all(torch.isfinite(torch.tensor(l)) for l in losses),
            }
        except Exception as error:  # noqa: BLE001
            print(f"  FAILED: {type(error).__name__}: {error}")
            results[label] = {"error": f"{type(error).__name__}: {error}"}

        del model
        torch.cuda.empty_cache()

    print("\n=== verdict ===")

    a = results.get("unfused (current)", {})
    b = results.get("deepspeed_evo", {})

    if "error" in b:
        print("  fused kernel unavailable here; keep the current path")
        verdict = "kernel_unavailable"
    else:
        speedup = a["seconds_per_step"] / b["seconds_per_step"]
        drift = max(abs(x - y) for x, y in zip(a["losses"], b["losses"]))
        rel = drift / max(1e-9, max(abs(x) for x in a["losses"]))

        print(f"  speedup        : {speedup:.2f}x")
        print(f"  peak memory    : {a['peak_gib']:.2f} -> {b['peak_gib']:.2f} GiB")
        print(f"  max loss drift : {drift:.4f}  ({rel * 100:.3f}% of loss)")

        # Same maths, different reduction order: agreement should be tight.
        verdict = "adopt" if (speedup > 1.15 and rel < 0.02) else "keep_current"
        print(f"  verdict        : {verdict}")

        results["speedup"] = speedup
        results["max_loss_drift"] = drift
        results["relative_drift"] = rel

    results["verdict"] = verdict

    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
