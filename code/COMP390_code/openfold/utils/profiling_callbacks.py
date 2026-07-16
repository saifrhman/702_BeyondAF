# openfold/utils/profiling_callbacks.py
from __future__ import annotations

import json
import os
import time
from typing import Any, Dict, Optional

import pytorch_lightning as pl

import csv, os


class EpochSummaryCSVLogger(pl.Callback):
    """
    Write ONE line per epoch with fixed columns (human-readable).
    """
    def __init__(self, out_dir: str, filename: str = "epoch_summary.csv"):
        super().__init__()
        self.path = os.path.join(out_dir, filename)
        os.makedirs(out_dir, exist_ok=True)
        self._header_written = False

        # 你想看的字段（按你要的顺序排）
        self.columns = [
            "epoch", "global_step",
            "train/loss_epoch", "train/lddt_ca", "train/drmsd_ca",
            "perf/data_time_s", "perf/compute_time_s",
        ]

    def _get_scalar(self, trainer: pl.Trainer, key: str):
        v = trainer.callback_metrics.get(key, None)
        if v is None:
            return ""
        try:
            return float(v.item()) if hasattr(v, "item") else float(v)
        except Exception:
            return ""

    def on_train_epoch_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule) -> None:
        # 只在 rank0 写，避免 DDP 多进程重复写
        if trainer.is_global_zero:
            write_header = (not os.path.exists(self.path)) or (os.path.getsize(self.path) == 0)
            with open(self.path, "a", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=self.columns)
                if write_header:
                    writer.writeheader()

                row = {
                    "epoch": int(trainer.current_epoch),
                    "global_step": int(trainer.global_step),
                }
                for k in self.columns:
                    if k in ("epoch", "global_step"):
                        continue
                    row[k] = self._get_scalar(trainer, k)

                writer.writerow(row)


class DataLoadComputeTimer(pl.Callback):
    """
    Roughly separates:
      - data_time: time between previous batch end and current batch start
      - compute_time: time spent inside the batch (train_step + backward/optim step dominated)

    Notes:
      - In Lightning, the next batch is typically prefetched before on_train_batch_start.
        So "data_time" here approximates dataloader latency / pipeline stalls between steps.
      - This is good enough to detect 'GPU waiting for data' (data_time comparable to compute_time).
    """
    def __init__(self, log_every_n_steps: int = 50):
        super().__init__()
        self.log_every_n_steps = log_every_n_steps
        self._prev_batch_end_t: Optional[float] = None
        self._batch_start_t: Optional[float] = None

    def on_train_batch_start(self, trainer: pl.Trainer, pl_module: pl.LightningModule, batch: Any, batch_idx: int) -> None:
        now = time.perf_counter()
        if self._prev_batch_end_t is None:
            data_time = 0.0
        else:
            data_time = now - self._prev_batch_end_t

        self._batch_start_t = now

        # log step-level
        if trainer.global_step % self.log_every_n_steps == 0:
            # logger=True only works if you have a logger; otherwise it is a no-op.
            pl_module.log("perf/data_time_s", data_time, on_step=True, on_epoch=False, prog_bar=False, logger=True)

    def on_train_batch_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule, outputs: Any, batch: Any, batch_idx: int) -> None:
        now = time.perf_counter()
        if self._batch_start_t is None:
            compute_time = 0.0
        else:
            compute_time = now - self._batch_start_t

        self._prev_batch_end_t = now

        if trainer.global_step % self.log_every_n_steps == 0:
            pl_module.log("perf/compute_time_s", compute_time, on_step=True, on_epoch=False, prog_bar=False, logger=True)
            # ratio helps spot stalls quickly
            ratio = (compute_time / (compute_time + 1e-12)) if (compute_time + 0.0) == 0 else compute_time / (compute_time + 1e-12)
            # Above ratio is not super meaningful; better log data/compute ratio:
            # keep it simple:
            # pl_module.log("perf/data_over_compute", data_time / (compute_time + 1e-12), ...)

    def on_train_epoch_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule) -> None:
        # If you want epoch-level averages, easiest is to rely on logger aggregations,
        # or extend this callback to keep running sums.
        pass


class EpochJSONLLogger(pl.Callback):
    """
    Writes one JSON line per epoch (append-only) so nothing gets overwritten by tqdm.
    """
    def __init__(self, out_dir: str, filename: str = "epoch_metrics.jsonl"):
        super().__init__()
        self.out_dir = out_dir
        self.filename = filename
        self.path = os.path.join(out_dir, filename)
        os.makedirs(out_dir, exist_ok=True)

    def _dump(self, trainer: pl.Trainer, stage: str) -> None:
        # trainer.callback_metrics contains epoch-level metrics (tensors + floats)
        row: Dict[str, Any] = {"stage": stage, "epoch": int(trainer.current_epoch), "global_step": int(trainer.global_step)}
        for k, v in trainer.callback_metrics.items():
            try:
                if hasattr(v, "item"):
                    row[k] = float(v.item())
                else:
                    row[k] = float(v)
            except Exception:
                # non-scalar / non-numeric -> skip
                continue

        with open(self.path, "a") as f:
            f.write(json.dumps(row) + "\n")

    def on_train_epoch_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule) -> None:
        self._dump(trainer, stage="train_epoch_end")

    def on_validation_epoch_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule) -> None:
        self._dump(trainer, stage="val_epoch_end")