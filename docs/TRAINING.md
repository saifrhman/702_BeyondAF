# OpenFold retraining on the PDBClean population

This document records what is being trained, on what, how, and what state the
run is in. It is the companion to [`docs/sequence_redundancy.md`](sequence_redundancy.md),
which explains how the training population was derived.

---

## 1. What is being trained

A full OpenFold model from **random initialisation**. No pretrained AlphaFold
or OpenFold weights are loaded at any point. `--resume_from_ckpt` is used only
to continue *this* run across Slurm wall-time boundaries.

This matters for interpreting the result: the experiment asks whether a model
trained on a geometrically and sequence-deduplicated population behaves
differently from one trained on the raw population. Initialising from released
weights would answer a different question — how quickly a converged model
adapts — and would make the comparison meaningless.

The consequence is cost. A randomly initialised model needs the full training
budget before its predictions are worth comparing, and the numbers below
should be read with that in mind.

---

## 2. The population

| | |
|---|---|
| Release | `PDBClean-…-dedup-v1-seqid100-cov100-exact-v1` |
| Derivation | geometric de-duplication, then sequence-redundancy resolution at identity 1.0 / coverage 1.0 with the exact post-pass |
| Chains | **142,056** — exactly one per distinct sequence |
| Training split | **138,056** |
| Validation split | **4,000** held out |
| Validation actually used | **3,969** (see §5) |

The split is not random with respect to redundancy. It was constructed so that
**no validation chain shares a geometric near-duplicate edge or an identical
sequence with any training chain** — 0 crossing edges, 0 shared sequences.
Length distributions were matched (medians 156 / 157). A split that ignored
this would leak: a validation chain whose near-duplicate sits in training is
not held out in any useful sense.

---

## 3. Optimisation

| Setting | Value |
|---|---|
| Config preset | `initial_training` |
| Precision | bf16 |
| GPUs | 1 × H100 80GB |
| Batch | 1, with `accumulate_grad_batches=8` → global batch 8 |
| Epoch length | 1,024 samples → **128 optimiser steps/epoch** |
| Budget | **400 epochs = 51,200 steps** |
| Seed | 20260101 |
| LR | max 1e-3, 1,000-step warmup, ×0.95 every 1,200 steps after step 25,600 |
| Validation | every 5 epochs, on `val/lddt_ca` |
| Crop | 256 residues (training only) |

The LR schedule is **not** upstream's. Upstream decays every 50,000 steps
after step 50,000, which inside a 51,200-step budget means the rate never
decays at all. The schedule above is pinned in `lr_schedule.json` in the run
root, and the training log prints which file it came from.

**Measured pace: ~56 min/epoch**, so the full 400 epochs is roughly **15.6
days** of H100 time, spread across 3-day Slurm allocations with requeue.

---

## 4. Checkpointing

Three checkpoint callbacks, each answering a different question:

| Callback | Retention | Purpose |
|---|---|---|
| rolling | latest 1 | resume a killed job |
| milestone | every 50 epochs, kept forever | training history; something to fall back on |
| best | top 3 by `val/lddt_ca` | "which model do I actually use" |

Each checkpoint is ~1.4 GB.

> **A trap worth knowing about.** PyTorch Lightning's
> `ModelCheckpoint._should_save_on_train_epoch_end()` returns `False` whenever
> `check_val_every_n_epoch != 1`, which silently moves *every* save to
> validation-end — including a checkpoint configured with `every_n_epochs=1`.
> Setting validation to every 5 epochs therefore turned the rolling resume
> point into a once-every-5-epochs one. The run completed four epochs, died
> inside its first validation, and had written nothing at all. The rolling and
> milestone callbacks now pass `save_on_train_epoch_end=True` explicitly; only
> the best-checkpoint callback is left to follow validation, which is where it
> belongs. See `scripts/openfold_training/openfold_src_checkpoint_timing.patch`.

---

## 5. The validation length cap

Training crops every chain to 256 residues. **OpenFold's validation path does
not crop** — it runs the full chain, and template triangular attention is
roughly cubic in length. One 4,174-residue chain (`8j07_g6`) therefore asked
for **98.31 GiB** on an 80 GB H100 and took the run down at its first
validation pass, six hours in.

The split's median length is 156, but the tail is long: p99 = 739, max = 4,174.

`scripts/openfold_training/cap_validation_chain_length.py` excludes the tail
rather than cropping it, because cropping would change what the reported lDDT
means. At a **768-residue cap**, 3,969 of 4,000 chains remain and the longest
kept is 761. The 31 excluded chains are listed in
`validation_length_cap.json`, and the uncapped cache is kept beside it so the
cap can be raised without rebuilding the split.

Only `val_mmcif_data_cache.json` is pruned: OpenFold enumerates validation
chains from the alignment directory and then drops any chain absent from the
cache, so every symlink stays in place.

> **Pruning the cache is necessary but was not sufficient.** The filter in
> `OpenFoldSingleDataset` runs only `if self.chain_data_cache is not None`, and
> the **monomer** `OpenFoldDataModule` never stored the validation cache at
> all: `train_openfold.py` passes `--val_mmcif_data_cache_path` inside
> `**vars(args)`, the constructor had no matching parameter, and `**kwargs`
> swallowed it silently. The eval dataset was therefore built with
> `chain_data_cache=None`, no filter ran, and all 4,000 chains loaded —
> including the 4,174-residue one. Three runs died at exactly the same
> `98.31 GiB` allocation before this was found. The multimer module already
> threaded the argument through correctly; only the monomer path was broken.
> Fixed in `scripts/openfold_training/openfold_src_val_cache_path.patch`.
>
> The lesson generalises: a filter that silently no-ops when its input is
> `None` will not tell you it did nothing. Verify the *effect* — chain count
> and longest chain actually loaded — not the presence of the filtering code.

---

## 6. Running it

```bash
# one-time: check the environment resolves
bash scripts/pdbclean_doctor.sh

# submit (queues to every GPU partition at or above the run's tier floor)
bash scripts/openfold_training/submit_train_full.sh \
     <run_name> <total_epochs> <epoch_len> <gpus> <accum>
```

A pinned `.resubmit_cmd` is written into the run root and replayed verbatim by
the watchdog, so a rescue always reproduces the original submission rather
than whatever the current shell happens to export.

**One Slurm detail that costs hours if missed.** The `#SBATCH` header names
`gpu-l40s` (`DefMemPerNode=UNLIMITED`) and the submitter overrides the
partition to `gpu-h100` (`DefMemPerCPU=21150`). Both reach the job
environment, and `srun` refuses to start when two of
`SLURM_MEM_PER_{CPU,GPU,NODE}` are set:

```
srun: fatal: SLURM_MEM_PER_CPU, SLURM_MEM_PER_GPU, and SLURM_MEM_PER_NODE
      are mutually exclusive.
```

The job dies in about a second — long enough to look like a scheduling
hiccup, short enough to leave no useful log. `gpu_train_full.sbatch` clears
all three before `srun`.

### The watchdog

`train_watchdog.sbatch` re-submits a stalled run and re-arms itself. Two
things it gets right that are easy to get wrong:

* Its stall threshold is derived from the epoch length
  (`STEPS_PER_EPOCH × 60 × 3`), not a fixed timeout. A fixed 1,800 s threshold
  against a ~56-minute epoch fires mid-epoch, every time, and queues a pile of
  duplicate rescues.
* It gives up after 10 consecutive rescues that make no progress, rather than
  resubmitting a broken job forever. Clearing `.watchdog_stalls` re-arms it.

---

## 7. Current status

As of **2026-09-21**, job `10614448` on gpu31:

* Restarted from random initialisation at 09:58 — the four epochs completed on
  18 September were lost, because no checkpoint had ever been written.
* Epoch 0 reproduced the earlier attempt exactly (loss `85.69731140136719`,
  lDDT-Cα `0.04852814972400665`, identical to 16 significant figures), which
  confirms seed 20260101 is deterministic. The lost epochs cost wall clock,
  not science.

| Epoch | Step | train/loss | train/lDDT-Cα |
|---|---|---|---|
| 0 | 128 | 85.697 | 0.0485 |
| 1 | 256 | 53.158 | 0.0692 |

The checkpoint committed under `checkpoints/` is **epoch 1 of 400**. It is an
early snapshot for provenance and recovery, **not a trained model** — an
lDDT-Cα of 0.069 is barely above random. Its manifest records the exact
population, optimisation settings and metrics it corresponds to.

Meaningful comparisons should use the milestone checkpoints (every 50 epochs)
and the `val/lddt_ca`-selected best checkpoints, once the run has covered
enough of its budget for those numbers to mean anything.
