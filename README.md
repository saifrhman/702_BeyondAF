# Beyond AlphaFold: Filling the Blind Spots in Protein Structure Prediction

COMP702 dissertation project — Saif Ur Rehman, University of Liverpool.

This repository implements **PDBClean**: a reproducible pipeline that takes an
immutable PDB snapshot, cleans it under an explicit quality protocol, computes
an exact geometric representation of every protein backbone, detects
geometrically duplicate chains, and publishes a deduplicated chain dataset with
complete provenance.

It also contains the downstream OpenFold work that consumes that dataset.

**Status labels used throughout this document**

| Label | Meaning |
|-------|---------|
| **[IMPLEMENTED]** | In the current code and exercised by the test suite. |
| **[FROZEN]** | A published COMP702 result. Immutable. |
| **[IN PROGRESS]** | Being worked on now. Not finished. |
| **[FUTURE]** | Planned. **Not implemented. No results exist.** |

---

## Quickstart — running this on your own Barkla account

Nothing in this repository names an individual account. Every path is derived
at run time from `$USER` and from wherever you cloned it, so the setup is a
clone, an environment, and a preflight check.

**1. Clone it somewhere with space.** The repository is ~6 GB and the working
data is much larger, so put it on fastscratch rather than in `$HOME`:

```bash
cd /mnt/fastscratch/users/$USER          # or $HOME/fastscratch
git clone https://github.com/saifrhman/702_BeyondAF.git
cd 702_BeyondAF
```

Large binaries (checkpoints, reference structures) are stored with Git LFS.
If `git lfs` is not installed, the clone still works but those files arrive as
small pointer text files:

```bash
git lfs install && git lfs pull         # only if you need the binaries
```

**2. Build the environment** from the pinned spec:

```bash
module load apps/anaconda3               # or however conda is provided
conda env create -p "$HOME/fastscratch/envs/bri_env_1.2.2" \
                 -f reproducibility/bri_environment.yml
```

**3. Check it resolves.** This is the step that saves you time — it reports
every path the pipeline will use and prints the command that fixes anything
missing. It writes nothing:

```bash
bash scripts/pdbclean_doctor.sh
```

```
PDBClean preflight
  user         abcd1234
  repository   /mnt/fastscratch/users/abcd1234/702_BeyondAF
Filesystem roots (derived from $USER; override any in the environment)
  OK    PDBCLEAN_SCRATCH_ROOT          /mnt/scratch/users/abcd1234
  OK    PDBCLEAN_FASTSCRATCH_ROOT      /mnt/fastscratch/users/abcd1234
PDBClean pipeline
  OK    python                         .../bri_env_1.2.2/bin/python (3.10.20)
  OK    pdbclean package               importable from $PDBCLEAN_REPO_ROOT/src
  OK    slurm                          /usr/bin/sbatch  partition nodes
  13 ok, 0 note, 0 missing
  Ready.
```

Fix every `MISS` before going further. A `NOTE` is only a warning — the
OpenFold rows are all optional, and the PDBClean pipeline runs without them.

**4. Look before you run.** Both commands are safe on a login node and change
nothing:

```bash
source config/pdbclean/pipeline_env.sh
pdbclean stages                          # the canonical pipeline, in order
pdbclean plan --config config/pdbclean/profiles/comp702_frozen_20260101.yaml
```

`plan` prints what each stage *would* do, which inputs it needs and what it
would publish. Read it before submitting anything.

**5. Run a stage.** Locally for the cheap stages, Slurm for the rest:

```bash
# local, single machine
pdbclean run --config <profile> --stage <name> --executor local

# Slurm array (partition and array width come from the environment)
pdbclean run --config <profile> --stage <name> --executor slurm
```

Two profiles ship with the repository:

| Profile | Publishes |
|---|---|
| `comp702_frozen_20260101.yaml` | geometry only — the frozen 499,770-chain population |
| `comp702_seqclust_20260101.yaml` | geometry, then sequence-redundancy resolution |

**6. Verify your install** (optional but quick):

```bash
pytest tests -q          # 1,234 tests, ~7 minutes
```

### If your account is laid out differently

Every variable takes the form `${VAR:-default}`, so exporting one beforehand
always wins. The two that matter most:

```bash
export PDBCLEAN_SCRATCH_ROOT=/some/big/disk/$USER
export PDBCLEAN_FASTSCRATCH_ROOT=/some/fast/disk/$USER
```

The full table is in [§5](#running-this-on-barkla-as-a-different-user). To
supply a completely separate environment file, export `PDBCLEAN_ENV_FILE` and
every wrapper will source that instead.

### Retraining OpenFold

That layer is optional and much heavier — a GPU allocation, an MSA store and
the OpenFold source. It has its own document:
[`docs/TRAINING.md`](docs/TRAINING.md).

---

## Contents

- [Quickstart — running this on your own Barkla account](#quickstart--running-this-on-your-own-barkla-account) — clone, environment, preflight
1. [Motivation](#1-motivation)
2. [The scientific pipeline](#2-the-scientific-pipeline)
3. [Stage-to-code map](#3-stage-to-code-map)
4. [Repository architecture](#4-repository-architecture)
5. [Installation](#5-installation)
    - [Running this on Barkla as a different user](#running-this-on-barkla-as-a-different-user) — every path is derived from `$USER`
6. [Running the pipeline](#6-running-the-pipeline)
7. [Configuration architecture](#7-configuration-architecture)
8. [Snapshot selection and preservation](#8-snapshot-selection-and-preservation)
9. [Bronze / Silver / Gold](#9-bronze--silver--gold)
10. [Provenance and reproducibility](#10-provenance-and-reproducibility)
11. [The web UI](#11-the-web-ui)
12. [Historical run workflow](#12-historical-run-workflow)
13. [Frozen COMP702 result](#13-frozen-comp702-result-frozen)
    - 13.1 [Threshold sensitivity](#131-threshold-sensitivity-implemented) — τ and Brain-threshold results
14. [Testing](#14-testing)
15. [OpenFold training view and retraining](#15-openfold-training-view-and-retraining)
16. [Current status](#16-current-status)
17. [Future work](#17-future-work)
18. [Authority order](#18-authority-order)

Deeper detail lives in [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md),
[`docs/CONFIGURATION.md`](docs/CONFIGURATION.md),
[`docs/PROVENANCE.md`](docs/PROVENANCE.md),
[`docs/REPOSITORY_MAP.md`](docs/REPOSITORY_MAP.md),
[`docs/sequence_redundancy.md`](docs/sequence_redundancy.md) and
[`docs/TRAINING.md`](docs/TRAINING.md). This README is
self-sufficient for understanding and operating the project; those documents
expand on it rather than repeating it.

---

## 1. Motivation

AlphaFold and its open reimplementations are trained and evaluated on the
Protein Data Bank. The PDB is highly redundant: the same backbone conformation
appears many times over. Redundancy inflates apparent accuracy, biases
training, and hides the cases where prediction actually fails — the blind
spots.

This project attacks that from the data side. Rather than measuring redundancy
by **sequence identity**, it measures it **geometrically and exactly**, using
the Backbone Rigid Invariant (BRI), and produces a deduplicated dataset in
which every removal decision is individually justified and auditable.

The deduplicated dataset is then intended as the training population for
OpenFold retraining — see [§15](#15-openfold-training-view-and-retraining),
[§16](#16-current-status) and [§17](#17-future-work)
for exactly how far that has and has not progressed.

---

## 2. The scientific pipeline

### 2.1 End-to-end flow [IMPLEMENTED]

```
snapshot resolution                 pin to a concrete immutable snapshot
        ↓
Bronze source inventory             PDB ID, S3 key, size, ETag
        ↓
deterministic parsing (Silver)      reconstructed on demand, not stored
        ↓
model selection                     model 1
        ↓
structural cleaning                 Protocol 3.2 rules Q001–Q006
        ↓
geometric validation                quarantine degenerate backbones
        ↓
complete BRI                        the final geometric representation
        ↓
precision-grid representation       exact integer representation units
        ↓
Brain                               9-D average BRI — filtering layer only
        ↓
exact chain-length buckets          different lengths are incomparable
        ↓
Brain candidate filtering           lossless prefilter (cKDTree)
        ↓
complete-BRI NN search              exact L∞ radius search (cover tree)
        ↓
complete-BRI L∞ classification      exact / near duplicate
        ↓
investigation and validation        Stages 11–13, review only
        ↓
direct-edge-safe representative     Stage 14
   selection
        ↓
retained Gold dataset               the published release
        ↓
sequence-redundancy resolution      Stage 14d, OPTIONAL -- off unless asked
        ↓
sequence-reduced release            a separate release identifier
```

The last two steps are optional and off by default. A run that does not enable
them publishes the geometry-only population; enabling
`sequence_clustering.enabled` publishes both, the second derived from the
first and under its own identifier. See
[`docs/sequence_redundancy.md`](docs/sequence_redundancy.md).

### 2.2 Complete BRI

The **complete Backbone Rigid Invariant** (BRI v1.2.2) is the final geometric
representation of a chain: an *m*×9 matrix of rigid-motion-invariant backbone
coordinates. Two chains are only ever compared through their complete BRI.

This is the project's own computational representation. The BRI *definition*
comes from the MATCH work; the exact integer arithmetic, the pipeline around
it, and every threshold choice are COMP702 engineering decisions and are not
prescribed by any paper.

### 2.3 Representation precision *p*

Complete BRI is represented on a configurable **precision grid**:

```
BRI_units = round(BRI / p)
```

The validated COMP702 default is **p = 0.001 Å**, at which one representation
unit is one milliångström and the rule is exactly `BRI_mA = round(1000 × BRI)`.

`p` is the precision at which geometry is *recorded*. It is **not** a duplicate
threshold. See [§7.3](#73-precision-p-versus-threshold-τ).

In the executable code, p = 0.001 Å is currently **structural**: `compute_bri`
ends with `numpy.around(..., 3)` to reproduce the pinned BRI v1.2.2
canonicalisation, and every downstream representation assumes exact integer
milliångströms. p is exposed as configuration so that a precision study is an
explicitly identified scientific configuration — but a run configured for a
different grid is refused by the production stages rather than silently
producing v1.2.2 output under a different label.

### 2.4 Brain

**Brain** is the 9-dimensional average BRI vector (MATCH Definition 5.1): the
column means of the BRI matrix excluding its first row. It is defined only for
chains with *m* ≥ 2; chains with *m* = 1 have no Brain vector.

Brain is the **filtering and indexing layer only**. It never classifies
duplicates.

### 2.5 Duplicate detection

1. **Exact chain-length grouping.** BRI matrices of different lengths are not
   comparable, so comparison only ever happens inside one exact-*m* bucket.
   This is a correctness requirement, not an optimisation.

2. **Brain filtering.** Within a bucket, a **lossless** prefilter — SciPy
   `cKDTree` with `p=inf, eps=0`, plus an exact integer post-filter — at
   L∞ ≤ 0.010 Å. Because the bucket shares a common denominator, the threshold
   is exactly `tau_units × (m − 1)` in integer sum units. `cKDTree` is used
   **only** here, never for final classification.

3. **Complete-BRI nearest-neighbour search.** The Elkin–Kurlin **compressed
   cover tree** performs the exact complete-BRI L∞ radius search over the Brain
   candidates. This is the production search engine.

4. **Classification.** Final classification is always complete-BRI L∞:

   | Class | Criterion |
   |-------|-----------|
   | exact duplicate | `d == 0` |
   | near duplicate | `d ≤ τ`, where τ = 0.010 Å, **inclusive** |

   The comparison is `≤`, not `<`. A pair at exactly 10 units *is* a near
   duplicate.

### 2.6 Redundancy resolution (Stage 14)

A near-duplicate **graph** is built over the detected pairs, but:

* a connected component is **not** a duplicate equivalence class — two chains
  in one component may be far apart geometrically;
* there is **no** transitive removal;
* **every removed chain must have its own direct `d ≤ τ` edge to the chain it
  was assigned to.**

Representative selection walks each component with a deterministic
quality-ordered greedy **direct-edge cover**. Chains with *m* = 1 are all
retained.

This removes redundancy **of shape only**. Two chains with an identical
sequence but different coordinates are not near-duplicates under this
criterion and both survive, so the retained population still carries sequence
redundancy by design — 499,770 retained chains hold 142,056 distinct
sequences. Removing that is the separate, optional job of
[§2.6.1](#261-sequence-redundancy-resolution-stage-14d-optional).

### 2.6.1 Sequence-redundancy resolution (Stage 14d, optional)

Clusters the Stage-14c survivors on their **retained (post-trimming)**
sequences with MMseqs2 and keeps one chain per cluster.

MMseqs2 decides cluster *membership only*. The survivor is chosen by the same
deterministic ranking Stage 14b applies to geometric components, with the
helpers imported from that entry point rather than copied — the clusterer's own
greedy set cover is order-dependent and was observed placing byte-identical
sequences under different representatives.

Two parameters carry the scientific weight:

* `min_seq_id` — 1.0 by default. Removes only chains indistinguishable from
  their representative; anything looser discards homologues a structure model
  has reason to see.
* `coverage` / `cov_mode` — bidirectional, so a short chain cannot be absorbed
  into a long one on a shared domain. **At 0.8 the operation is not exact**: it
  merges length variants, and in the executed 0.8 run it merged `2OLO` into
  `2OLN` and deleted one of the two observed conformations of a fold switcher.
  At 1.0, with the exact post-pass enabled, the result is exactly one chain per
  distinct sequence.

`docs/sequence_redundancy.md` records the threshold sweep, the validation
gates, the determinism check and the population counts.

### 2.7 What the result is, and is not

The output is **a geometrically deduplicated PDB chain dataset under COMP702
representative policy v1**.

It is *not* a claim that any removed chain corresponds to an incorrect or
invalid experiment.

---

## 3. Stage-to-code map

The canonical scientific vocabulary is **Prerequisites A–C, then Stage 1
through Stage 14**, with Stage 14 realised by its subdivisions 14a, 14b, 14c
and the optional 14d. Prerequisites are lettered so they can never be mistaken
for scientific stages. This is the single canonical table; other sections
reference it rather than repeating it.

| Canonical stage | Purpose | Primary module / producer | Slurm wrapper | Output | Layer |
|---|---|---|---|---|---|
| **Prerequisite A** — Snapshot resolution | Pin to a concrete immutable snapshot | `pdbclean.snapshot_selection` | — | pinned identity in provenance | snapshot |
| **Prerequisite B** — Bronze source manifest | Immutable source inventory | `scripts/pdbclean/create_manifest.py` | — | `source_manifest.parquet` | bronze |
| **Prerequisite C** — Silver parsed representation | Deterministic parsing | `pdbclean.mmcif_parser` | — | *not persisted by design* | silver |
| **Stage 1** — Structural cleaning | Protocol 3.2 rules Q001–Q006 | `pdbclean.cleaning` via `scripts/pdbclean/run_quality_task.py` | `scripts/pdbclean/run_quality_array.sbatch`, `submit_quality_pipeline.sh` | `quality/merged/accepted.parquet` | gold |
| **Stage 2** — Geometric validation | Quarantine degenerate backbones | `pdbclean.geometric_validation`, `scripts/pdbclean/finalize_geometric_validation.py` | `run_geometric_validation_array.sbatch` | `finalized/eligible.parquet` | gold |
| **Stage 3** — Complete BRI | Compute complete BRI | `pdbclean.bri`, `scripts/pdbclean/finalize_bri.py` | `run_bri_array.sbatch` | `bri/finalized/bri.parquet` | gold |
| **Stage 4** — BRI numerical representation | Represent BRI on the precision grid | *same producer as Stage 3* (`numpy.around(..., 3)` inside `compute_bri`); integer conversion in `pdbclean.full_bri_compare` | *same* | *same artefact* | gold |
| **Stage 5** — Brain | 9-D average BRI | `pdbclean.brain`, `pdbclean.brain_finalize_cli` | — | `brain/finalized/brain.parquet` | gold |
| **Stage 6** — Exact chain-length grouping | Partition into exact-*m* buckets | `pdbclean.length_buckets_cli` | — | `finalized/bucket_index.parquet` | gold |
| **Stage 7** — Brain candidate filtering | Lossless same-length prefilter | `pdbclean.brain_prefilter`, `pdbclean.brain_prefilter_production` | — | `finalized/candidates.parquet` | gold |
| **Stage 8** — Complete-BRI NN search | Exact L∞ radius search | `pdbclean.compressed_cover_tree`, `pdbclean.full_bri_nn_production` | — | `finalized/candidate_near_duplicates.parquet` | gold |
| **Stage 9** — Complete-BRI distance representation | Authoritative distance representation | *same producer as Stage 8* — the search emits represented distances | *same* | *same artefact* | gold |
| **Stage 10** — Duplicate classification | Exact / near classification | `pdbclean.duplicate_classification`, `..._production` | — | `finalized/candidate_classifications.parquet` | gold |
| **Stage 11** — Acta-style downstream investigation | Review pass | *not orchestrated* | — | `acta_downstream_investigation_v2/` | gold |
| **Stage 12** — Scientific validation gates | Validation evidence | *not orchestrated* | — | `acta_manual_review_manifest_v2/` | gold |
| **Stage 13** — Detailed investigation / review | Manual review subset | *not orchestrated* | — | `acta_detailed_review_v2/` | gold |
| **Stage 14a** — Geometric redundancy graph | Build the near-duplicate graph | `scripts/build_stage14_geometric_graph.py` | `task_scripts/run_stage14_geometric_graph.sbatch` | `stage14_geometric_graph/` | gold |
| **Stage 14b** — Representative selection | Direct-edge cover | `scripts/select_stage14_representatives.py` | `task_scripts/run_stage14_representatives.sbatch` | `representative_mapping.parquet` | gold |
| **Stage 14c** — Final Gold release | Publish the retained dataset | `scripts/build_stage14_final_release.py` | `task_scripts/run_stage14_final_release.sbatch` | `data/retained_chains.parquet` | gold |
| **Stage 14d** — Sequence-redundancy resolution *(optional)* | Cluster the survivors by sequence, keep one per cluster | `pdbclean.sequence_clustering_production` | *single job* | `sequence_clustering/` + its own release | gold |

**Shared producers.** Stages 3 and 4 share one implementation, because
`compute_bri` applies the precision grid at the point of computation. Stages 8
and 9 share one implementation, because the search emits its distances already
represented. Both canonical identities are preserved everywhere — in the
registry, the UI and provenance — because they are distinct scientific
concepts.

**Stage 14a/b/c** are *engineering subdivisions* of scientific Stage 14, not
separate scientific stages. **Stage 14d** is a fourth, and differs from the
others in being **optional**: it is disabled by default, and a run that leaves
it off publishes the geometry-only population and nothing else. When enabled it
consumes the Stage-14c release read-only and publishes a second release under
its own identifier, so the two populations coexist rather than one replacing
the other. See [`docs/sequence_redundancy.md`](docs/sequence_redundancy.md).

**Stages 11–13** are investigation and validation passes. They are **not** on
the release path and are **never** a deletion relation. `Stage 13` in
particular is a manual review subset and must never be used as the global
Stage-14 deletion set.

`pdbclean stages` prints this table from the registry that the planner and UI
both use.

---

## 4. Repository architecture

| Path | Status | Role |
|------|--------|------|
| `src/pdbclean/` | ACTIVE | The pipeline package: the scientific implementation plus the configuration, orchestration, provenance, inspection and UI layers. |
| `src/pdbclean/ui/` | ACTIVE | Stdlib HTTP server and static assets for the web UI. No framework dependency. |
| `config/pdbclean/` | mixed | `protocol_3_2_comp702_v1.yaml` and `stage14_representative_policy_v1.yaml` are **FROZEN and byte-immutable** (their SHA256s are embedded in frozen provenance). `profiles/` and `pipeline_env.sh` are ACTIVE. |
| `scripts/pdbclean/` | ACTIVE | Per-stage entry points and Slurm array wrappers for Stages 1–2 and the finalisers. |
| `scripts/` (top level) | mixed | `build_stage14_*.py` and `select_stage14_representatives.py` are ACTIVE (Stages 14a–c). `01_prepare_pdb707k.py`, `extract_pdb707k_bri_vectors.py`, `analyze_2olo_*.py` are LEGACY PDB707K/COMP390 analysis. |
| `scripts/openfold_training/` | IN PROGRESS | OpenFold relaxation and BRI LAI work. Not part of PDBClean. |
| `task_scripts/` | ACTIVE | Slurm wrappers for Stages 14a–c and the scientific regression harness. Every argument is derived from the resolved run configuration. |
| `tests/pdbclean/` | TESTING | The full suite, including the scientific regression layer. |
| `docs/` | DOCUMENTATION | Architecture, configuration, provenance, repository map, the pipeline specification, the development status log and [`sequence_redundancy.md`](docs/sequence_redundancy.md) (Stage 14d). |
| `docs/provenance/` | FROZEN | Release provenance for the 20260101 publication and the Acta review. |
| `outputs/pdbclean/<snapshot>/<protocol>/` | GENERATED / FROZEN | Stage outputs. The 20260101 tree is frozen. Gitignored (large). |
| `outputs/releases/` | FROZEN | Published Gold releases. Immutable. |
| `outputs/runs/` | GENERATED | Run provenance directories. Append-only. |
| `outputs/snapshot_store/` | GENERATED | Durable snapshot preservation — see [§8.2](#82-durable-preservation-versus-hot-cache-implemented). |
| `outputs/snapshot_cache/` | GENERATED | Disposable working cache. Safe to delete; see [§8.2](#82-durable-preservation-versus-hot-cache-implemented). |
| `reports/` | mixed | Acta review CSVs and evidence (FROZEN); `molstar_exact_duplicate_examples/` is ACTIVE and wired into the UI. |
| `reproducibility/` | FROZEN | Pinned environment exports and `bri_version.txt`. |
| `tools/` | FROZEN | The pinned BRI v1.2.2 reference implementation used by the differential gate. |
| `reference/acta_2025/` | DOCUMENTATION | Wlodawer et al., *Acta Cryst D* 2025 (doi 10.1107/S2059798325001883). |
| `data/` | LEGACY | PDB707K-era inputs. The 126 MB sequence table itself lives in `~/COMP702_BeyondAF/data/`, which is where `scripts/01_prepare_pdb707k.py` reads it from. |
| `code/COMP390_code/` | POINTER | Minhao's COMP390 dissertation work is **not** vendored here. It lives in `~/COMP702_BeyondAF/code/COMP390_code/`, which is what `config/comp702_paths.sh` and the audit scripts resolve to. Only `LEGACY.md` and one sbatch that exists nowhere else are kept. |
| `sbatch/` | LEGACY / HISTORICAL | COMP390-era batch scripts. See its `LEGACY.md`. |
| `logs/` | GENERATED | Slurm output. |

See [`docs/REPOSITORY_MAP.md`](docs/REPOSITORY_MAP.md) for the module-level
breakdown.

---

## 5. Installation

### Environment

The pipeline needs Python ≥ 3.10 with `numpy`, `scipy`, `pyarrow`, `pandas`,
`gemmi` and `PyYAML`. On Barkla the pinned environment is
`$PDBCLEAN_FASTSCRATCH_ROOT/envs/bri_env_1.2.2` — that is, under your own
account — and is exported in `reproducibility/bri_environment.yml`:

```bash
conda env create -p "$PDBCLEAN_FASTSCRATCH_ROOT/envs/bri_env_1.2.2" \
                 -f reproducibility/bri_environment.yml
```

```bash
# editable install, provides the `pdbclean` entry point
pip install -e .

# with test dependencies
pip install -e ".[test]"
```

### Running this on Barkla as a different user

Nothing in the executable surface of this repository names an individual
account. Every path is derived at run time from `$USER` and from where the
clone happens to sit, by a single file:

```
config/pdbclean/pipeline_env.sh
```

Each shell and Slurm wrapper sources it on the line after its `#SBATCH`
block, so cloning the repository and running it is the whole setup:

```bash
git clone <url> && cd COMP702_pdbclean_pipeline
bash scripts/pdbclean_doctor.sh
```

`pdbclean_doctor.sh` resolves every path the pipeline will use, reports
whether it exists, and for anything missing prints the command that fixes it.
It writes nothing, and exits non-zero while the pipeline cannot yet run. Run
it first; it is faster than discovering a missing directory inside a GPU
allocation.

**What is derived, and what you can override.** Every variable below follows
`${VAR:-<default>}`, so exporting it beforehand always wins. Setting none of
them is the expected case on Barkla.

| Variable | Default | What it is |
|---|---|---|
| `PDBCLEAN_REPO_ROOT` | the clone's own location | Derived from the env file's path, so a clone anywhere works |
| `PDBCLEAN_SCRATCH_ROOT` | `$HOME/scratch`, else `/mnt/scratch/users/$USER` | Bulk storage |
| `PDBCLEAN_FASTSCRATCH_ROOT` | `$HOME/fastscratch`, else `/mnt/fastscratch/users/$USER` | Working storage |
| `PDBCLEAN_PYTHON` | the pinned env's interpreter if present, else `python` | A bare `python` on a compute node has none of the dependencies |
| `PDBCLEAN_CONDA_ENV` | `$PDBCLEAN_FASTSCRATCH_ROOT/envs/bri_env_1.2.2` | Rebuild from `reproducibility/bri_environment.yml` |
| `PDBCLEAN_OUTPUT_ROOT` / `_RELEASE_ROOT` / `_RUN_ROOT` | under `$PDBCLEAN_REPO_ROOT/outputs/` | Created on first run |
| `PDBCLEAN_SLURM_PARTITION` | `nodes` | Change for a different cluster |
| `OPENFOLD_SRC`, `OPENFOLD_TRAIN_ENV`, `OPENFOLD_RUNS_ROOT`, `OPENFOLD_MSA_STORE`, `OPENFOLD_MMCIF_DIR`, `MMSEQS_BIN` | under the two scratch roots | Retraining only; the PDBClean pipeline does not need them |

Two notes on the Barkla-specific parts. `$HOME/scratch` and `$HOME/fastscratch`
are the conventional symlinks, but the `/mnt/...` paths are used directly where
they are absent, so both layouts work. And `TMPDIR`, `TRITON_CACHE_DIR` and
`TORCH_EXTENSIONS_DIR` are redirected onto fastscratch, because left on a
compute node's local `/tmp` they fill the disk and the job dies mid-epoch.

To keep the repository on one filesystem and the working data on another, or
to run the whole thing outside Barkla, point the two roots wherever you like:

```bash
export PDBCLEAN_SCRATCH_ROOT=/some/big/disk/$USER
export PDBCLEAN_FASTSCRATCH_ROOT=/some/fast/disk/$USER
bash scripts/pdbclean_doctor.sh
```

An entirely separate environment file can be supplied instead, and the
wrappers will source that rather than the default:

```bash
export PDBCLEAN_ENV_FILE=/path/to/my_pipeline_env.sh
```

One thing is deliberately *not* rewritten: the JSON manifests and logs under
`outputs/` and `logs/` still contain absolute paths under the account that
produced them. Those are provenance records of runs that actually happened —
rewriting them would falsify the record. They are inputs to no code path.

### Without installing

Every command also works straight from the source tree:

```bash
PYTHONPATH=src python -m pdbclean.cli <subcommand>
```

Both forms are used interchangeably below; `pdbclean X` and
`PYTHONPATH=src python -m pdbclean.cli X` are equivalent.

### Verify

```bash
pdbclean stages          # prints the canonical pipeline
pdbclean config          # prints the resolved configuration
pytest tests -q          # runs the full suite
```

### External tools

* **Slurm** (`sbatch`) for HPC execution.
* **MMseqs2** for Stage 18 MSA generation — downstream, see
  [§15](#15-openfold-training-view-and-retraining).
* Mol\* is loaded in the browser from a CDN by the pair viewer; no local
  install is required.

---

## 6. Running the pipeline

### Inspect (safe on a login node)

```bash
pdbclean snapshots                       # list available snapshots, newest first
pdbclean snapshots --limit 10
pdbclean config                          # resolved configuration + per-value sources
pdbclean config --config config/pdbclean/profiles/comp702_frozen_20260101.yaml
pdbclean plan                            # what would run, what would be reused
pdbclean stages                          # the canonical stage chain and its gates
pdbclean stage-command --stage redundancy_graph --shell   # one stage's exact argv
```

### Dry-run behaviour

**`pdbclean run` is dry-run by default.** It resolves configuration, pins the
snapshot, prints the fully resolved scientific configuration and the plan, asks
for confirmation, creates the run directory and writes provenance — and then
prints the command each outstanding stage *would* execute without executing
anything. Nothing runs until you choose an executor.

### Interactive run

```bash
pdbclean run --interactive
```

Presents the snapshot menu (position 1 is "latest complete snapshot
[default]", so pressing Enter takes the default), then the resolved
configuration, then the plan, then a confirmation prompt.

### Non-interactive / config-driven run

```bash
pdbclean run --config config/pdbclean/profiles/comp702_frozen_20260101.yaml --yes
pdbclean run --snapshot 2026-04-15 --yes
pdbclean run --set duplicate_search.near_duplicate_threshold_angstrom=0.005 --yes
pdbclean run --plan-only            # resolve and plan, create nothing
```

### Choosing which population a run publishes

Two profiles ship, and the only scientific difference between them is whether
Stage 14d runs:

```bash
# geometry only -- the frozen 499,770-chain population
pdbclean run --config config/pdbclean/profiles/comp702_frozen_20260101.yaml --yes

# geometry, then sequence-redundancy resolution -- publishes a second release
pdbclean run --config config/pdbclean/profiles/comp702_seqclust_20260101.yaml --yes
```

With the stage disabled the planner reports it as *not applicable* rather than
as work forever outstanding, so a geometry-only run plans clean. The switch is
a single key:

```bash
pdbclean run --set sequence_clustering.enabled=true \
             --set sequence_clustering.coverage=1.0 \
             --set sequence_clustering.exact_post_pass=true --yes
```

### Local execution

```bash
pdbclean run --executor local --yes
```

### Slurm / HPC execution

```bash
pdbclean run --executor slurm --yes        # submits each outstanding stage with sbatch
```

or submit the Stage-14 wrappers directly:

```bash
sbatch task_scripts/run_stage14_geometric_graph.sbatch
sbatch task_scripts/run_stage14_representatives.sbatch
sbatch task_scripts/run_stage14_final_release.sbatch
sbatch task_scripts/run_stage14_regression.sbatch
```

Submit from the repository root, or export `PDBCLEAN_REPO_ROOT` — Slurm copies
batch scripts into its spool directory, so the repository cannot be inferred
from `BASH_SOURCE` inside a running job. Select a different profile with
`PDBCLEAN_PROFILE=<file> sbatch ...`.

**Barkla policy.** Login nodes are for lightweight inspection and orchestration
only. All heavy work runs on compute nodes through `sbatch`; batch scripts are
never executed directly. Arrays use a small physical-worker array that strides
the logical work rather than one enormous array.

### Inspect duplicates from the CLI

```bash
pdbclean duplicates --summary
pdbclean duplicates --exact-only --limit 20
pdbclean duplicates --pdb-id 1a0t --chain A
pdbclean duplicates --relationship removed --limit 50
pdbclean duplicates --min-distance 5 --max-distance 10
pdbclean duplicates --json --limit 5
```

### Inspect historical runs

```bash
pdbclean status                                   # list every recorded run
pdbclean status run-20260821T033704Z-ad0801a7     # one run, with its stage table
pdbclean status <run-id> --json
```

For the full stage-by-stage drill-down, use the UI —
see [§12](#12-historical-run-workflow).

### Reproduce the frozen 2026-01-01 configuration

```bash
pdbclean plan --config config/pdbclean/profiles/comp702_frozen_20260101.yaml
```

This resolves to scientific hash
`25b8e62a87cb90797af41cd4149dfd4280e3a7aed99428e70fa97117c5bababa` and reports
every stage as reusable against the frozen outputs.

### Start the UI

```bash
pdbclean ui                                      # http://127.0.0.1:8765/
PYTHONPATH=src python -m pdbclean.cli ui         # equivalent, no install

pdbclean ui --port 9000
pdbclean ui --no-browser                         # do not open a browser
pdbclean ui --config config/pdbclean/profiles/comp702_frozen_20260101.yaml
```

* Default bind address: **127.0.0.1** (loopback only).
* Default port: **8765**.

**Remote use.** The server binds to loopback, so reach it over an SSH tunnel
from your workstation:

```bash
ssh -N -L 8765:127.0.0.1:8765 USER@BARKLA_LOGIN_HOST
```

then open `http://127.0.0.1:8765/` locally. Replace `USER` and
`BARKLA_LOGIN_HOST` with your own credentials.

**Stop / restart.** Press `Ctrl-C` in the terminal running `pdbclean ui`. To
restart on a different port, stop it and start again with `--port`.

### Open duplicate pairs in Mol\*

In the UI, open **Duplicate Explorer**, filter to the pairs you want, and use
the **View pair** link on any row that has a prepared scene. Prepared scenes
live in `reports/molstar_exact_duplicate_examples/`.

Mol\* is for **human inspection only**. It never determines whether two chains
are duplicates — that comes from the complete-BRI L∞ calculation alone.

---

## 7. Configuration architecture

### 7.1 Precedence

```
validated built-in defaults  →  configuration profile  →  explicit CLI/UI override
```

A later layer replaces an earlier one; mappings merge recursively. The resolver
records which layer supplied **every leaf**, so provenance can answer "where
did this value come from" for any parameter.

The **fully resolved** configuration — not the input file — is what a run
executes, and it is resolved exactly **once**, at run creation, then persisted.
Downstream work loads that document rather than re-resolving it.

### 7.2 Two hashes

| Hash | Covers | Meaning |
|------|--------|---------|
| `resolved_config_sha256` | the whole canonical document | exact run identity |
| `scientific_config_sha256` | the scientific projection only | *scientific* run identity |

Both are **host-independent**. Infrastructure (`storage`, `execution`,
`observability`) and dataset expectation gates are excluded from the scientific
projection, so choosing `--executor slurm` does not look like a different
experiment — while changing a threshold, a rule, the model scope, the snapshot,
the representative policy or the representation precision always does.

**Theme is not configuration.** The UI light/dark preference is stored in the
browser's `localStorage`, never sent to the backend, and cannot affect either
hash.

### 7.3 Precision *p* versus threshold *τ*

These are two **distinct experimental axes** and are never conflated:

| | Question it asks |
|---|---|
| **p** — representation precision | How sensitive is BRI-based redundancy detection to the numerical precision at which backbone geometry is represented? |
| **τ** — near-duplicate threshold | How sensitive are redundancy relationships, and downstream model behaviour, to the geometric near-duplicate threshold? |

Both are configurable, independently, and each produces a distinct scientific
identity.

The τ question has been answered at the frozen snapshot, together with a
measurement of what the Brain prefilter threshold costs — results, tables and
reproduction commands in
[§13.1](#131-threshold-sensitivity-implemented). The *p* question remains open,
because the executable stages implement only p = 0.001 Å; see
[§17](#17-future-work) item D.

**Grid compatibility.** A threshold must be an exact whole number of
representation units: `τ / p` must be an integer. `p = 0.001, τ = 0.010` → 10
units; `p = 0.002, τ = 0.010` → 5 units; `p = 0.003, τ = 0.010` is **rejected**
with a clear error rather than silently rounded.

### 7.4 Fixed by definition

Refused as overrides in both the CLI and the UI: `brain.dimension = 9`
(Definition 5.1), `duplicate_search.metric = L_infinity`,
`final_classification_basis = complete_BRI`, and both operators
`less_than_or_equal`. The Brain prefilter may also never be tighter than the
classifier — that would break its losslessness guarantee.

Full detail, including the superseded historical `geometric_search` block:
[`docs/CONFIGURATION.md`](docs/CONFIGURATION.md).

---

## 8. Snapshot selection and preservation

### 8.1 Selection

The default is the **latest complete snapshot**. Any historical snapshot stays
reproducible by naming it; both date forms are accepted everywhere:

```bash
pdbclean plan                              # latest complete
pdbclean plan --snapshot 2026-04-15
pdbclean plan --snapshot 20260415
```

`latest_complete` is only a **selection mode**. It resolves *once*, at the
start of a run, to a concrete identity:

```
latest_complete  →  2026-04-15
```

From that point the run's identity is the concrete snapshot. Provenance retains
both facts — `selection_mode: latest_complete` and
`resolved_snapshot: 20260415` — and a resumed or reproduced run never
re-resolves "latest" and silently switches snapshots.

### 8.2 Durable preservation versus hot cache [IMPLEMENTED]

Upstream and local caches expire. A completed scientific run must not become
unreproducible because of that, so storage separates two layers:

```
PDB snapshot source (S3)
        ↓  resolve to a concrete snapshot identity
immutable Bronze manifest
        ↓  preserve by content identity  →  DURABLE
outputs/snapshot_store/objects/<content-id>
outputs/snapshot_store/snapshots/<YYYYMMDD>.manifest.json
        ↓  materialise for computation  →  DISPOSABLE
outputs/snapshot_cache/<YYYYMMDD>/
        ↓
PDBClean stages
```

* The **durable layer** is content-addressed. An object unchanged between two
  snapshots is preserved **once** and referenced by both snapshot manifests,
  rather than stored twice. Identity comes from the verified provenance the
  pipeline already uses — S3 key, byte size, ETag, and a content hash where one
  has been computed. Filename-only identity is never sufficient and is refused.
* The **hot layer** is optimised for parsing throughput and may be deleted and
  rebuilt at any time.

Snapshot manifests are immutable: an existing manifest is never rewritten.

Availability is reported as one of `remote_available`, `hot`, `preserved`,
`materialised`, `verified`, or `unknown`. Preservation state is an
**operational** fact: it records how reproducible a run is and can never alter
the scientific snapshot identity.

Both roots are configuration (`storage.durable_snapshot_root`,
`storage.hot_cache_root`) and default to paths under `outputs/`; nothing
assumes a particular cluster path is writable.

> **No bulk data has been promoted into the durable store.** The architecture,
> provenance and configuration support exist; physically preserving existing
> snapshot data is a separate, explicitly approved operation.

---

## 9. Bronze / Silver / Gold

| Layer | Contents | Persisted | Regenerable |
|-------|----------|-----------|-------------|
| **Bronze** | Immutable source inventory for one snapshot: PDB ID, S3 key, compressed size, ETag, manifest timestamp. No scientific filtering. | yes | from the archive |
| **Silver** | The deterministic parsed representation. | **no, by design** | always, from Bronze identity |
| **Gold** | Everything scientifically derived: accepted/rejected chains, geometric validation, complete BRI, Brain, buckets, candidates, distances, classifications, the graph, the representative mapping, and the retained-chain release. | yes | by re-running the stages |

Silver is deliberately not persisted: it is reconstructed on demand from the
immutable Bronze object identity by the versioned parser, so the archive is
never stored twice. The planner reports it as `not_applicable` rather than
missing, and its gate is verified transitively — every Gold chain record
carries the source key and ETag it was parsed from.

**What is immutable:** published releases, the frozen 20260101 outputs, the two
frozen configuration YAMLs, and all historical provenance. **What gates
publication:** every stage must reach `validation_pass` before anything
downstream may start, and a release is only published once every gate has
passed.

---

## 10. Provenance and reproducibility

Every run gets its own directory, created **before** any work starts:

```
outputs/runs/run-<UTC stamp>-<8 hex>/
    run.json             the run record (atomic writes)
    events.jsonl         append-only event log
    resolved_run.yaml    the canonical configuration this run executes
    resolved_run.json    the same document, canonical JSON
```

Recorded: run ID and timestamps; snapshot selection mode **and** resolved
snapshot; the full resolved configuration with per-leaf value sources and both
hashes; git branch, commit and dirty state; Python and library versions; the
BRI implementation version; per-stage canonical identity, status, validation
verdict, input/output counts, paths, checksums and Slurm job IDs; runtime
environment (hostname, `$TMPDIR`, Slurm IDs) kept **separate** from the
canonical configuration; and the final release path with artefact hashes.

`run.json` is written atomically; `events.jsonl` is only ever appended to.
Historical provenance is never overwritten — creating a run refuses outright if
the directory already exists.

See [`docs/PROVENANCE.md`](docs/PROVENANCE.md) for the complete field list.

---

## 11. The web UI

`pdbclean ui` serves a restrained, information-dense research interface over
exactly the same backend the CLI drives. A UI-configured run and a CLI-configured
run produce the same `resolved_run.yaml` and execute the same commands; the
test suite asserts this against a live server.

| View | Contents |
|------|----------|
| **Run configuration** | Snapshot selector (with resolved identity, preservation status and hot-cache status); structural filtering (Q005 minimum backbone distance, minimum N–CA–C angle); **BRI representation precision p**; Brain filtering threshold; complete-BRI threshold τ; the fully resolved configuration with per-value sources; and the exact `resolved_run.yaml` a run would execute. |
| **Pipeline** | The canonical stage chain with status, validation verdict, purpose, parameters, counts, runtime, output paths, checksums and Slurm job IDs. |
| **Duplicate Explorer** | Filter by PDB ID, chain, classification, Stage-14 relationship, chain length and distance. It filters and displays; it never re-classifies, and its counts come from the stage summaries. |
| **Gold release** | Shown only once every gate has passed. Nothing is displayed for a run that has not completed. |
| **Runs** | Every run's provenance, with the full historical drill-down of [§12](#12-historical-run-workflow). |
| **Method** | The canonical pipeline and the Bronze/Silver/Gold lifecycle, generated from the same registry the planner uses. Each stage expands to its full scientific description. |

**Stage identity.** Every view labels stages canonically — `Stage 5 — Brain`,
never a bare `Brain` — and always in canonical pipeline order, never
alphabetically.

**Light / dark theme.** A toggle in the header switches between a professional
light theme and the existing dark theme. It respects the OS colour-scheme
preference when no explicit choice has been made, persists the choice locally,
switches without reloading or losing UI state, and is applied before first paint
so there is no flash. It is a viewing preference only: it never enters the
scientific configuration and cannot change either hash.

### 11.1 Artefact Viewer

Every supported artefact path shown anywhere in the UI — Gold release, Pipeline,
Runs, stage details, Duplicate Explorer — is clickable and opens the same
viewer. There is one implementation, not one per page.

| Format | View |
|--------|------|
| `.parquet` | schema, column types, row count, row-group count, writer, and the **actual rows**, paginated |
| `.json` | formatted, parsed |
| `.csv` / `.tsv` | paginated table |
| `.yaml`, `.txt`, `.log`, `_SUCCESS` | read-only text, bounded |
| anything else | metadata, path and hash only — never rendered |

The viewer header shows the artefact's provenance: run, stage, snapshot, both
configuration hashes, producing script, validation status, SHA256, size and row
count, so it is always clear which run and stage produced the table on screen.

**Paging is server-side and bounded.** A page is read from the file's row
groups and the read stops once the page is filled; the whole table is never
materialised. Page size is capped at 200 rows. Page 1 of the 3.24-million-row
classification table returns in tens of milliseconds. Search and sort operate
over a bounded scan window and the viewer says so rather than implying the
whole table was sorted.

**Download original** returns the recorded bytes unmodified — the same SHA256
as the artefact on disk, with a proper `Content-Disposition` filename. A
separate *Export view as CSV* is offered for the rows currently displayed; it is
labelled a convenience export, and the original Parquet remains the
authoritative scientific artefact.

**Path safety.** The viewer is not a filesystem browser. It reads only from the
configured output, release, run, durable-snapshot, hot-cache and Mol* report
roots. Absolute paths outside those roots, `../` traversal, symlinks escaping a
root, and repository source files are all refused with `403`.

### 11.2 Duplicate Explorer and Mol\*

Each row exposes two distinct things, deliberately kept apart:

* **Open source table** — the authoritative Parquet the row came from, opened in
  the Artefact Viewer. This is the scientific evidence.
* **View in Mol\*** — visual inspection of the pair. This is not evidence.

Scenes are generated **on demand** from the run's own snapshot and cached; no
pre-generated `.mvsj` is required. Available views are side-by-side, superposed
(Kabsch on paired backbone atoms, for display only), chains-only and deposited
context.

The representation follows the chain: many detected duplicates are only two or
three residues long, and a cartoon ribbon draws nothing for a dipeptide, so
chains shorter than 12 residues render as ball-and-stick and longer ones as
cartoon. The viewer reports its progress (resolving source, loading structures,
rendering) and, if a load fails, keeps the Mol\* host mounted and shows the
reason with a Retry button rather than collapsing to an empty panel. The cache is keyed by run, snapshot, pair and view, is safe to delete,
and is never a scientific artefact.

**Source availability is not Gold retention.** These are different questions
and the code keeps them apart:

| Question | Answer comes from |
|---|---|
| *Can we inspect this deposited structure?* | the run's **source** layer |
| *Was this chain kept in the deduplicated training dataset?* | Stage-14 retention |

Stage-14 removal means "not retained in the geometrically deduplicated Gold
training population". It does **not** mean the deposited structure was lost. A
removed chain is exactly as inspectable as a retained one; its retained/removed
status is shown as metadata and never gates visualisation.

**Snapshot correctness.** Structures resolve in this order — hot cache →
durable store → prepared examples → the run's own **Bronze source manifest**.
The manifest gives the snapshot-scoped object key and ETag, e.g.
`20260101/pub/pdb/data/structures/divided/mmCIF/ac/7acj.cif.gz`, so a missing
structure is materialised **on demand, one entry at a time**, and its ETag is
verified against the run's manifest before display. A mismatch is refused
rather than shown. Nothing undated is ever fetched, so a historical run cannot
silently display a newer revision.

All 1,072,751 near-duplicate pairs in the frozen release span 7,259 distinct
entries, **100% of which are present in that run's Bronze manifest** — so every
pair is source-resolvable without bulk-copying the snapshot.

**Chain namespaces.** PDBClean's canonical identity is `label_asym_id`, which
is what BRI was computed on and what the viewer selects. Deposited files also
carry `auth_asym_id`, and the two differ often — for 73.5% of removed chains in
the frozen release (`7acr` is label `Z`, auth `W`). Both are resolved from the
run's own cleaning output and displayed; the viewer never assumes they match.

Where a structure genuinely cannot be resolved the row gives the specific
reason — `source not materialised`, `source manifest missing`,
`not in snapshot`, `source fetch failed`, `chain not found` — never a bare
"no prepared scene". The pair page has an expandable **Source provenance**
panel showing run, snapshot, both chain identifiers, source key, ETag and the
local file, so the structure on screen is traceable to the run.

**Mol\* is inspection only and never determines duplicate classification.**
Complete-BRI L∞ remains authoritative; the viewer shows the recorded
classification, distance, representative relation, run and snapshot, and says so
on screen.

---

## 12. Historical run workflow

```
Runs  →  select a historical run
      →  ordered prerequisites and Stage 1–14 timeline
      →  click any stage
      →  identity, status, configuration, inputs, outputs, validation,
         execution provenance, reuse, and its artefacts
      →  optionally: Open in Duplicate Explorer  →  View in Mol*
```

Each stage panel shows, where the run actually recorded it: canonical stage
number and name, implementation and substage, Bronze/Silver/Gold layer and
scientific purpose; status and validation verdict; the scientific parameters
that stage used (including representation precision, thresholds, model scope
and snapshot identity); input counts and upstream stages; output counts, paths
and checksums; the validation gate; Slurm job IDs, node, timings, entry point,
git commit and both configuration hashes; and whether output was newly
generated or reused, with the reason.

Fields the run did not record display as **not recorded**. Nothing is invented.

Every stage also carries its full **scientific description** — purpose,
method, input, output, scientific role, implementation notes and method
references — shared across runs and shown alongside that run's actual values.
The descriptions are written to be detailed enough to explain the methodology
directly from the UI.

Artefacts are clickable and open in the Artefact Viewer ([§11.1](#111-artefact-viewer)):
Stage 10's classification table, Stage 14b's representative mapping, Stage 14c's
retained-chain dataset and removed-chain audit, and so on — each showing the
actual records, with the original file downloadable.

**Historical inspection is strictly read-only.** Opening a run or a stage never
modifies `run.json`, appends an event, modifies outputs or manifests,
re-resolves "latest", changes the snapshot identity, recomputes hashes,
launches a job, or changes any validation status.

Terminology: prerequisites are **not** scientific stage numbers; canonical
stage numbers are the scientific identities; execution ordinals are internal
orchestration order only; Stage 14a/b/c are subdivisions of Stage 14; and
Stages 3/4 and 8/9 share producers while keeping separate identities. See
[§3](#3-stage-to-code-map) for the canonical table.

---

## 13. Frozen COMP702 result [FROZEN]

Snapshot **2026-01-01**, protocol `protocol3.2-comp702-v1`, model 1,
p = 0.001 Å, τ = 0.010 Å inclusive.

### Population

| Quantity | Count |
|----------|-------|
| Canonical eligible chains | 578,524 |
| Chains with a defined Brain (*m* ≥ 2) | 577,760 |
| Chains with *m* = 1 (all retained) | 764 |
| Length buckets | 1,308 |

### Pairs

| Quantity | Count |
|----------|-------|
| Brain candidate pairs | 3,240,429 |
| **Tested pairs** | **3,531,895** |
| Near duplicates (`d ≤ 10 units`) | 1,072,751 |
| — of which exact (`d = 0`) | 17,373 |
| — of which non-zero near | 1,055,378 |
| Not near duplicates | 2,459,144 |

### Deduplication

| Quantity | Count |
|----------|-------|
| Graph edges | 1,068,256 |
| **Removed chains** | **78,754** |
| **Retained chains** | **499,770** |

> **A duplicate pair count is not a removed-chain count.** 1,072,751 near
> duplicate *pairs* led to 78,754 *chain* removals. One retained representative
> can absorb many pairs, and a pair between two retained chains removes
> nothing.

Every one of the 78,754 removals was independently audited to have its own
direct `d ≤ 10 units` edge to its assigned representative.

Release: `outputs/releases/PDBClean-20260101-protocol3.2-comp702-v1-dedup-v1`

| Artefact | SHA256 |
|----------|--------|
| retained dataset | `8ae52ad96586c2552f74083b480350973c86bdcca41ae1f30f7353472d769c8b` |
| removed-chain audit | `4cb3bea6c6a61f27de60818d097cf72c0c047f603d13f76c8286bbae647d3360` |
| release manifest | `1e6d6b249b6530fb501351fe6bd8d78647d3dad549db67e3fde486c2e3f8b918` |

#### Releases derived from it

Stage 14d publishes under its own identifier and never writes into the release
above, which remains byte-identical to its publication. Three exist:

| Release suffix | Parameters | Retained | Note |
|---|---|---|---|
| `-seqid100-v1` | identity 1.0, coverage 0.8 | 84,156 | **Not exact** — merges length variants; 518 removals below 50% identity, and `2OLO` merged into `2OLN`. |
| `-seqid100-cov100-v1` | identity 1.0, coverage 1.0 | 143,226 | Exact: every removal at 100% identity. 1,170 survivors still share a sequence (MMseqs2 non-transitivity). |
| `-seqid100-cov100-exact-v1` | identity 1.0, coverage 1.0, exact post-pass | **142,056** | Exactly one chain per distinct sequence. **This is the population currently training.** |

See [`docs/sequence_redundancy.md`](docs/sequence_redundancy.md) for the sweep
that motivates coverage 1.0 and the per-release validation gates.
| `_SUCCESS` | `945c6c34358b127ea07365384f6f50429af315a26d9878233e4638cacf34c400` |

Every artefact in this table is clickable in the UI's Gold release page and
opens in the Artefact Viewer — `retained_chains.parquet`,
`removed_chain_audit.parquet`, `representative_mapping.parquet` and the
near-duplicate edge tables all show their actual records, with the original file
downloadable.

**This release is immutable**, and the pipeline wrappers refuse to overwrite
it. These counts are the acceptance gates for *this snapshot only* — they are
never generic expectations, and a different snapshot derives its own counts.

Full detail: `docs/PDBCLEAN_2026_FINDINGS_AND_DECISIONS.md` and
`docs/provenance/pdbclean_20260101_dedup_v1.json`.

### 13.1 Threshold sensitivity [IMPLEMENTED]

[§7.3](#73-precision-p-versus-threshold-τ) poses τ and the Brain threshold as
research questions. Both have now been run at the frozen snapshot, holding
everything else fixed, so the sensitivity is measured rather than argued.

**τ = 0.005 Å against the frozen τ = 0.010 Å.** The Brain threshold moves with
τ, because the prefilter is only sound when it is at least as loose as the
classifier.

| Quantity | τ = 0.005 Å | τ = 0.010 Å (frozen) | Change |
|----------|------------:|---------------------:|-------:|
| Brain candidate pairs | 1,805,248 | 3,240,429 | −44.3% |
| Brain participating chains | 200,729 | 305,747 | −34.3% |
| Brain components | 43,715 | 56,517 | −22.7% |
| Near duplicates (total) | 931,644 | 1,072,751 | −13.2% |
| — graph edges (*m* ≥ 2) | 930,588 | 1,068,256 | −12.9% |
| — *m* = 1 pairs | 1,056 | 4,495 | −76.5% |
| Edge components | 20,249 | 20,789 | −2.6% |
| — cliques | 18,967 | 20,494 | −7.5% |
| — **non-cliques** | **1,282** | **295** | **×4.35** |
| Chains touching an edge | 96,714 | 99,854 | −3.1% |
| Chains with no edge | 481,046 | 477,906 | +0.7% |
| Representatives (*m* ≥ 2) | 21,303 | 21,100 | +1.0% |
| **Removed chains** | **75,411** | **78,754** | **−4.2%** |
| **Retained chains** | **503,113** | **499,770** | **+0.7%** |

Four things follow, and the third is the one that matters most.

**The edge sets nest exactly.** The τ = 0.005 Å edge set is the τ = 0.010 Å
edge set restricted to `d ≤ 5` units — not approximately, but bin for bin
across the whole distance histogram, 930,588 edges on both sides:

| `d` (mÅ) | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---:|---:|---:|---:|---:|---:|
| edges | 17,364 | 26,392 | 103,142 | 469,193 | 189,238 | 125,259 |

That is a monotonicity check on the whole search stack. An exact radius search
at a smaller radius *must* return a subset, and two independent runs of the
compressed cover tree over 578,524 chains agree to the single edge.

**The retained set barely moves.** Halving τ removes 3,343 fewer chains — 4.2%
of the removals, but only 0.7% of the retained population. The headline
deduplication result is not balanced on the choice of τ.

**Non-transitivity gets worse as τ tightens.** Non-clique components rise from
295 to 1,282. Fewer edges make components sparser, so more of them contain
chains that are *not* mutually near-duplicate. This is direct empirical support
for the Stage-14 rule in [§2.6](#26-redundancy-resolution-stage-14): a connected
component is not a duplicate equivalence class, and a design that removed
transitively would have removed chains with no direct edge to their
representative in 1,282 components at τ = 0.005 Å, against 295 at the frozen
threshold. Tightening the threshold makes transitive removal *more* wrong, not
less.

**The Brain prefilter is cheap only because it is tight.** Holding τ = 0.010 Å
and loosening the Brain threshold to 0.10 Å:

| Brain threshold | Candidate pairs | Relative cost |
|-----------------|----------------:|--------------:|
| 0.01 Å (frozen) | 3,240,429 | 1× |
| 0.10 Å | 70,651,044 | **21.8×** |

Both settings are *sound* — the prefilter is lossless at any threshold at least
as large as τ, so the final classification is identical either way. The 21.8×
is purely the cost of handing the cover tree more candidates to reject. This
run covers the prefilter stage only; the downstream stages were not executed,
because the answer they would produce is already known to be unchanged.

Reproducing either study is a configuration override, not a code change:

```bash
# tau study
pdbclean run --config config/pdbclean/protocol_3_2_comp702_v1.yaml \
  --set duplicate_search.near_duplicate_threshold_angstrom=0.005 \
  --set brain_filter.threshold_angstrom=0.005 \
  --executor slurm

# Brain prefilter cost study
pdbclean run --config config/pdbclean/protocol_3_2_comp702_v1.yaml \
  --set brain_filter.threshold_angstrom=0.100 \
  --executor slurm
```

τ must land exactly on the representation grid *p*, so `τ / p` is an integer —
see [§7.3](#73-precision-p-versus-threshold-τ). `0.005 / 0.001 = 5` is
accepted; a value like `0.0035` is rejected at configuration time rather than
silently rounded.

Outputs: `outputs/pdbclean_tau0p005/` and `outputs/pdbclean_brain0p10/`.

---

## 14. Testing

```bash
pytest tests -q
```

The suite includes a scientific regression layer that pins the method: the
validated defaults; inclusive `≤` semantics at, below and above the threshold;
Brain's dimension, first-row exclusion and undefinedness at *m* = 1; the Brain
prefilter against a brute-force oracle at several thresholds and lengths; the
frozen 20260101 summaries as oracles (population and pair accounting, the
direct-edge guarantee, components not treated as equivalence classes, all 764
*m* = 1 chains retained, Stage 13 not used as the Stage-14 edge set); both
frozen YAMLs asserted byte-identical by SHA256; precision and grid validation;
snapshot pinning and preservation; canonical stage ordering and identity;
historical-run read-only guarantees; and UI/CLI equivalence against a live
server.

The OpenFold training-view layer adds its own tests: the run encoding and its
round trip; the arithmetic placing Gold `label_seq_id` values on OpenFold's
seqres index, including a non-standard numbering origin; refusal — rather than
silent damage — when a lineage falls outside the deposited sequence, duplicates
a residue, or names an absent chain; and, for the alignment index, a known-good
index corrupted in each way it could plausibly be wrong (chain pointed at the
wrong MSA, truncated byte range, dangling reference, missing chain) with the
verifier required to fail on each. A verifier that cannot fail is not evidence.

Tests that read the large gitignored frozen outputs skip cleanly when absent.

A heavier regression re-runs the real Stage 14a/b entry points on a compute
node and compares every artefact against the frozen release:

```bash
sbatch task_scripts/run_stage14_regression.sbatch
```

It writes only to a scratch regression root; the frozen release is read-only.

---

## 15. OpenFold training view and retraining

This section covers the wiring between the frozen PDBClean Gold dataset and
OpenFold training. The deduplication science in §2–§13 is unchanged by it: what
follows only decides *which residues OpenFold reads* and *how it finds the
corresponding MSA*.

### 15.1 The retained-chain training view [IMPLEMENTED]

OpenFold builds every structural feature from `mmcif_object.chain_to_seqres`:
`make_mmcif_features` takes `num_res = len(chain_to_seqres[chain_id])` and
`get_atom_coords` then walks `range(num_res)`. That sequence is the *deposited*
polymer read from `_entity_poly_seq`, and it includes residues Protocol 3.2
removed.

The MSAs, however, were searched on each chain's `retained_sequence`. Feeding
OpenFold the deposited sequence therefore produces structure features of one
length and MSA features of another — and **the pipeline does not raise**. On
`102l_A` it silently produced `aatype (165, 21)` against `msa (1648, 163)`.
Training on that would learn against misaligned evolutionary signal, which is
worse than a crash.

The invariant now enforced for every training example is:

```
OpenFold projected sequence == Gold retained_sequence == MSA query row
```

character for character.

[`src/pdbclean/openfold_training_view.py`](src/pdbclean/openfold_training_view.py)
rebuilds one chain's `chain_to_seqres` and `seqres_to_structure` restricted to
Gold's `retained_label_seq_ids` and re-indexed to `0..L-1`. Coordinates, the
header, the model, the loss and the MSA semantics are untouched. The residue set
is **not** recomputed: it is taken from the same `retained_label_seq_ids` lineage
that `geometric_validation.reconstruct_retained_backbone_chain` uses, so the
project has one definition of "retained", not two.

OpenFold places Gold's residue ids on its own index as
`seq_idx = label_seq_id − min(_entity_poly_seq.num)`. That origin is **read from
the deposited file**, not assumed to be 1: an off-by-one there would shift every
residue against its coordinates while every length stayed correct, and nothing
downstream would notice.

### 15.2 Compact indexes [IMPLEMENTED]

Both indexes exist to avoid materialising hundreds of thousands of files against
a filesystem already near its inode quota.

**Alignment index** — `build_alignment_index.py`. OpenFold resolves an entry as
`open(join(alignment_dir, db))`, `seek(start)`, `read(size)`. Nothing requires
`db` to be a *packed* database, so each entry points at its own existing
`<sequence_sha256>.a3m` with `start=0`. That is the identical contract at **zero
duplicated bytes and zero new inodes**; packing would have copied 161 GiB.
Equivalence to a packed entry was confirmed against OpenFold's own
`_parse_msa_data` by reading the same MSA both ways.

**Projection index** — `build_projection_index.py`. Every one of the 499,770
lineages is a single contiguous run, so a chain's retained residues cost two
integers rather than a residue list (28 MB rather than hundreds).

Each projection entry also carries a truncated SHA256 of the retained sequence,
and the dataloader adapter verifies the projected sequence against it on every
sample. A length check is not sufficient — see §15.4. That digest is also the MSA
store's content address, so structure and MSA are tied together at load time.

### 15.3 Caches and the training population [IMPLEMENTED]

`build_chain_data_cache.py` drives OpenFold's train chain-data cache from the
Gold manifest rather than from a directory listing, because upstream's
`generate_chain_data_cache.py` records `chain_to_seqres` for every chain of every
entry. Using that here would have described a different population in two ways:
the sampler's length-based sampling probability and `max_single_aa_prop` filter
would be computed over residues Protocol 3.2 removed, and chains Gold did not
retain would be described. `seq` is therefore the retained sequence. Only
`resolution` and `release_date` come from the deposited file, read with
OpenFold's own `_get_header` so the values match what OpenFold would compute.

The cache is not optional: `OpenFoldDataset.looped_samples` indexes
`chain_data_cache[chain_id]` directly, so a missing cache is a `TypeError` on the
first batch, and a chain missing *from* it is dropped from training with only a
log line.

### 15.4 Eleven excluded chains [IMPLEMENTED]

Full-population validation found 11 chains whose projected sequence differs from
the retained sequence **at identical length**: `1aw8_B`, `1aw8_E`, `6rxh_B`,
`6v24_A`, `8au0_A`, `8au0_C`, `8tx9_A`, `8tx9_B`, `8tx9_D`, `9ixd_A`, `9ixf_A`.

The cause is upstream and is not a PDBClean defect. Entries with point
microheterogeneity carry two `_entity_poly_seq` rows for one residue number
(MET/MSE at 101 in `8tx9`, THR/AEI at 19 in `6v24`, PYR/SER at 1 in `1aw8`).
OpenFold's `_get_protein_chains` appends every row, so its seqres runs one
residue long and every index past that point shifts. **OpenFold cannot build
correct features for these chains with or without the projection.**

They are excluded rather than repaired: repairing would mean changing OpenFold's
core polymer parsing for every structure to fix a handful, and the alternative to
exclusion is training a chain against an MSA built for a different sequence. The
exclusion is recorded with its chain list in
`training_view/excluded_chains.json`, so the trained population is auditable and
never silently smaller than Gold.

**Gold population 499,770 → trainable population 499,759.**

> **Note.** The figures in this section describe the *geometry-only*
> population. The run now training uses the Stage-14d population instead —
> 142,056 chains, of which 138,056 train and 4,000 are held out. The held-out
> set is chosen so that no validation chain shares a complete-BRI near-duplicate
> edge *or* an identical sequence with any training chain, which the redundancy
> graph already knows and so costs a lookup rather than a search. Checkpoints
> are selected on `val/lddt_ca` every 5 epochs alongside the milestone
> schedule. Its projection lives in
> `outputs/openfold_training/20260101/pdbclean-final-v1/`.

### 15.5 Validation [IMPLEMENTED]

Heavy validation runs as Slurm CPU array jobs; none of it runs on a login node.

| Check | Result |
|-------|--------|
| MSA corpus complete, unduplicated, correctly attributed | 142,056 / 142,056, 0 missing, 0 orphan, 0 misattributed |
| Alignment index resolves to each chain's own MSA | 499,770 chains, 0 uncovered, 2,000 read back |
| Projected sequence == retained sequence (full population) | 499,759 / 499,770 |
| Terminal-trimmed chains covered | 8,841 |
| `auth_chain_id != label_chain_id` chains covered | 249,033 |
| Structure/MSA feature widths agree | 1,600 / 1,600 |
| Dataloader preflight through the real feature pipeline | 12 / 12, dims agree after cropping |

Counting entries proves nothing on its own here: an off-by-one in the chain
mapping still gives every chain *an* MSA. The population check therefore follows
each entry exactly as OpenFold will and compares the delivered sequence to the
chain's own Gold sequence.

### 15.6 The OpenFold input-view adapter

OpenFold itself is **not** vendored into this repository. A working copy is used
at `/mnt/fastscratch/users/sgsrehm1/openfold_src`, taken from commit
`da89cd28446abcd7be95459b7b349dedaee666c0` of `saifrhman/702_BeyondAF`. It lives
on fastscratch because `$HOME` is at its inode quota.

The changes to it are confined to the input view and to operational limits:

* `data_modules.py` — a `projection_index` argument on `OpenFoldSingleDataset`,
  the projection and digest check in `_parse_mmcif`, and
  `train_projection_index_path` on `OpenFoldDataModule`. A chain with no
  projection entry raises rather than falling back to the deposited polymer.
* `train_openfold.py` — a `--train_projection_index_path` flag, and
  `OPENFOLD_SAVE_TOP_K` / `OPENFOLD_MILESTONE_EVERY_N_EPOCHS` so checkpoint
  retention can be bounded. Upstream keeps every epoch's checkpoint
  (`save_top_k=-1`); at ~1.5 GB each that would exceed the filesystem quota long
  before the step budget is reached.

No model architecture, loss, MSA semantics, BRI science, duplicate science or
retained-chain identity is modified.

### 15.7 Reported conflict: MSA sharing between sequence-identical chains

§16 records the frozen Stage-18 policy as *"an MSA is **not** shared between
sequence-identical retained chains."* **The implemented corpus does share
them.** 499,770 chains resolve to 142,056 distinct MSAs, so 357,714 chains reuse
an MSA generated for an identical sequence.

Per §18, this is reported rather than silently resolved.

The scientific argument for sharing is that an MSA is a deterministic function of
the query sequence and the reference database, so two identical retained
sequences searched separately yield the same alignment; sharing avoids 3.5× the
search compute and 161 GiB of duplicated storage for no change in result. That
argument was **not** ratified as a change to the frozen policy before the corpus
was built, so the policy text and the artefact disagree and the decision is
outstanding.

Nothing about the *provenance* is ambiguous: the store is content-addressed by
`sha256(retained_sequence)`, and every chain's mapping is verified.

### 15.8 Retraining status [IN PROGRESS]

A from-scratch training run is in progress on the trainable population. It uses
random initialisation — **no pretrained AlphaFold or OpenFold weights are
loaded** — because the scientific goal is a new model trained on the
geometry-cleaned dataset.

Configuration: AF2 `initial_training` preset, crop 256, 128 MSA clusters, bf16,
1×A100, per-GPU batch 1 with gradient accumulation 8 (**global batch 8**), Adam
with the AlphaFold LR schedule (1,000-step warmup, decay from 50,000).
Continuation across the wall clock uses Lightning's SLURM requeue, which restores
global step, optimizer and scheduler, so the step counter never resets.

**This is not a reproduction of AlphaFold's training protocol**, and must not be
described as one:

* budget ~409,600 samples ≈ 51,200 optimizer steps — roughly **5%** of AF2's
  initial-training sample budget, and **less than one pass** over the 499,759
  chains;
* **no templates** — the alignment index carries only `.a3m`, so
  `make_template_features` takes its empty branch and `max_template_date` is
  inert;
* MSAs are UniRef30 via MMseqs2 only, not AF2's UniRef90 + BFD + Mgnify stack;
* no distillation set;
* global batch 8 versus AF2's 128;
* no fine-tuning stage and no held-out validation set.

**No trained model, accuracy, generalisation or failure-mode result is claimed.**

---

## 16. Current status

### Completed and frozen

* **Prerequisites A–C and Stages 1–14** — the PDBClean geometric
  deduplication pipeline, complete and frozen for snapshot 2026-01-01. See
  [§13](#13-frozen-comp702-result-frozen). **[FROZEN]**
* Configuration, provenance, orchestration, UI, Duplicate Explorer, Mol\*
  integration and historical-run inspection. **[IMPLEMENTED]**
* **Threshold sensitivity studies.** τ = 0.005 Å run to completion against the
  frozen τ = 0.010 Å, and the Brain prefilter cost measured at 0.10 Å. See
  [§13.1](#131-threshold-sensitivity-implemented). **[IMPLEMENTED]**

### In progress

These are roadmap numbers for the OpenFold preparation chain. They are distinct
from the *registered* pipeline stages in [§4](#4-repository-architecture): the
last registered stage is **Stage 14d**, sequence-redundancy resolution, which
defines the population this chain then prepares. See
[`docs/sequence_redundancy.md`](docs/sequence_redundancy.md).

The chain begins at 16, not 15. **Stage 15 is deliberately unassigned**, held
for sequence-redundancy resolution should it be promoted from a Stage-14
subdivision to a scientific stage in its own right. The gap is intentional and
not a missing entry.

* **Stage 16 — OpenFold training-population preparation.** **[IMPLEMENTED]**
  499,770 retained chains, 118,197 source entries.
* **Stage 17 — exact snapshot/source materialisation and retained-chain
  training-view preparation.** **[IMPLEMENTED]** All 118,197 mmCIFs
  materialised; the retained-chain training view is validated over the full
  population. See [§15](#15-openfold-training-view-and-retraining).
* **Stage 18 — fresh MMseqs2 alignment/MSA generation.** **[IMPLEMENTED]**
  142,056 MSAs, complete, unduplicated and correctly attributed.

  The frozen policy for this stage: MMseqs2; the query for each chain is that
  chain's exact PDBClean `retained_sequence`; a fresh MSA is generated for
  **every** one of the 499,770 retained chains; MSAs are **not** reused from
  RODA/OpenProteinSet, from COMP390, or from a previous snapshot; an MSA is
  **not** shared between sequence-identical retained chains. Reference
  databases are inputs, not downloaded MSAs. Note that
  `UniRef30_2021_03.tar.gz` is HH-suite format and must not be used as the
  MMseqs UniRef database.

  This work is out of scope for the pipeline productisation and was not
  modified by it.

* `scripts/openfold_training/` — prediction relaxation and BRI LAI work.
  **[IN PROGRESS]**

### Not yet done

* **Stage 19** — full structure ↔ alignment coverage validation.
  **[IMPLEMENTED]** Every retained chain resolves to its own MSA and to a
  projected structure of the same length; see
  [§15.5](#155-validation-implemented).
* **Stage 20** — OpenFold dataloader smoke test. **[IMPLEMENTED]** Real chains
  traverse the full feature pipeline with agreeing dimensions.
* **Stage 21** — GPU training smoke test. **[IMPLEMENTED]** Eight optimizer
  steps on real data with finite loss and gradients; checkpoint written and
  reloaded; peak 9.45 GiB.
* **Stage 22** — full OpenFold retraining. **[IN PROGRESS]** A from-scratch run
  is executing; see [§15.8](#158-retraining-status-in-progress) for its budget
  and for why it is not an AlphaFold-protocol reproduction.
* **Stage 23** — new checkpoint generation and downstream prediction /
  evaluation. **[FUTURE]**

**No training or evaluation results exist.** No trained model is published, and
no accuracy, generalisation or failure-mode result is claimed anywhere in this
repository. A training run being in progress is not a result.

---

## 17. Future work

Everything in this section is **[FUTURE]**. None of it is implemented, and none
of these questions has been answered.

**A. OpenFold retraining and evaluation.** Train OpenFold on the geometrically
deduplicated dataset, produce new checkpoints, and evaluate generalisation and
failure behaviour.

**B. Controlled comparison.** Where compute allows, a matched OpenFold training
run with and without geometric redundancy removal, to isolate the effect of
deduplication itself.

**C. Threshold studies.** *Partly answered* — see
[§13.1](#131-threshold-sensitivity-implemented). τ = 0.005 Å has been run
against the frozen τ = 0.010 Å, and the Brain prefilter cost measured at
0.10 Å. What remains is the *downstream* half of the question: whether a model
trained on the τ = 0.005 Å retained set behaves differently from one trained on
the frozen set. Since the two sets differ by only 3,343 chains out of ~500,000,
that comparison needs a design that can resolve a small effect, and it is not
worth GPU time until **A** and **B** are done.

**D. Representation-precision studies.** Vary p while retaining 0.001 Å as the
validated default, and study how stable redundancy relationships are under
coarser or finer representation. Note that the executable stages currently
implement only p = 0.001 Å; another grid requires changing the BRI
canonicalisation, which is a scientific decision.

**E. Sequence versus geometry.** Compare MMseqs sequence redundancy against BRI
geometric redundancy — how much do they agree, and where do they diverge?

**F. Structural novelty / out-of-distribution behaviour.** Evaluate model
behaviour as a function of geometric distance from, and density of, the
training data.

**G. Confidence and failure behaviour.** Investigate calibration, confidence
and failure modes on structures that are geometrically distant from training
data.

**H. Snapshot generalisation.** Run the pipeline reproducibly over additional
PDB snapshots and compare redundancy structure across releases.

**I. Snapshot storage optimisation.** Incremental and content-addressed
preservation at scale, hot-cache materialisation strategies, and storage
efficiency measurements.

**J. Training weighting.** A future architecture may combine geometric and
sequence information when weighting training examples. **No weighting function
is currently defined, frozen or implied**, and none should be inferred from
anything in this repository.

**K. Additional model backends.** AlphaFold or ColabFold integration if
feasible. OpenFold remains the current active training target.

---

## 18. Authority order

When sources disagree, resolve in this order:

1. current executable code and task scripts;
2. generated manifests, validation reports and run logs;
3. frozen Git commits and provenance;
4. the latest explicit project decisions;
5. README and documentation;
6. the COMP702 proposal;
7. scientific papers, including Wlodawer et al., *Acta Cryst D* 2025
   (doi 10.1107/S2059798325001883), in `reference/acta_2025/`;
8. the COMP390 dissertation material in `~/COMP702_BeyondAF/code/COMP390_code/`.

Where documentation conflicts with a frozen production artefact, the conflict
is reported rather than silently resolved. One such conflict — the superseded
historical `geometric_search` block — is recorded in
[`docs/CONFIGURATION.md`](docs/CONFIGURATION.md) §5.
