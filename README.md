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

## Contents

1. [Motivation](#1-motivation)
2. [The scientific pipeline](#2-the-scientific-pipeline)
3. [Stage-to-code map](#3-stage-to-code-map)
4. [Repository architecture](#4-repository-architecture)
5. [Installation](#5-installation)
6. [Running the pipeline](#6-running-the-pipeline)
7. [Configuration architecture](#7-configuration-architecture)
8. [Snapshot selection and preservation](#8-snapshot-selection-and-preservation)
9. [Bronze / Silver / Gold](#9-bronze--silver--gold)
10. [Provenance and reproducibility](#10-provenance-and-reproducibility)
11. [The web UI](#11-the-web-ui)
12. [Historical run workflow](#12-historical-run-workflow)
13. [Frozen COMP702 result](#13-frozen-comp702-result-frozen)
14. [Testing](#14-testing)
15. [Current status](#15-current-status)
16. [Future work](#16-future-work)
17. [Authority order](#17-authority-order)

Deeper detail lives in [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md),
[`docs/CONFIGURATION.md`](docs/CONFIGURATION.md),
[`docs/PROVENANCE.md`](docs/PROVENANCE.md) and
[`docs/REPOSITORY_MAP.md`](docs/REPOSITORY_MAP.md). This README is
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
OpenFold retraining — see [§15](#15-current-status) and [§16](#16-future-work)
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
```

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

### 2.7 What the result is, and is not

The output is **a geometrically deduplicated PDB chain dataset under COMP702
representative policy v1**.

It is *not* a claim that any removed chain corresponds to an incorrect or
invalid experiment.

---

## 3. Stage-to-code map

The canonical scientific vocabulary is **Prerequisites A–C, then Stage 1
through Stage 14**. Prerequisites are lettered so they can never be mistaken
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

**Shared producers.** Stages 3 and 4 share one implementation, because
`compute_bri` applies the precision grid at the point of computation. Stages 8
and 9 share one implementation, because the search emits its distances already
represented. Both canonical identities are preserved everywhere — in the
registry, the UI and provenance — because they are distinct scientific
concepts.

**Stage 14a/b/c** are *engineering subdivisions* of scientific Stage 14, not
separate scientific stages.

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
| `docs/` | DOCUMENTATION | Architecture, configuration, provenance, repository map, the pipeline specification and the development status log. |
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
| `data/` | LEGACY | The PDB707K sequence table. |
| `code/COMP390_code/` | LEGACY / HISTORICAL | Minhao's COMP390 dissertation work, retained in full. See its `LEGACY.md`. |
| `sbatch/` | LEGACY / HISTORICAL | COMP390-era batch scripts. See its `LEGACY.md`. |
| `logs/` | GENERATED | Slurm output. |

See [`docs/REPOSITORY_MAP.md`](docs/REPOSITORY_MAP.md) for the module-level
breakdown.

---

## 5. Installation

### Environment

The pipeline needs Python ≥ 3.10 with `numpy`, `scipy`, `pyarrow`, `pandas`,
`gemmi` and `PyYAML`. On Barkla the pinned environment is
`~/fastscratch/envs/bri_env_1.2.2`, exported in `reproducibility/`.

```bash
# editable install, provides the `pdbclean` entry point
pip install -e .

# with test dependencies
pip install -e ".[test]"
```

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
* **MMseqs2** for Stage 17 MSA generation — downstream, see
  [§15](#15-current-status).
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
| **Method** | The canonical pipeline and the Bronze/Silver/Gold lifecycle, generated from the same registry the planner uses. |

**Stage identity.** Every view labels stages canonically — `Stage 5 — Brain`,
never a bare `Brain` — and always in canonical pipeline order, never
alphabetically.

**Light / dark theme.** A toggle in the header switches between a professional
light theme and the existing dark theme. It respects the OS colour-scheme
preference when no explicit choice has been made, persists the choice locally,
switches without reloading or losing UI state, and is applied before first paint
so there is no flash. It is a viewing preference only: it never enters the
scientific configuration and cannot change either hash.

**Mol\*.** Prepared MolViewSpec scenes are integrated for visual inspection of
duplicate pairs. **Mol\* is inspection only and never determines duplicate
classification** — complete-BRI L∞ remains authoritative, and the viewer says
so on screen.

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

Artefacts are browsable in place: JSON is formatted, CSV/TSV is previewed as a
bounded table, Parquet shows schema, row count, key/value metadata and a
bounded row sample, text and logs open in a read-only viewer, and anything else
shows metadata, path and hash only. Previews are capped so a large dataset is
never loaded into the browser.

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
| `_SUCCESS` | `945c6c34358b127ea07365384f6f50429af315a26d9878233e4638cacf34c400` |

**This release is immutable**, and the pipeline wrappers refuse to overwrite
it. These counts are the acceptance gates for *this snapshot only* — they are
never generic expectations, and a different snapshot derives its own counts.

Full detail: `docs/PDBCLEAN_2026_FINDINGS_AND_DECISIONS.md` and
`docs/provenance/pdbclean_20260101_dedup_v1.json`.

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

Tests that read the large gitignored frozen outputs skip cleanly when absent.

A heavier regression re-runs the real Stage 14a/b entry points on a compute
node and compares every artefact against the frozen release:

```bash
sbatch task_scripts/run_stage14_regression.sbatch
```

It writes only to a scratch regression root; the frozen release is read-only.

---

## 15. Current status

### Completed and frozen

* **Prerequisites A–C and Stages 1–14** — the PDBClean geometric
  deduplication pipeline, complete and frozen for snapshot 2026-01-01. See
  [§13](#13-frozen-comp702-result-frozen). **[FROZEN]**
* Configuration, provenance, orchestration, UI, Duplicate Explorer, Mol\*
  integration and historical-run inspection. **[IMPLEMENTED]**

### In progress

* **Stage 15 — OpenFold training-population preparation.** **[IN PROGRESS]**
* **Stage 16 — exact snapshot/source materialisation and retained-chain
  training-view preparation.** **[IN PROGRESS]**
* **Stage 17 — fresh MMseqs2 alignment/MSA generation.** **[IN PROGRESS]**

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

* **Stage 18** — full structure ↔ alignment coverage validation. **[FUTURE]**
* **Stage 19** — OpenFold dataloader smoke test. **[FUTURE]**
* **Stage 20** — GPU training smoke test. **[FUTURE]**
* **Stage 21** — full OpenFold retraining. **[FUTURE]**
* **Stage 22** — new checkpoint generation and downstream prediction /
  evaluation. **[FUTURE]**

**No training or evaluation results exist.** No model has been trained on the
deduplicated dataset, and no accuracy, generalisation or failure-mode result is
claimed anywhere in this repository.

---

## 16. Future work

Everything in this section is **[FUTURE]**. None of it is implemented, and none
of these questions has been answered.

**A. OpenFold retraining and evaluation.** Train OpenFold on the geometrically
deduplicated dataset, produce new checkpoints, and evaluate generalisation and
failure behaviour.

**B. Controlled comparison.** Where compute allows, a matched OpenFold training
run with and without geometric redundancy removal, to isolate the effect of
deduplication itself.

**C. Threshold studies.** Vary the complete-BRI near-duplicate threshold τ and
observe the effect on the duplicate graph, the retained set and downstream
model behaviour.

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

## 17. Authority order

When sources disagree, resolve in this order:

1. current executable code and task scripts;
2. generated manifests, validation reports and run logs;
3. frozen Git commits and provenance;
4. the latest explicit project decisions;
5. README and documentation;
6. the COMP702 proposal;
7. scientific papers, including Wlodawer et al., *Acta Cryst D* 2025
   (doi 10.1107/S2059798325001883), in `reference/acta_2025/`;
8. the COMP390 dissertation material in `code/COMP390_code/`.

Where documentation conflicts with a frozen production artefact, the conflict
is reported rather than silently resolved. One such conflict — the superseded
historical `geometric_search` block — is recorded in
[`docs/CONFIGURATION.md`](docs/CONFIGURATION.md) §5.
