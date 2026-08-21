# Configuration

Every scientific parameter of the PDBClean pipeline is explicit, resolvable
from configuration, and recorded in provenance. Nothing scientific is decided
by a hard-coded literal in a script or by the name of a directory.

This document describes how a value is chosen. It does not describe the
scientific method — for that, see `docs/pdbclean/pipeline_spec.md`.

---

## 1. Precedence

Values resolve in one fixed order, weakest first:

| # | Layer | Where it comes from |
|---|-------|---------------------|
| 1 | `builtin_default` | `pdbclean.defaults.VALIDATED_DEFAULTS` |
| 2 | `config_file:<path>` | `--config <file>`, or the UI profile selector |
| 3 | `override:<origin>` | `--set key=value` (CLI), or a UI form field |

A later layer replaces a value supplied by an earlier one. Mappings merge
recursively; scalars and lists are replaced wholesale. A profile may name a
`base_config:`, which is applied *before* the profile itself, so the profile
remains the stronger of the two.

The resolver records which layer supplied **every leaf**, so provenance can
answer "where did this value come from?" for any parameter:

```
$ pdbclean config --config config/pdbclean/profiles/comp702_frozen_20260101.yaml
  Snapshot identity        20260101   [config_file:...comp702_frozen_20260101.yaml]
  Brain dimension          9
  ...
```

A value shown without a bracketed origin came from the validated built-in
defaults.

### The resolved configuration is what runs

The input file is not what executes. The **fully resolved** configuration is:

* projected onto the Protocol 3.2 stage schema (`to_protocol_config`);
* written to the run directory as `resolved_run.yaml` and `resolved_run.json`;
* hashed, and both hashes recorded in `run.json`.

The CLI, the UI and the Slurm wrappers all build stage arguments from the same
resolved configuration through `pdbclean.cli.stage_command`, so they cannot
drift apart.

---

## 2. Canonical configuration vs runtime provenance

A run resolves its configuration **once**, at creation, and persists it as
`resolved_run.yaml` / `resolved_run.json`. Downstream work — a later stage, a
Slurm job on another node, the UI inspecting the run — *loads* that document
(`load_resolved_config`, `RunProvenance.resolved_config`) rather than
re-resolving it. A run's identity therefore cannot drift with the host that
happens to execute it.

Host facts are kept out of the canonical document by design.
`storage.temporary_root` is stored as the template `${TMPDIR}/pdbclean` — a
*policy* ("use this node's own scratch") that reads the same on every node.
Variables in `RUNTIME_ENVIRONMENT_VARIABLES` (`TMPDIR`, `TMP`, `TEMP`,
`HOSTNAME`, `SLURMD_NODENAME`, the `SLURM_*` identifiers) are deliberately left
unexpanded; every other `${VAR}` expands as normal.

The template is resolved to a concrete path in exactly two places, both of
which are execution, not configuration:

* `to_protocol_config()` — the handoff that actually launches a stage;
* `runtime_environment()` — recorded into provenance as what the job ran with.

So the run record separates the two cleanly:

```
resolved_config      storage.temporary_root: "${TMPDIR}/pdbclean"   (canonical)
runtime.created_on   TMPDIR=/scratch/node010/999999, hostname, SLURM_JOB_ID
runtime.observed[]   one entry per stage execution, with the host it ran on
```

### Two hashes

| Hash | Covers | Use |
|------|--------|-----|
| `resolved_config_sha256` | the complete canonical document | exact run identity |
| `scientific_config_sha256` | the scientific projection only | *scientific* run identity |

Both are host-independent. Verified on Barkla: the frozen profile resolves to
`60482c15…` / `25b8e62a…` on a login node and on compute node `node008` alike.

The scientific projection contains `release`, `selection`, `quality_rules`,
`post_cleaning_geometric_validation`, `bri`, `brain`, `brain_filter`,
`duplicate_search`, `graph`, `representative_selection`, the snapshot
**identity** (`snapshot_id`, `bucket_url`) and `defaults_version`.

It excludes `storage`, `execution`, `observability` and `expectations`, so that
choosing `--executor slurm` instead of `local`, or gating a run on known
dataset counts, does not make it look like a different experiment — while a
change to a threshold, a rule, the model scope, the snapshot or the
representative policy always does.

Two runs with the same `scientific_config_sha256` are the same experiment.

---

## 3. What is configurable, and what is not

### Configurable without editing source

| Key | Default | Meaning |
|-----|---------|---------|
| `snapshot.mode` | `latest_complete` | how the snapshot is chosen |
| `snapshot.snapshot_id` | `null` | explicit `YYYYMMDD` identity |
| `selection.models.model_id` | `1` | model scope |
| `quality_rules.backbone_distance.minimum_distance_angstrom` | `0.01` | Q005 |
| `post_cleaning_geometric_validation.minimum_triangle_angle_degrees` | `3.0` | N–CA–C gate |
| `bri.representation_precision_angstrom` | `0.001` | representation precision **p** (see §3.1) |
| `brain_filter.threshold_angstrom` | `0.010` | Brain prefilter radius |
| `duplicate_search.near_duplicate_threshold_angstrom` | `0.010` | near-duplicate criterion |
| `representative_selection.minimum_deduplicated_chain_length` | `2` | shortest deduplicated chain |
| `storage.*`, `execution.*`, `observability.*` | see defaults | infrastructure only |

Changing any of the scientific rows produces a different
`scientific_config_sha256`, so existing stage output is **not** reused and the
downstream stages are re-planned. That is the intended behaviour: a different
configuration is a different experiment, not a variant of the same one.

### Fixed by definition — refused as overrides

| Key | Fixed at | Why |
|-----|----------|-----|
| `brain.dimension` | `9` | Definition 5.1 |
| `duplicate_search.metric` | `L_infinity` | the method's metric |
| `duplicate_search.final_classification_basis` | `complete_BRI` | Brain filters, complete BRI classifies |
| `duplicate_search.operator` | `less_than_or_equal` | the criterion is inclusive |
| `brain_filter.operator` | `less_than_or_equal` | as above |

`validate_resolved_config` rejects a configuration that changes any of these,
in the CLI and in the UI alike, before a run identity is created.

Two further guards are enforced:

* **Grid.** A threshold must land exactly on the *configured* precision grid
  (see §3.1). At the validated p = 0.001 Å, `0.0105 Å` is refused.
* **Losslessness.** `brain_filter.threshold_angstrom` may not be *smaller* than
  `duplicate_search.near_duplicate_threshold_angstrom`. The Brain filter is a
  lossless prefilter; a smaller radius would silently discard true complete-BRI
  near duplicates.

### 3.1 Representation precision p versus threshold tau

`bri.representation_precision_angstrom` (**p**, default `0.001`) is the grid on
which complete BRI is *represented*.
`duplicate_search.near_duplicate_threshold_angstrom` (**tau**, default `0.010`)
is what a distance is *compared to*. They are different numbers with different
meanings and are never conflated.

They are two independent experimental axes:

| Axis | The question it asks |
|------|----------------------|
| **p** | How sensitive is BRI-based redundancy detection to the numerical precision at which backbone geometry is represented? |
| **tau** | How sensitive are redundancy relationships, and downstream model behaviour, to the geometric near-duplicate threshold? |

Both are configurable; each independently changes the scientific identity.

**Grid compatibility.** A threshold must be an exact whole number of
representation units -- `tau / p` must be an integer:

| p | tau | units | |
|---|-----|-------|---|
| 0.001 | 0.010 | 10 | validated default; 1 unit = 1 mA |
| 0.002 | 0.010 | 5 | accepted |
| 0.005 | 0.010 | 2 | accepted |
| 0.003 | 0.010 | -- | **rejected**, never silently rounded |

`quantise_angstrom_to_units` performs this check deterministically and raises
with the offending ratio rather than rounding.

**Unit terminology.** Only the validated grid has an established physical name.
At p = 0.001 A one unit is one milliangstrom and the pipeline says "mA"; on any
other grid it says "representation unit". Historical frozen fields such as
`d_bri_mA` are correct for the runs that produced them and are never renamed.

**What p does not change.** Exposing p as configuration did not change the
executable method. p = 0.001 A remains *structural* in the production code:
`bri.compute_bri` ends with `numpy.around(..., 3)` to reproduce the pinned
BRI v1.2.2 canonicalisation, `full_bri_compare` scales by 1000,
`brain.compute_brain` requires canonical 3-decimal input, and `full_bri_nn`
requires exact int64 milliangstroms.

A run configured for any other p is therefore accepted as a **distinct
scientific configuration** -- it gets its own `scientific_config_sha256`, and
`PRECISION_DEPENDENT_STAGES` prevents it from reusing any artefact from Stage 3
onward -- but the production stages refuse to execute it. Implementing another
grid means changing the BRI canonicalisation, which is a scientific decision,
not a configuration change.

### Dataset-version expectation gates

`expectations.*` are assertions about a *particular dataset*, not parameters of
the method. The COMP702 figures — 578,524 canonical, 577,760 with a Brain,
1,068,256 edges, 499,770 retained, 78,754 removed — are the known-correct
counts for the 2026-01-01 experiment. They are **regression gates for that
snapshot only** and are never generic defaults:

| | frozen 2026-01-01 profile | a new / latest snapshot |
|---|---|---|
| `expectations.*` | the nine known counts | all `null` |
| Stage 13/15 argv | `--expected-…` gates | `--no-expectation-gate` |
| m ≥ 2 population | asserted against the gate | derived from *its own* Stage-6 `brain_defined_chain_count` |
| Failure mode | drift from 2026 fails the run | counts differ from 2026 freely |

The built-in defaults declare every expectation as `null`, and
`test_validated_defaults_carry_no_frozen_counts` asserts that no 2026 figure
appears anywhere in them. A new snapshot's counts are not known in advance, so
requiring them to match 2026 would be wrong.

A gate is never dropped *silently*: a stage invoked with no expectations must
say so explicitly with `--no-expectation-gate`, otherwise it refuses to run.
Where an expectation is supplied it is asserted against the value observed from
the data — a gate checks a result, it never substitutes for measuring one.

---

## 4. Snapshot selection

The default is **the latest complete snapshot**. A historical snapshot stays
reproducible by naming it:

```bash
pdbclean plan                                     # latest complete
pdbclean plan --snapshot 2026-01-01               # explicit, ISO form
pdbclean plan --snapshot 20260101                 # explicit, archive form
pdbclean plan --config config/pdbclean/profiles/comp702_frozen_20260101.yaml
pdbclean snapshots                                # list, newest first
```

Both `YYYY-MM-DD` and `YYYYMMDD` are accepted everywhere and normalise to the
archive identity. An unquoted `snapshot_id: 20260101` in YAML parses as an
integer; it is coerced back to the eight-digit string, because it names the
same snapshot.

Precedence for the snapshot specifically: an explicit `--snapshot` (or the UI
picker) beats `snapshot.snapshot_id` in the configuration file, which beats the
`latest_complete` default. If `snapshot.mode` requests an explicit snapshot but
no identity is set, the run is refused — it never silently falls back to
"latest".

Before any processing begins, the snapshot is resolved to a **concrete,
verified identity** and pinned into the configuration, which is what gets
hashed and written to provenance. The interactive picker offers position 1 as
"latest complete snapshot [default]", so pressing Enter takes the default.

---

## 5. The historical `geometric_search` block

`config/pdbclean/protocol_3_2_comp702_v1.yaml` contains a `geometric_search`
block with `query_radius: 1.0` and `retain_when: less_than 1.0`. `docs/pdbclean/pipeline_spec.md` §9–10 and `docs/pdbclean/development_status.md` §24 describe the same design.

**That is the historical PDB707K design, and it is superseded.** The executable
2026 pipeline, every generated manifest, and the frozen release all use the
inclusive complete-BRI criterion `d_bri_mA <= 10` (0.010 Å). Under the project's
authority order, executable code and generated manifests outrank documentation
and legacy configuration.

`geometric_search` is therefore **deliberately absent** from
`VALIDATED_DEFAULTS`: it is not applied, and it is not silently re-interpreted.
The two frozen YAML files are byte-immutable — their SHA256s are embedded in
frozen provenance — so all new configuration lives in new files under
`config/pdbclean/profiles/`.

---

## 6. Infrastructure

`config/pdbclean/pipeline_env.sh` derives every path from the repository
location and allows every value to be overridden from the environment. It
contains no scientific parameters.

| Variable | Default |
|----------|---------|
| `PDBCLEAN_REPO_ROOT` | derived from the file's own location |
| `PDBCLEAN_PYTHON` | `python` |
| `PDBCLEAN_CONDA_ENV` | `$HOME/fastscratch/envs/bri_env_1.2.2` |
| `PDBCLEAN_OUTPUT_ROOT` | `$PDBCLEAN_REPO_ROOT/outputs/pdbclean` |
| `PDBCLEAN_RELEASE_ROOT` | `$PDBCLEAN_REPO_ROOT/outputs/releases` |
| `PDBCLEAN_RUN_ROOT` | `$PDBCLEAN_REPO_ROOT/outputs/runs` |
| `storage.durable_snapshot_root` | `outputs/snapshot_store` (configuration, not env) |
| `storage.hot_cache_root` | `outputs/snapshot_cache` (configuration, not env) |
| `PDBCLEAN_LOG_ROOT` | `$HOME/fastscratch/pdbclean_logs` |
| `PDBCLEAN_PROTOCOL_CONFIG` | `config/pdbclean/protocol_3_2_comp702_v1.yaml` |
| `PDBCLEAN_REPRESENTATIVE_POLICY` | `config/pdbclean/stage14_representative_policy_v1.yaml` |
| `PDBCLEAN_PROFILE` | `config/pdbclean/profiles/comp702_frozen_20260101.yaml` |

**Slurm caveat.** Slurm copies a batch script into its spool directory, so
`BASH_SOURCE` cannot locate the repository from inside a running job. Task
scripts resolve `PDBCLEAN_REPO_ROOT` from the environment or from
`SLURM_SUBMIT_DIR`, with an explicit guard. Submit from the repository root, or
export `PDBCLEAN_REPO_ROOT`.
