# Provenance

Every run is reconstructable from its own record. Provenance is written
*before* any work begins, it is append-only, and historical provenance is never
overwritten.

---

## 1. The run directory

```
outputs/runs/run-20260820T134501Z-3f9ab2c1/
    run.json                 # the run record, rewritten atomically
    events.jsonl             # append-only event log
    resolved_run.yaml        # the configuration this run executes
    resolved_run.json        # the same document, canonical JSON
```

The run id is `run-<UTC timestamp>-<8 hex>`. Creating a run refuses outright if
the directory already exists:

```
RunProvenanceError: Run directory already exists, refusing to overwrite
historical provenance: outputs/runs/run-...
```

`run.json` is written atomically (temp file + rename), so a crash mid-write
cannot leave a truncated record. `events.jsonl` is only ever appended to.

---

## 2. What `run.json` records

### Identity

| Field | Meaning |
|-------|---------|
| `run_id` | `run-<UTC stamp>-<8 hex>` |
| `created_at`, `updated_at`, `status_updated_at` | UTC ISO-8601 |
| `status` | `created` → … → `complete` / failed |
| `run_directory` | absolute path |

### Configuration

| Field | Meaning |
|-------|---------|
| `resolved_config` | the canonical document, resolved once at creation |
| `resolved_config_sha256` | hash of the whole document |
| `scientific_config_sha256` | hash of the scientific projection (see `docs/CONFIGURATION.md` §2) |
| `resolved_config_yaml`, `resolved_config_json` | paths to the persisted copies |
| `config_layers` | e.g. `["builtin_default", "config_file:…", "override:cli"]` |
| `config_file`, `config_overrides` | what the operator actually supplied |
| `config_value_sources` | **per-leaf** origin: `{"brain.dimension": "builtin_default", …}` |
| `defaults_version` | version of `VALIDATED_DEFAULTS` |

### Dataset

| Field | Meaning |
|-------|---------|
| `snapshot` | the pinned identity: `snapshot_id`, `display`, `layout`, `source_prefix`, `sample_mmcif_key`, `selection_mode` |
| `snapshot.selection_mode` | how the snapshot was chosen (`latest_complete`, `fixed`, ...). `latest_complete` resolves **once**; a resumed run never re-resolves it. |
| `snapshot_durability` | availability (`preserved` / `materialised` / ...) and whether the run is reproducible without a cache. Operational only; never part of the scientific identity. |
| `inputs` | per-input path, byte size and SHA256 |

### Code and environment

| Field | Meaning |
|-------|---------|
| `git.branch`, `git.commit` | HEAD at run time |
| `git.working_tree_dirty`, `git.dirty_path_count` | whether the tree was clean |
| `environment.python_version`, `python_executable`, `hostname`, `platform` | interpreter and host |
| `environment.numpy_version`, `scipy_version`, `pyarrow_version`, `gemmi_version`, `yaml_version` | scientific library versions |
| `environment.bri_implementation` | from `reproducibility/bri_version.txt` |
| `environment.slurm` | job id, array job/task id, partition, nodelist |

### Stages

One entry per stage:

| Field | Meaning |
|-------|---------|
| `stage_id`, `title`, `layer` | identity |
| `canonical_stage` | the authoritative Stage 1-14 label (`"Stage 7"`, `"Stage 14a"`, ...) |
| `status`, `validation` | the lifecycle states (see `docs/ARCHITECTURE.md` §5) |
| `started_at`, `finished_at`, `runtime_seconds` | timing |
| `input_count`, `output_count` | from the stage's own summary |
| `output_path`, `manifest_path`, `summary_path` | where the artefacts are |
| `checksums` | SHA256 of the primary outputs |
| `slurm_job_ids` | every job that contributed |
| `attempts`, `reused` | restart history |
| `scientific_parameters` | the resolved values this stage consumed |
| `messages` | human-readable notes, including any incompatibility found |

### Runtime (execution provenance)

Host facts live here, never in the configuration:

| Field | Meaning |
|-------|---------|
| `runtime.created_on` | the host the run was created on: hostname, `$TMPDIR`, Slurm identifiers, and the concrete storage paths its templates resolve to |
| `runtime.observed[]` | one entry per stage execution, recording the host that stage actually ran on |

This is what keeps the canonical configuration identical across machines: the
run stores `storage.temporary_root: "${TMPDIR}/pdbclean"` -- a policy that reads
the same everywhere -- while the node-local path it expanded to is recorded as
an execution fact. Both configuration hashes are therefore host-independent.

### Verdicts and release

| Field | Meaning |
|-------|---------|
| `validation` | named gates with `PASS` / `FAIL` and their evidence |
| `release` | `release_path`, and per-artefact `path` / `bytes` / `sha256` |

---

## 3. The event log

`events.jsonl` is one JSON object per line, each with `at` and `event`:

```jsonl
{"at":"…","event":"run_created","resolved_config_sha256":"…","scientific_config_sha256":"…","snapshot_id":"20260101"}
{"at":"…","event":"plan","stages":15,"to_run":0,"reusable":14}
{"at":"…","event":"stage_updated","stage_id":"complete_bri","status":"running"}
{"at":"…","event":"stage_updated","stage_id":"complete_bri","status":"execution_complete","slurm_job_ids":["10284723"]}
{"at":"…","event":"validation","name":"complete_bri","verdict":"PASS"}
{"at":"…","event":"run_status","status":"complete"}
```

Lines are only ever appended. To see what a run did and when, read the log in
order; to see the final state, read `run.json`.

---

## 4. Inspecting runs

```bash
pdbclean status                    # list every run under the run root
pdbclean status --run <run-id>     # one run, with its stage table
pdbclean status --run <run-id> --json
```

The UI's **Runs** page shows the same records, including both hashes and the
per-leaf configuration sources, and lets any stage artefact be opened in the
Artefact Viewer.

**Inspection is read-only.** Opening a run, a stage or an artefact never writes
`run.json`, appends an event, re-resolves a snapshot, recomputes a hash, changes
a validation status or launches a job. Generated Mol* scenes are written to a
disposable cache outside every release and run directory, and are not
provenance.

---

## 5. Frozen provenance

`docs/provenance/pdbclean_20260101_dedup_v1.json` records the frozen COMP702
20260101 release. It, and the release directory it describes, are immutable.

Two configuration files are byte-immutable because their SHA256s are embedded
in that frozen provenance:

* `config/pdbclean/protocol_3_2_comp702_v1.yaml`
* `config/pdbclean/stage14_representative_policy_v1.yaml`

`tests/pdbclean/test_scientific_invariants.py` asserts both hashes on every test
run. All configuration added by the productisation work lives in new files under
`config/pdbclean/profiles/`.
