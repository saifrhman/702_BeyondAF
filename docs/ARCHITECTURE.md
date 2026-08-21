# Architecture

How the PDBClean pipeline is put together. The scientific method it implements
is specified in `docs/pdbclean/pipeline_spec.md`; this document describes the
software around it.

---

## 1. One backend, three front ends

```
                 ┌──────────────────────────────────────────┐
   interactive   │                                          │
   CLI ─────────▶│                                          │
                 │   pdbclean.runconfig                     │
   non-inter-    │     resolve_run_config()                 │
   active CLI ──▶│         │                                │
   / config      │         ▼                                │
                 │   ResolvedRunConfig  ──▶ resolved_run.yaml
   Slurm  ──────▶│         │                                │
                 │         ▼                                │
   Web UI ──────▶│   pdbclean.pipeline                      │
                 │     plan_pipeline()  ──▶ PipelinePlan     │
                 │         │                                │
                 │         ▼                                │
                 │   pdbclean.cli.stage_command()           │
                 │         │                                │
                 │         ▼                                │
                 │   the EXISTING stage entry points        │
                 └──────────────────────────────────────────┘
```

There is exactly one implementation of the science. The CLI, the UI and the
Slurm wrappers are all thin layers that resolve a configuration and then invoke
the same stage entry points with arguments derived from it. No front end
contains a scientific decision.

`tests/pdbclean/test_ui_cli_equivalence.py` runs the real HTTP server and
compares its resolved configuration, its `resolved_run.yaml` and its plan
against the CLI's, so a divergence fails the build.

---

## 2. Modules

| Module | Responsibility |
|--------|----------------|
| `pdbclean.defaults` | The validated built-in defaults. Single source of truth for every scientific default. |
| `pdbclean.runconfig` | Layered resolution, per-leaf source tracking, validation, hashing, projection onto the Protocol 3.2 stage schema. |
| `pdbclean.snapshot_selection` | Discovery, the interactive picker, and pinning a run to a concrete snapshot identity. |
| `pdbclean.stage_registry` | Static description of every stage: canonical Stage 1-14 identity, layer, purpose, directory template, counts, dependencies, validation gate, entry point. |
| `pdbclean.pipeline` | Inspects artefacts, decides reuse/run/blocked, plans the run, and executes through a pluggable executor. |
| `pdbclean.run_provenance` | Owns a run directory: `run.json`, the append-only `events.jsonl`, and the resolved configuration. |
| `pdbclean.duplicates` | The read-only Duplicate Explorer over the published pair tables. |
| `pdbclean.run_inspection` | Read-only historical-run timeline and per-stage audit detail. |
| `pdbclean.artefacts` | Bounded, read-only previews of stage artefacts (JSON, CSV/TSV, Parquet, text). |
| `pdbclean.snapshot_store` | Durable content-addressed preservation and disposable hot materialisation. |
| `pdbclean.molstar_scenes` | The validated Mol* scene core, extracted verbatim from the preparation script. |
| `pdbclean.molstar_service` | On-demand, snapshot-correct scene generation with a disposable cache. |
| `pdbclean.source_index` | Bronze source-object and chain-namespace resolution for a run's snapshot. |
| `pdbclean.cli` | The `pdbclean` command. |
| `pdbclean.ui` | A stdlib `ThreadingHTTPServer` over exactly the modules above. |

Everything else under `src/pdbclean/` is the existing scientific implementation
(BRI, Brain, cleaning, geometric validation, the prefilter, the nearest-neighbour
search, classification). None of it was given a new scientific behaviour by the
productisation work.

---

## 3. Bronze / Silver / Gold

| Layer | Contents | Persisted |
|-------|----------|-----------|
| **Bronze** | Immutable source inventory for one snapshot: PDB ID, S3 key, compressed size, ETag, manifest timestamp. No scientific filtering. | yes |
| **Silver** | The deterministic parsed representation. | **no, by design** |
| **Gold** | Everything scientifically derived: accepted/rejected chains, geometric validation, complete BRI, Brain, length buckets, candidate pairs, complete-BRI distances, classifications, the redundancy graph, the representative mapping, and the final retained-chain release. | yes |

Silver is deliberately not persisted. It is reconstructed on demand from the
immutable Bronze object identity (S3 key, byte size, ETag) by the versioned
parser, so the archive is never stored twice. The planner reports the Silver
stage as `not_applicable` rather than pretending it is missing, and the gate is
verified transitively: every Gold chain record carries the source key and ETag
it was parsed from.

---

## 4. The registry: prerequisites + canonical Stage 1–14

The authoritative scientific vocabulary is **Stage 1 through Stage 14**, fixed
by `docs/PDBCLEAN_2026_FINDINGS_AND_DECISIONS.md` and by the
`pdbclean_stage<N>_*` schema names embedded in the frozen artefacts. The
immutable source inventory exists *before* Stage 1.

The registry has 15 entries: **three prerequisites plus the canonical stages.**
Its `ordinal` is an execution-order index only; the `canonical_stage` label is
the authoritative identity and is what the CLI, the UI and provenance all show.

| ordinal | canonical | registry entry | layer |
|---|---|---|---|
| 1 | prerequisite | `snapshot` | snapshot |
| 2 | prerequisite | `bronze_source_manifest` | bronze |
| 3 | prerequisite | `silver_parse` | silver |
| 4 | **Stage 1** | `structural_cleaning` | gold |
| 5 | **Stage 2** | `geometric_validation` | gold |
| 6 | **Stage 3–4** | `complete_bri` | gold |
| 7 | **Stage 5** | `brain` | gold |
| 8 | **Stage 6** | `length_buckets` | gold |
| 9 | **Stage 7** | `candidate_filtering` | gold |
| 10 | **Stage 8–9** | `complete_bri_nn` | gold |
| 11 | **Stage 10** | `duplicate_classification` | gold |
| 12 | Stage 14 input | `downstream_metadata` | gold |
| 13 | **Stage 14a** | `redundancy_graph` | gold |
| 14 | **Stage 14b** | `representative_selection` | gold |
| 15 | **Stage 14c** | `gold_release` | gold |

Stage 3–4 and Stage 8–9 are single implementations: canonical Stage 4 is the
exact milliångström representation of the Stage 3 BRI, and canonical Stage 9 is
the distance representation produced by the Stage 8 search. Neither is a
separate executable step, and neither has been collapsed away — the labels are
carried explicitly.

### Stages 11, 12 and 13 are deliberately not orchestrated

Canonical Stages 11–13 were investigation and validation passes feeding the
Acta-style review, not steps on the release path:

| canonical | role | frozen output |
|---|---|---|
| Stage 11 | investigation | `acta_downstream_investigation_v2` |
| Stage 12 | validation | `acta_manual_review_manifest_v2` |
| Stage 13 | investigation | `acta_detailed_review_v2` |

They are declared in `stage_registry.NON_ORCHESTRATED_CANONICAL_STAGES` so the
registry never silently omits part of the canonical vocabulary, and no
orchestrated entry claims those labels.

**This matters.** The near-duplicate graph is the first step of canonical
Stage 14 — it is labelled `Stage 14a`, never `Stage 13`. Canonical Stage 13 is
the manual detailed-review subset and is explicitly **not** the global Stage-14
deletion relation. `test_the_graph_is_stage_14_not_stage_13` pins this.

`pdbclean stages` prints each entry's canonical label, purpose and validation
gate; the registry is the single source for the CLI, the planner, the UI's
Method page and the per-stage provenance record.

The scientific sequence these stages implement is:

> clean → exact chain-length grouping → complete BRI → Brain → Brain filtering
> → fast NN on complete BRI → exact complete-BRI L∞ classification →
> redundancy resolution

Brain is the filtering and indexing layer only. The final classification is
always complete-BRI L∞.

---

## 5. Stage lifecycle and validation gates

A stage moves through an explicit vocabulary that the UI renders distinctly:

```
pending → running → execution_complete → validating → validation_pass → complete
                                              │
                                              ├─▶ validation_fail
                                              └─▶ partial          (output, no marker)

blocked          (an upstream gate has not passed)
not_applicable   (Silver)
```

A stage may only start once every stage it depends on has reached
`validation_pass`. Execution finishing is not the same event as validation
passing, and the two are never collapsed.

---

## 6. Restartability

Reuse is decided from **configuration and stage summaries**, never from
directory names.

For each stage, `inspect_stage`:

1. locates the stage directory from the resolved configuration;
2. reads its `global_summary.json`;
3. compares the summary's recorded parameters against the resolved
   configuration through the stage's `compatibility` map;
4. checks the success marker and the primary output.

If the summary disagrees with the configuration — a different snapshot, a
different protocol, a different threshold — the stage is reported as
`execution_complete` / `validation_fail` and its action is `run`. Output
produced under a different scientific configuration is never silently reused,
even when it sits in a correctly named directory.

If the summary agrees and the marker and output are present, the action is
`reuse`. Everything downstream of a non-passing stage is `blocked`.

---

## 7. Execution

`plan_pipeline` decides *what* to do; an `Executor` decides *how*:

| Executor | Behaviour |
|----------|-----------|
| `DryRunExecutor` | the default: prints the argv, executes nothing |
| `LocalExecutor` | runs the command in-process on this host |
| `SlurmExecutor` | `sbatch --parsable`, records the job id in provenance |

Barkla login nodes are for lightweight inspection and orchestration only.
Anything heavy — parsing, BRI, the searches, validation — runs on compute nodes
via `sbatch`. Batch scripts are never executed directly.

---

## 8. Provenance

Every run gets its own directory under `storage.run_root`, created **before**
any work starts. See `docs/PROVENANCE.md`.

---

## 8a. Snapshot durability

Storage separates a **durable, content-addressed** layer that keeps a completed
run reproducible after any cache expires, from a **disposable hot cache**
optimised for parsing throughput. An object unchanged between two snapshots is
preserved once and referenced by both snapshot manifests; identity comes from
the verified S3 key, size, ETag and content hash, never from a filename.
Snapshot manifests are immutable. Both roots are configuration
(`storage.durable_snapshot_root`, `storage.hot_cache_root`).

Availability is reported as `remote_available`, `hot`, `preserved`,
`materialised`, `verified` or `unknown`. It is an operational fact and can
never alter the scientific snapshot identity. See README §8.2.

---

## 8b. Historical-run inspection

`pdbclean.run_inspection` assembles a run's own records into the canonical
scientific timeline and, per canonical stage, the identity, status,
configuration, inputs, outputs, validation, execution provenance, reuse and
artefacts that run recorded. Fields a run did not record read as
`not recorded`; nothing is invented.

It is strictly read-only: no write, no event append, no snapshot
re-resolution, no hash recomputation, no job submission. Artefact previews are
bounded by row count and byte count, so a large dataset is never materialised
into a browser.

---

## 9. The Duplicate Explorer and Mol\*

The Explorer is a paginated, filterable, read-only view over the published pair
tables. It prefers the Stage-11 classifications and falls back to the Stage-10
near-duplicate tables, joining the Stage-14 representative mapping to report
each pair's relationship (`removed` / `retained` / `unaffected`). It filters and
displays; it never re-classifies, and the counts it reports come from the stage
summaries the pipeline already wrote.

Mol\* is wired in for **human inspection only**, and now works on demand for
any detected pair whose structures are locally available, rather than only for
pre-generated example scenes.

`pdbclean.molstar_scenes` holds the scene core — mmCIF backbone extraction,
Kabsch superposition, MolViewSpec builders — extracted verbatim in behaviour
from `reports/molstar_exact_duplicate_examples/prepare_scenes.py`, which is a
script that executes at import. A test asserts the extracted core still
reproduces the frozen `.mvsj` transform matrices and `metrics.json` figures, so
the two cannot drift.

**Source availability is separate from Gold retention.** Stage-14 removal
excludes a chain from the deduplicated training population; it does not delete
the deposited structure. Resolution therefore goes through the run's *source*
layer and never consults `retained_chains.parquet`.

`pdbclean.source_index` supplies that layer from the run's own records: the
Bronze manifest gives each entry's snapshot-scoped object key, ETag and size,
and the accepted-chain table gives both chain namespaces
(`label_asym_id` / `auth_asym_id`, which differ for most removed chains). Both
use Arrow filter pushdown, so one lookup does not load a 246k-row manifest.

`pdbclean.molstar_service` resolves each structure from *this run's* snapshot
(hot cache, then durable store, then prepared examples, then materialised on
demand from the Bronze source object with its ETag verified — never an undated
"current entry" fetch),
builds the requested view, and caches the result under a key covering run,
snapshot, pair and view. The cache is derived data: deleting it loses nothing,
and generating a scene changes no configuration hash and touches no release.

A rendered view never determines whether two chains are duplicates — that comes
from the complete-BRI calculation alone, and both the service payload and the
viewer say so.

---

## 9a. Artefact inspection

`pdbclean.artefacts` provides one bounded, read-only access layer used by every
UI surface that shows an artefact: metadata and hashing, format detection,
Parquet schema, and paginated table access that reads row groups until the
requested page is filled rather than materialising the file. Page size is
capped; search and sort run over a bounded scan window and report that they did.

The server exposes it behind a strict allowlist covering only the configured
output, release, run, durable-snapshot, hot-cache and Mol* report roots.
`Path.resolve()` collapses traversal and follows symlinks *before* the
containment check, so neither can escape. Repository source files are not
readable. Downloads return the recorded bytes unmodified.
