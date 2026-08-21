# Repository map

What lives where, and which parts are current, frozen, or historical.

Nothing in this repository has been deleted or moved by the productisation
work. Legacy material is labelled, not removed.

---

## Status legend

| Marker | Meaning |
|--------|---------|
| **CURRENT** | Active implementation; edit this. |
| **FROZEN** | Immutable published artefacts or provenance. Read only. |
| **LEGACY** | Historical COMP390 / PDB707K material, retained for the record. |
| **WIP** | Work in progress; see the README for its status. |

---

## Top level

| Path | Status | Contents |
|------|--------|----------|
| `src/pdbclean/` | CURRENT | The pipeline package: the scientific implementation plus the configuration, orchestration, provenance and UI layers. |
| `scripts/pdbclean/` | CURRENT | Per-stage entry points and their Slurm array wrappers for Stages 4–12. |
| `scripts/` (top level) | mixed | Stage 13–15 entry points (`build_stage14_*.py`, `select_stage14_representatives.py`) are CURRENT. `01_prepare_pdb707k.py`, `extract_pdb707k_bri_vectors.py`, `build_pdb707k_bri_manifest.py` and `analyze_2olo_*.py` are LEGACY PDB707K/COMP390 analysis. |
| `scripts/openfold_training/` | WIP | OpenFold training and relaxation scripts. Not part of the PDBClean pipeline. |
| `task_scripts/` | CURRENT | Slurm wrappers for Stages 13–15 and the scientific regression harness. All derive their arguments from the resolved run configuration. |
| `config/pdbclean/` | mixed | See below. |
| `tests/pdbclean/` | CURRENT | The full test suite, including the scientific regression tests. |
| `docs/` | CURRENT | This document, plus `ARCHITECTURE.md`, `CONFIGURATION.md`, `PROVENANCE.md`, the pipeline specification and the development status log. |
| `docs/provenance/` | FROZEN | Release provenance for the frozen 20260101 publication and the Acta review. |
| `outputs/pdbclean/20260101/protocol3.2-comp702-v1/` | FROZEN | The frozen Stage 1–14 outputs. Gitignored (large); read only. |
| `outputs/releases/PDBClean-20260101-…-dedup-v1/` | FROZEN | The frozen Gold release. Immutable. |
| `outputs/runs/` | CURRENT | Run provenance directories. Append-only. |
| `outputs/snapshot_store/` | CURRENT | Durable content-addressed snapshot preservation. |
| `outputs/snapshot_cache/` | CURRENT | Disposable hot materialisation. Safe to delete. |
| `reports/` | mixed | Acta review CSVs and evidence (FROZEN), plus `molstar_exact_duplicate_examples/` (CURRENT, wired into the UI). |
| `reproducibility/` | FROZEN | Pinned environment exports and `bri_version.txt`. |
| `tools/` | FROZEN | The pinned BRI v1.2.2 reference implementation used by the differential gate. |
| `reference/acta_2025/` | reference | Wlodawer et al., *Acta Cryst D* 2025 (doi 10.1107/S2059798325001883). |
| `data/` | LEGACY | `checked_PDB707K_cleaned_chains_sequences_19Feb2025.csv` — the PDB707K sequence table. |
| `code/COMP390_code/` | LEGACY | Minhao's COMP390 dissertation work. Retained in full; see `code/COMP390_code/LEGACY.md`. |
| `sbatch/` | LEGACY | COMP390-era PDB707K / AlphaFold / OpenFold batch scripts. See `sbatch/LEGACY.md`. |
| `logs/` | — | Slurm output. Not tracked as a deliverable. |

---

## `config/pdbclean/`

| Path | Status | Notes |
|------|--------|-------|
| `protocol_3_2_comp702_v1.yaml` | **FROZEN, byte-immutable** | Its SHA256 is embedded in frozen provenance. Contains a historical `geometric_search` block that is superseded — see `docs/CONFIGURATION.md` §5. |
| `stage14_representative_policy_v1.yaml` | **FROZEN, byte-immutable** | Its SHA256 is embedded in frozen provenance and in the release manifest. |
| `acta_downstream_investigation_v{1,2}.yaml` | FROZEN | Configuration for the Acta review passes. |
| `profiles/comp702_frozen_20260101.yaml` | CURRENT | Reproduces the frozen 20260101 run, including its expectation gates. |
| `pipeline_env.sh` | CURRENT | Infrastructure paths only. No scientific parameters. |

Because the two frozen YAMLs cannot change, **all** configuration added by the
productisation work lives in new files.

---

## `src/pdbclean/` in detail

### Scientific implementation (unchanged in method)

```
cleaning.py                        quality rules Q001–Q006
geometric_validation*.py           post-cleaning geometric gate
bri*.py, full_bri*.py              complete BRI
brain.py, brain_*.py               Brain (9-D average BRI)
brain_prefilter*.py                lossless same-length prefilter
compressed_cover_tree.py           exact L∞ radius search
duplicate_classification*.py       exact/near classification
downstream_metadata*.py            entry metadata
manifest.py, snapshot.py, gold.py  Bronze/Gold artefacts
```

### Productisation layer

```
defaults.py            the validated built-in defaults
runconfig.py           layered resolution, validation, hashing
snapshot_selection.py  discovery, picker, pinning
stage_registry.py      prerequisites + canonical Stage 1-14
pipeline.py            planning, reuse decisions, executors
run_provenance.py      run directories and the event log
duplicates.py          the Duplicate Explorer
run_inspection.py      read-only historical-run timeline and stage audit
artefacts.py           bounded, read-only artefact previews
snapshot_store.py      durable preservation + disposable hot cache
cli.py                 the `pdbclean` command
ui/                    the web UI (stdlib HTTP server + static assets)
```

---

## Entry points

```bash
pdbclean snapshots            # list available snapshots
pdbclean config               # show the resolved configuration and its sources
pdbclean plan                 # what would run, what would be reused
pdbclean run                  # execute (dry-run by default)
pdbclean status               # run provenance
pdbclean stages               # describe the stage chain
pdbclean stage-command        # print one stage's argv
pdbclean duplicates           # query the detected duplicate pairs
pdbclean ui                   # serve the web UI
```

Slurm:

```bash
sbatch task_scripts/run_stage14_geometric_graph.sbatch
sbatch task_scripts/run_stage14_representatives.sbatch
sbatch task_scripts/run_stage14_final_release.sbatch
sbatch task_scripts/run_stage14_regression.sbatch
```

Submit from the repository root, or export `PDBCLEAN_REPO_ROOT` — Slurm copies
batch scripts into its spool directory, so the repository cannot be inferred
from `BASH_SOURCE` inside a running job.
