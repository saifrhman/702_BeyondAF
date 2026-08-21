# LEGACY — COMP390-era batch scripts

This directory holds the original PDB707K / AlphaFold / OpenFold batch scripts
from the COMP390 era:

| Script | Purpose |
|--------|---------|
| `01_prepare_pdb707k.sbatch` | prepare the PDB707K dataset |
| `02_run_2olo_alphafold.sbatch` | AlphaFold run for the 2OLO case study |
| `03_run_2olo_openfold_and_project.sbatch` | OpenFold run and projection |
| `smoke_2olo_trained_ep32.sbatch` | smoke test of a trained checkpoint |

**Status: LEGACY.** They are retained as a historical record and are not part
of the current COMP702 PDBClean pipeline. They are not invoked by any current
entry point, and their paths and environment assumptions are not maintained.

They implement the earlier PDB707K design, in which geometric retention used
`L∞ < 1.0`. That design is **superseded** by the 2026 pipeline's inclusive
complete-BRI criterion `d_bri_mA <= 10` (0.010 Å) — see
`docs/CONFIGURATION.md` §5.

Current pipeline wrappers live in `task_scripts/` and `scripts/pdbclean/`, and
derive every argument from the resolved run configuration.
