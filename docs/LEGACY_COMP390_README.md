> **LEGACY DOCUMENT.** This is the previous repository README, describing the
> COMP390 script-cleanup work and its central path configuration. It is kept
> verbatim as a historical record. It does **not** describe the current COMP702
> PDBClean pipeline — for that, see the repository `README.md`,
> `docs/ARCHITECTURE.md` and `docs/REPOSITORY_MAP.md`.

---

# 702 BeyondAF

Cleaned and centralised codebase for the COMP702 dissertation project:

**Beyond AlphaFold: Filling the Blind Spots in Protein Structure Predictions**

This repository contains cleaned task scripts from the COMP390/OpenFold reproduction workflow, adapted for Saif Ur Rehman's COMP702 project on the University of Liverpool Barkla2 HPC cluster.

The main cleanup goal was to remove user-specific hard-coded paths, remove fixed node locks, replace legacy paths with a central config file, and make the scripts easier to reproduce.

## Project status

The current repository contains cleaned scripts for:

- OpenFold data preparation
- OpenFold training and checkpoint handling
- OpenFold inference
- Packing and archive utilities
- Relaxation runs
- pLDDT extraction and analysis
- Relaxation analysis
- Display-selection analysis

Final cleanup audit status:

    Live task script files checked: 48
    Non-ASCII text: none
    Minhao hard-coded paths: none
    Fixed gpu32/gpu15/gpu13/node037 locks: none
    Final audit: bad=0

## Main folders

The important folders in this repository are:

- `config/`
  - Contains the central path configuration file.

- `code/COMP390_code/task_scripts/`
  - Contains the cleaned task scripts, organised into Task 01 to Task 09 folders.

- `code/COMP390_code/task_scripts/Task_01_DataPrep_sbatch/`
  - Data preparation, RODA download, alignment extraction, cache generation, and OpenFold subset setup.

- `code/COMP390_code/task_scripts/Task_02_Training_sbatch/`
  - OpenFold training and optional data-staging scripts.

- `code/COMP390_code/task_scripts/Task_03_Inference_sbatch/`
  - OpenFold inference, repeated inference, checkpoint sweeps, and comparison scripts.

- `code/COMP390_code/task_scripts/Task_04_Packing_sbatch/`
  - Scripts for packing OpenFold subset data, test-set data, and alignment data.

- `code/COMP390_code/task_scripts/Task_05_LegacyMisc_sbatch/`
  - Legacy relaxation script retained and cleaned for reproducibility.

- `code/COMP390_code/task_scripts/Task_06_TrainingScripts/`
  - Python wrapper for OpenFold GPU training.

- `code/COMP390_code/task_scripts/Task_07_pLDDT_Analysis/`
  - pLDDT extraction and pLDDT/error relationship analysis.

- `code/COMP390_code/task_scripts/Task_08_Relaxation_Analysis/`
  - Relaxation success/failure analysis by file, epoch, repeat, and chain.

- `code/COMP390_code/task_scripts/Task_09_DisplaySelection/`
  - Display-selection analysis for choosing good chains and epochs to visualise.

## Central path config

All cleaned scripts use the central config file:

    config/comp702_paths.sh

Before running scripts manually, source the config:

    source ~/COMP702_BeyondAF/config/comp702_paths.sh

Important variables include:

    COMP702_ROOT
    COMP390_ROOT
    OPENFOLD_CODE_DIR
    OPENFOLD_ENV
    COMP702_DATA_ENV
    OPENFOLD_DATA_ROOT
    OPENFOLD_RUNS_ROOT
    OPENFOLD_INFERENCE_FASTA_DIR
    OPENFOLD_INFERENCE_CIF_DIR
    OPENFOLD_INFERENCE_ALIGN_DIR
    OPENFOLD_INFERENCE_OUT_ROOT
    OPENFOLD_INFERENCE_RUNS_ROOT
    OPENFOLD_EP32_CKPT

## What is intentionally not included

The following are intentionally excluded from Git:

    data/
    outputs/
    logs/
    envs/
    code/COMP390_code/openfold/
    code/COMP390_code/checkpoints/
    code/COMP390_code/Result/
    *.ckpt
    *.pt
    *.pth
    *.npz
    *.pkl
    *.pickle
    *.tar
    *.tar.gz
    *.tar.zst
    *.zip

This keeps the repository lightweight and avoids uploading large datasets, model checkpoints, generated outputs, logs, and environment folders.

## Notes on filenames

Some filenames still contain historical labels such as `gpu32` or `node037`.

Examples:

    openfold_subset_gpu32.sbatch
    extract_alignment_gpu32.sbatch
    pack_node037.sbatch
    run_4bp8_100pred_ep32_gpu32.sbatch

These names are historical. The cleaned script contents no longer force a specific node such as `gpu32`, `gpu15`, `gpu13`, or `node037`.

## Running scripts

Most `.sbatch` scripts are intended for Barkla2 Slurm usage:

    sbatch path/to/script.sbatch

Do not submit jobs blindly. Check the required inputs first, especially:

    echo "$OPENFOLD_INFERENCE_FASTA_DIR"
    echo "$OPENFOLD_INFERENCE_CIF_DIR"
    echo "$OPENFOLD_INFERENCE_ALIGN_DIR"
    echo "$OPENFOLD_EP32_CKPT"

For Python analysis scripts, activate the relevant environment first, then run:

    python path/to/script.py

## Re-running the cleanup audit

Use this from the project code root:

    cd ~/COMP702_BeyondAF/code/COMP390_code
    TASKS="$PWD/task_scripts"

    bad=0

    echo "===== Bash syntax check for all .sbatch files ====="
    while IFS= read -r f; do
      echo "Checking sbatch: ${f#$TASKS/}"
      bash -n "$f" || bad=1
    done < <(find "$TASKS" -type f -name "*.sbatch" | sort)

    echo
    echo "===== Python syntax check for all .py files ====="
    while IFS= read -r f; do
      echo "Checking python: ${f#$TASKS/}"
      python -m py_compile "$f" || bad=1
    done < <(find "$TASKS" -type f -name "*.py" | sort)

    echo
    echo "===== Non-ASCII check ====="
    if LC_ALL=C grep -RInP '[^\x00-\x7F]' "$TASKS" --include="*.sbatch" --include="*.py"; then
      bad=1
    else
      echo "No non-ASCII text found in live task scripts."
    fi

    echo
    echo "===== Minhao path / fixed node check ====="
    if grep -RIn "/users/sgmwu14\|sgmwu14\|/mnt/data2/users/sgmwu14\|/users/%u/fastscratch\|localscratch\|#SBATCH -w gpu32\|#SBATCH -w gpu15\|#SBATCH -w gpu13\|#SBATCH -w node037" "$TASKS" --include="*.sbatch" --include="*.py"; then
      bad=1
    else
      echo "No Minhao hard-coded paths or fixed node locks found in live task scripts."
    fi

    echo
    echo "bad=$bad"

Expected result:

    bad=0

## Current development notes

This repository is currently a cleaned reproduction base, not a fully runnable packaged release.

Before running inference or training scripts, the relevant OpenFold source code, model parameter files, checkpoints, FASTA files, mmCIF files, and alignment directories must exist on the HPC filesystem.

The next step is to run one small smoke test before launching larger training, inference, or sweep jobs.
