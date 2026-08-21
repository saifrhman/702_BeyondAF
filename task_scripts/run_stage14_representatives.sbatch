#!/bin/bash
#SBATCH --job-name=stage14_reps
#SBATCH --partition=nodes
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=logs/stage14_reps_%j.out
#SBATCH --error=logs/stage14_reps_%j.err

# Stage 14: representative selection.
#
# Every path and every scientific argument is derived from the resolved run
# configuration, so this wrapper cannot drift from the CLI or the UI.
#
#   sbatch task_scripts/run_stage14_representatives.sbatch
#   PDBCLEAN_PROFILE=<profile.yaml> sbatch task_scripts/...
#
# Submit from the repository root, or export PDBCLEAN_REPO_ROOT.

set -euo pipefail

source "${PDBCLEAN_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}/task_scripts/stage14_common.sh"

pdbclean_banner "STAGE-14 REPRESENTATIVE SELECTION"

GRAPH="$(pdbclean_stage_argument representative_selection --graph-dir)"
OUT="$(pdbclean_stage_argument representative_selection --output-dir)"

if [[ ! -e "$GRAPH/_SUCCESS" ]]; then
    echo "ERROR: upstream Stage-14 graph has not passed: $GRAPH" >&2
    exit 1
fi

if [[ -e "$OUT" ]]; then
    echo "ERROR: Stage-14 representative output already exists: $OUT" >&2
    echo "Refusing to overwrite a completed stage." >&2
    exit 1
fi

pdbclean_run_stage representative_selection

echo
echo "Finished: $(date --iso-8601=seconds)"
echo "STAGE-14 REPRESENTATIVE JOB: PASS"
