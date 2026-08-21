#!/bin/bash
#SBATCH --job-name=stage14_graph
#SBATCH --partition=nodes
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/stage14_graph_%j.out
#SBATCH --error=logs/stage14_graph_%j.err

# Stage 13: near-duplicate redundancy graph.
#
# Every path and every scientific argument is derived from the resolved run
# configuration, so this wrapper cannot drift from the CLI or the UI.
#
#   sbatch task_scripts/run_stage14_geometric_graph.sbatch
#   PDBCLEAN_PROFILE=<profile.yaml> sbatch task_scripts/...
#
# Submit from the repository root, or export PDBCLEAN_REPO_ROOT.

set -euo pipefail

source "${PDBCLEAN_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}/task_scripts/stage14_common.sh"

pdbclean_banner "STAGE-14 NEAR-DUPLICATE GRAPH"

OUT="$(pdbclean_stage_argument redundancy_graph --output-dir)"

if [[ -e "$OUT/_SUCCESS" ]]; then
    echo "ERROR: finalized Stage-14 graph already exists: $OUT" >&2
    echo "Refusing to overwrite a completed stage." >&2
    exit 1
fi

pdbclean_run_stage redundancy_graph

echo
echo "Finished: $(date --iso-8601=seconds)"
echo "STAGE-14 GRAPH JOB: PASS"
