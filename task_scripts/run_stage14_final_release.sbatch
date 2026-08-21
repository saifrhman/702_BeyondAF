#!/bin/bash
#SBATCH --job-name=pdbclean_release
#SBATCH --partition=nodes
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/stage14_release_%j.out
#SBATCH --error=logs/stage14_release_%j.err

# Gold release: the deduplicated retained-chain publication.
#
# Every path and every scientific argument is derived from the resolved run
# configuration, so this wrapper cannot drift from the CLI or the UI.
#
#   sbatch task_scripts/run_stage14_final_release.sbatch
#   PDBCLEAN_PROFILE=<profile.yaml> sbatch task_scripts/...
#
# Submit from the repository root, or export PDBCLEAN_REPO_ROOT.

set -euo pipefail

source "${PDBCLEAN_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}/task_scripts/stage14_common.sh"

pdbclean_banner "STAGE-14 FINAL RELEASE"

RELEASE="$(pdbclean_stage_argument gold_release --output-dir)"

if [[ -e "$RELEASE" ]]; then
    echo "ERROR: release already exists: $RELEASE" >&2
    echo "A published release is immutable; refusing to overwrite it." >&2
    exit 1
fi

pdbclean_run_stage gold_release

echo
echo "Finished: $(date --iso-8601=seconds)"
echo "FINAL RELEASE JOB: PASS"
