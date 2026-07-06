#!/bin/bash -l
#SBATCH -J of_roda_dl
#SBATCH -p long
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/of_roda_dl_%j.log

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

echo "Job $SLURM_JOB_ID started on $(hostname) at $(date)"
echo "USER=$USER"

# AWS should be available in the job environment.
# If it is not available, load the relevant module or activate a conda environment before running this script.
which aws
aws --version

SRC="$RODA_S3_SRC"
DST="$RODA_ALIGNMENT_DST"

mkdir -p "$DST"

echo "Syncing from: $SRC"
echo "Syncing to:   $DST"

# Use --only-show-errors to avoid very large Slurm logs.
aws s3 sync \
  "$SRC" \
  "$DST" \
  --no-sign-request \
  --only-show-errors

echo "Job $SLURM_JOB_ID finished at $(date)"
