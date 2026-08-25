#!/bin/bash -l
#SBATCH --job-name=of_proj_index
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00

set -euo pipefail

REPO=/mnt/fastscratch/users/sgsrehm1/COMP702_pdbclean_pipeline
BASE="$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1"
OUT="$BASE/training_view"

export TMPDIR=/mnt/fastscratch/users/sgsrehm1/openfold_cache/tmp
mkdir -p "$TMPDIR" "$OUT"

PYTHON=/mnt/fastscratch/users/sgsrehm1/envs/bri_env_1.2.2/bin/python

echo "host        : $(hostname)"
echo "started     : $(date -Is)"
echo "pdbclean git: $(git -C "$REPO" rev-parse HEAD)"

"$PYTHON" "$REPO/scripts/openfold_training/build_projection_index.py" \
    --training-chains "$BASE/openfold_training_chains.parquet" \
    --output "$OUT/projection_index.json" \
    --manifest "$OUT/projection_index_manifest.json" \
    --expected-chains 499770

echo "finished    : $(date -Is)"
