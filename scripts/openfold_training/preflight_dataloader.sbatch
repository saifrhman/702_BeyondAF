#!/bin/bash -l
#SBATCH --job-name=of_preflight
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00

set -euo pipefail

REPO=/mnt/fastscratch/users/sgsrehm1/COMP702_pdbclean_pipeline
BASE="$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1"
VIEW="$BASE/training_view"
ENV=/mnt/fastscratch/users/sgsrehm1/envs/openfold_train

. "$ENV/etc/conda/activate.d/10-cache-dirs.sh"
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

echo "host    : $(hostname)"
echo "started : $(date -Is)"

"$ENV/bin/python" "$REPO/scripts/openfold_training/preflight_dataloader.py" \
    --mmcif-dir /users/sgsrehm1/scratch/COMP702_openfold_20260101/pdbclean-dedup-v1/mmcif \
    --msa-store /mnt/fastscratch/users/sgsrehm1/COMP702_openfold_msa/20260101/pdbclean-dedup-v1/msas \
    --alignment-index "$BASE/alignment_index/pdbclean_dedup_v1.index" \
    --projection-index "$VIEW/projection_index.json" \
    --train-filter "$VIEW/train_filter_view.txt" \
    --samples 12 \
    --report "$VIEW/dataloader_preflight.json"

echo "finished: $(date -Is)"
