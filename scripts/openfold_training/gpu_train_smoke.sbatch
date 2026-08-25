#!/bin/bash -l
#SBATCH --job-name=of_train_smoke
#SBATCH --partition=gpu-l40s-low
#SBATCH --gres=gpu:l40s:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00

# A real training smoke test on the real PDBClean population.
#
# Nothing here is synthetic: the chains, structures, MSAs, alignment index,
# retained-chain projection and chain cache are the production artefacts. Only
# the *budget* is small -- a handful of optimizer steps -- so that the dataloader,
# the forward pass, the loss, the backward pass, the optimizer step and
# checkpointing are each exercised once against genuine data.
#
# The model is initialised randomly. No pretrained AlphaFold or OpenFold weights
# are loaded: --resume_from_ckpt is deliberately not passed, because the
# scientific goal is a new model trained on the geometry-cleaned dataset.

set -euo pipefail

REPO=/mnt/fastscratch/users/sgsrehm1/COMP702_pdbclean_pipeline
BASE="$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1"
VIEW="$BASE/training_view"
OFSRC=/mnt/fastscratch/users/sgsrehm1/openfold_src

MMCIF=/users/sgsrehm1/scratch/COMP702_openfold_20260101/pdbclean-dedup-v1/mmcif
MSAS=/mnt/fastscratch/users/sgsrehm1/COMP702_openfold_msa/20260101/pdbclean-dedup-v1/msas

RUN_ROOT=/mnt/fastscratch/users/sgsrehm1/openfold_runs/smoke_${SLURM_JOB_ID}

ENV=/mnt/fastscratch/users/sgsrehm1/envs/openfold_train
. "$ENV/etc/conda/activate.d/10-cache-dirs.sh"

export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

mkdir -p "$RUN_ROOT"

echo "=== identity ==="
echo "host              : $(hostname)"
echo "slurm job         : ${SLURM_JOB_ID}"
echo "started           : $(date -Is)"
echo "pdbclean commit   : $(git -C "$REPO" rev-parse HEAD)"
echo "openfold source   : $OFSRC (from commit da89cd28446abcd7be95459b7b349dedaee666c0)"
echo "run root          : $RUN_ROOT"

echo
echo "=== gpu ==="
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader

echo
echo "=== dataset identity (hashes) ==="
for f in \
    "$BASE/alignment_index/pdbclean_dedup_v1.index" \
    "$VIEW/projection_index.json" \
    "$VIEW/train_chain_data_cache.json" \
    "$VIEW/train_filter_view.txt" \
    "$BASE/openfold_training_chains.parquet"
do
    if [[ -f "$f" ]]; then
        echo "$(sha256sum "$f")"
    else
        echo "MISSING: $f"
        exit 1
    fi
done

echo
echo "=== torch / cuda ==="
"$ENV/bin/python" - <<'PY'
import torch
print("torch          :", torch.__version__)
print("cuda build     :", torch.version.cuda)
print("cuda available :", torch.cuda.is_available())
print("device         :", torch.cuda.get_device_name(0))
print("capability     :", torch.cuda.get_device_capability(0))
PY

echo
echo "=== training ==="

cd "$OFSRC"

set -x
"$ENV/bin/python" train_openfold.py \
    "$MMCIF" \
    "$MSAS" \
    "$MMCIF" \
    "$RUN_ROOT" \
    2026-01-01 \
    --config_preset initial_training \
    --alignment_index_path "$BASE/alignment_index/pdbclean_dedup_v1.index" \
    --train_projection_index_path "$VIEW/projection_index.json" \
    --train_chain_data_cache_path "$VIEW/train_chain_data_cache.json" \
    --train_filter_path "$VIEW/train_filter_view.txt" \
    --experiment_config_json "$REPO/scripts/openfold_training/gpu_smoke_config.json" \
    --train_epoch_len 8 \
    --max_epochs 1 \
    --gpus 1 \
    --num_nodes 1 \
    --precision bf16 \
    --seed 20260101 \
    --log_every_n_steps 1 \
    --num_sanity_val_steps 0 \
    --checkpoint_every_epoch
set +x

echo
echo "=== checkpoints written ==="
find "$RUN_ROOT" -name '*.ckpt' -o -name '*.pt' | head -20

echo "finished          : $(date -Is)"
