#!/bin/bash -l
#SBATCH --job-name=of_attn_bench
#SBATCH --partition=gpu-l40s-low
#SBATCH --gres=gpu:l40s:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00

set -euo pipefail

REPO=/mnt/fastscratch/users/sgsrehm1/COMP702_pdbclean_pipeline
BASE="$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1"
VIEW="$BASE/training_view"
OFSRC=/mnt/fastscratch/users/sgsrehm1/openfold_src
ENV=/mnt/fastscratch/users/sgsrehm1/envs/openfold_train

. "$ENV/etc/conda/activate.d/10-cache-dirs.sh"
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"
# torch JIT-compiles CUDA extensions by shelling out to ninja, so the env bin
# must be on PATH -- sourcing the cache-dirs hook alone does not do that, and
# the failure surfaces only as "Ninja is required to load C++ extensions".
export PATH="$ENV/bin:$PATH"


# The DeepSpeed evoformer op is JIT-built from CUTLASS, which ships in the repo.
export CUTLASS_PATH="$OFSRC/cutlass"

echo "host        : $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo "CUTLASS_PATH: $CUTLASS_PATH"
echo "extensions  : $TORCH_EXTENSIONS_DIR"

cd "$OFSRC"

CKPT=$(find /mnt/fastscratch/users/sgsrehm1/openfold_runs/pdbclean_dedup_v1_scratch \
       -name '*.ckpt' -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -1 | cut -d' ' -f2-)
echo "checkpoint  : ${CKPT:-<none>}"

"$ENV/bin/python" "$REPO/scripts/openfold_training/benchmark_attention_kernel.py" \
    --mmcif-dir /users/sgsrehm1/scratch/COMP702_openfold_20260101/pdbclean-dedup-v1/mmcif \
    --msa-store /mnt/fastscratch/users/sgsrehm1/COMP702_openfold_msa/20260101/pdbclean-dedup-v1/msas \
    --alignment-index "$BASE/alignment_index/pdbclean_dedup_v1.index" \
    --projection-index "$VIEW/projection_index.json" \
    --train-filter "$VIEW/train_filter_view.txt" \
    ${CKPT:+--checkpoint "$CKPT"} \
    --steps 4 \
    --report "$VIEW/attention_kernel_benchmark.json"
