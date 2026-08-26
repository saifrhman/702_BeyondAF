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

# Host C++ compiler for the JIT build. Prepending $ENV/bin above is required for
# ninja, but it also shadows /usr/bin/gcc with conda's -- and this env has
# gcc_impl_linux-64 without gxx_impl_linux-64, so $ENV/libexec/gcc ships cc1,
# cc1obj and cc1objplus but no cc1plus. Any C++ source then dies with
# "gcc: fatal error: cannot execute 'cc1plus'", which reads like a missing
# compiler but is really a half-installed one winning the PATH race.
#
# The node's own gcc is complete, so point the build at it explicitly rather
# than installing a compiler into a conda env that a training job is currently
# running out of. The extension links against conda's libstdc++ (15.2) at run
# time, which is newer than gcc 11's, and that direction is ABI-compatible.
export CC=/usr/bin/gcc
export CXX=/usr/bin/g++
export NVCC_PREPEND_FLAGS="-ccbin ${CXX}"

if ! "$CXX" -x c++ -c /dev/null -o /dev/null 2>/dev/null; then
    echo "FATAL: $CXX cannot compile C++ on $(hostname); no cc1plus." >&2
    echo "       Not writing a verdict -- training keeps the unfused path." >&2
    exit 3
fi

echo "host c++    : $("$CXX" --version | head -1)"

# The DeepSpeed evoformer op is JIT-built from CUTLASS, which ships in the repo.
export CUTLASS_PATH="$OFSRC/cutlass"

echo "host        : $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo "CUTLASS_PATH: $CUTLASS_PATH"
echo "extensions  : $TORCH_EXTENSIONS_DIR"

cd "$OFSRC"

# Benchmark against the run's real weights, not random ones, so the loss
# agreement is measured where the model actually is. The live checkpoint is
# rewritten every epoch by the training job, though, and reading a 1.5 GB file
# mid-write yields a torn copy. Snapshot it, prove the snapshot loads, and fall
# back to random initialisation rather than trusting a half-written file -- the
# kernel comparison stays valid either way, because both paths get identical
# weights from an identical seed.
LIVE=$(find /mnt/fastscratch/users/sgsrehm1/openfold_runs/pdbclean_dedup_v1_scratch \
       -name '*.ckpt' -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -1 | cut -d' ' -f2-)

CKPT=""
if [[ -n "$LIVE" ]]; then
    SNAP="$TMPDIR/attn_bench_$SLURM_JOB_ID.ckpt"
    cp "$LIVE" "$SNAP" 2>/dev/null || true

    if "$ENV/bin/python" -c "
import sys, torch
torch.load(sys.argv[1], map_location='cpu', weights_only=False)['state_dict']
" "$SNAP" 2>/dev/null; then
        CKPT="$SNAP"
        echo "checkpoint  : $LIVE (snapshot verified)"
    else
        rm -f "$SNAP"
        echo "checkpoint  : $LIVE was mid-write; using random init instead"
    fi
fi

[[ -n "$CKPT" ]] || echo "checkpoint  : <none -- random initialisation>"
trap 'rm -f "$TMPDIR/attn_bench_$SLURM_JOB_ID.ckpt"' EXIT

"$ENV/bin/python" "$REPO/scripts/openfold_training/benchmark_attention_kernel.py" \
    --mmcif-dir /users/sgsrehm1/scratch/COMP702_openfold_20260101/pdbclean-dedup-v1/mmcif \
    --msa-store /mnt/fastscratch/users/sgsrehm1/COMP702_openfold_msa/20260101/pdbclean-dedup-v1/msas \
    --alignment-index "$BASE/alignment_index/pdbclean_dedup_v1.index" \
    --projection-index "$VIEW/projection_index.json" \
    --train-filter "$VIEW/train_filter_view.txt" \
    ${CKPT:+--checkpoint "$CKPT"} \
    --steps 4 \
    --report "$VIEW/attention_kernel_benchmark.json"
