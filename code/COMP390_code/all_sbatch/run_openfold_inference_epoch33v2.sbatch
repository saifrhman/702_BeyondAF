#!/bin/bash -l
#SBATCH -J of_infer3_ep33
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH -w gpu32
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/fastscratch/of_infer3_ep33_%j.log
#SBATCH --error=/users/%u/fastscratch/of_infer3_ep33_%j.err

set -euo pipefail

########################################
# Activate environment
########################################

set +u
source /mnt/data2/users/sgmwu14/tools/miniforge3/etc/profile.d/conda.sh
conda activate openfold_env
set -u

export PYTHONPATH=/users/sgmwu14/openfold:$PYTHONPATH
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

export TRITON_CACHE_DIR=/users/sgmwu14/fastscratch/triton_cache
mkdir -p "$TRITON_CACHE_DIR"

########################################
# Paths
########################################

FASTA_DIR=/users/sgmwu14/scratch/openfold_inference/fasta_dir
CIF_DIR=/users/sgmwu14/scratch/taskcif
ALIGN_DIR=/users/sgmwu14/scratch/openfold_inference/alignments_cf

# 新的总输出目录
RUN_ROOT=/users/sgmwu14/scratch/openfold_inference_ep33_rerun4

AF_OUT=${RUN_ROOT}/alphafold_output
OF_OUT=${RUN_ROOT}/openfold_output
MY_OUT=${RUN_ROOT}/your_model_output_ep33

MY_CKPT=/users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/27-280000.ckpt

mkdir -p "$RUN_ROOT" "$AF_OUT" "$OF_OUT" "$MY_OUT"

########################################
# 1) AlphaFold parameters (unrelaxed only)
########################################

echo "[$(date)] ===== Running AlphaFold params (unrelaxed only) ====="

python3 /users/sgmwu14/openfold/run_pretrained_openfold.py \
  "$FASTA_DIR" \
  "$CIF_DIR" \
  --output_dir "$AF_OUT" \
  --use_precomputed_alignments "$ALIGN_DIR" \
  --config_preset model_1 \
  --jax_param_path /users/sgmwu14/openfold/openfold/resources/params/params_model_1.npz \
  --model_device cuda:0 \
  --skip_relaxation

########################################
# 2) OpenFold official model (unrelaxed only)
########################################

echo "[$(date)] ===== Running OpenFold official (unrelaxed only) ====="

python3 /users/sgmwu14/openfold/run_pretrained_openfold.py \
  "$FASTA_DIR" \
  "$CIF_DIR" \
  --output_dir "$OF_OUT" \
  --use_precomputed_alignments "$ALIGN_DIR" \
  --config_preset model_1 \
  --openfold_checkpoint_path /users/sgmwu14/openfold/openfold/resources/openfold_params/finetuning_2.pt \
  --model_device cuda:0 \
  --skip_relaxation

########################################
# 3) Your trained model (ep33, unrelaxed only)
########################################

echo "[$(date)] ===== Running Your Model ep33 (unrelaxed only) ====="

python3 /users/sgmwu14/openfold/run_pretrained_openfold.py \
  "$FASTA_DIR" \
  "$CIF_DIR" \
  --output_dir "$MY_OUT" \
  --use_precomputed_alignments "$ALIGN_DIR" \
  --config_preset model_1 \
  --openfold_checkpoint_path "$MY_CKPT" \
  --model_device cuda:0 \
  --skip_relaxation

echo "[$(date)] ===== ALL UNRELAXED INFERENCE FINISHED ====="

########################################
# 4) Put each unrelaxed file into its own folder,
#    then run 5x relaxation in the same folder
########################################

echo "[$(date)] ===== Organising outputs and running 5x relaxation ====="

for MODEL_DIR in "$AF_OUT" "$OF_OUT" "$MY_OUT"; do
    MODEL_NAME=$(basename "$MODEL_DIR")
    PRED_DIR="${MODEL_DIR}/predictions"

    echo "[$(date)] ----- Processing model: ${MODEL_NAME} -----"

    find "$PRED_DIR" -maxdepth 1 -type f -name "*unrelaxed*.pdb" | while read -r pdb; do
        fname=$(basename "$pdb")
        base=$(basename "$pdb" .pdb)

        # 为每个蛋白/结果建立单独文件夹
        TARGET_DIR="${PRED_DIR}/${base}"
        mkdir -p "$TARGET_DIR"

        # 把 unrelaxed 文件移进去（如果还没移）
        if [ ! -f "${TARGET_DIR}/${fname}" ]; then
            mv "$pdb" "${TARGET_DIR}/"
        fi

        # 更新路径
        pdb="${TARGET_DIR}/${fname}"

        echo "[$(date)] Processing ${MODEL_NAME}: ${base}"
        echo "[$(date)] Unrelaxed file: ${pdb}"

        for i in 1 2 3 4 5; do
            out_pdb="${TARGET_DIR}/${base}_relax${i}.pdb"
            out_dbg="${TARGET_DIR}/${base}_relax${i}_debug.txt"

            echo "[$(date)]   Relax run ${i} -> ${out_pdb}"

            python3 - <<EOF || echo "[$(date)] WARNING: relax failed for ${MODEL_NAME} ${base} run ${i}"
from openfold.np import protein
from openfold.np.relax.relax import AmberRelaxation

with open(r"$pdb", "r") as f:
    pdb_str = f.read()

prot = protein.from_pdb_string(pdb_str)

relaxer = AmberRelaxation(
    max_iterations=0,
    tolerance=2.39,
    stiffness=10.0,
    exclude_residues=[],
    max_outer_iterations=20,
    use_gpu=True
)

relaxed_pdb, debug_data, violations = relaxer.process(prot=prot)

with open(r"$out_pdb", "w") as f:
    f.write(relaxed_pdb)

with open(r"$out_dbg", "w") as f:
    f.write(str(debug_data) + "\\n")
    f.write("violations_sum=" + str(violations.sum()) + "\\n")
EOF
        done
    done
done

echo "[$(date)] ===== ALL 5x RELAXATION FINISHED ====="