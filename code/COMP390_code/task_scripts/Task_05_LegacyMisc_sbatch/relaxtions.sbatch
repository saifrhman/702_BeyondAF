#!/bin/bash -l
#SBATCH -J of_relax_all
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/of_relax_all_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/of_relax_all_%j.err

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

echo "======================================"
echo "Job ID: ${SLURM_JOB_ID:-manual}"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
echo "======================================"

module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
conda activate "$OPENFOLD_ENV"

export PYTHONPATH="$OPENFOLD_CODE_DIR:${PYTHONPATH:-}"

ROOT="${OPENFOLD_RELAX_ROOT:-$OPENFOLD_RUNS5_ROOT}"

echo "[$(date)] Starting relaxation"
echo "[$(date)] ROOT=$ROOT"

if [[ ! -d "$ROOT" ]]; then
    echo "[$(date)] ERROR: directory not found: $ROOT"
    exit 1
fi

find "$ROOT" -type f -name "*_unrelaxed.pdb" | while read -r pdb; do
    relaxed_pdb="${pdb/_unrelaxed/_relaxed}"

    if [[ -f "$relaxed_pdb" ]]; then
        echo "[$(date)] Skip already relaxed file: $relaxed_pdb"
        continue
    fi

    echo "[$(date)] Relaxing: $pdb"

    PDB_IN="$pdb" PDB_OUT="$relaxed_pdb" python3 - <<'PYRELAX'
import os
from openfold.np import protein
from openfold.np.relax.relax import AmberRelaxation

pdb_in = os.environ["PDB_IN"]
pdb_out = os.environ["PDB_OUT"]

try:
    with open(pdb_in, "r") as f:
        pdb_str = f.read()

    prot = protein.from_pdb_string(pdb_str)

    relaxer = AmberRelaxation(
        max_iterations=0,
        tolerance=2.39,
        stiffness=10.0,
        exclude_residues=[],
        max_outer_iterations=20,
        use_gpu=True,
    )

    relaxed_pdb_str, _, violations = relaxer.process(prot=prot)

    with open(pdb_out, "w") as f:
        f.write(relaxed_pdb_str)

    print("Saved:", pdb_out)
    print("Violations:", float(violations.sum()))

except Exception as e:
    print("ERROR processing:", pdb_in)
    print(str(e))
PYRELAX

done

echo "[$(date)] All relaxation finished"
