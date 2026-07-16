#!/bin/bash -l
#SBATCH -J of_relax_all
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH -w gpu32
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH --output=/users/%u/fastscratch/of_relax_all_%j.log
#SBATCH --error=/users/%u/fastscratch/of_relax_all_%j.err

set -euo pipefail

########################################
# Activate environment
########################################
set +u
source /mnt/data2/users/sgmwu14/tools/miniforge3/etc/profile.d/conda.sh
conda activate openfold_env
set -u

export PYTHONPATH=/users/sgmwu14/openfold:$PYTHONPATH

########################################
# Root directory
########################################
ROOT=/users/sgmwu14/scratch/openfold_inference_runs5times

echo "[$(date)] ===== START RELAXATION ====="
echo "[$(date)] ROOT = $ROOT"

########################################
# Check root exists
########################################
if [[ ! -d "$ROOT" ]]; then
    echo "[$(date)] ERROR: directory not found: $ROOT"
    exit 1
fi

########################################
# Loop over all unrelaxed pdb files
########################################
find "$ROOT" -type f -name "*_unrelaxed.pdb" | while read -r pdb; do

    relaxed_pdb="${pdb/_unrelaxed/_relaxed}"

    if [[ -f "$relaxed_pdb" ]]; then
        echo "[$(date)] Skip (already relaxed): $relaxed_pdb"
        continue
    fi

    echo "[$(date)] Relaxing: $pdb"

    python3 - <<EOF
from openfold.np import protein
from openfold.np.relax.relax import AmberRelaxation

try:
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

    relaxed_pdb_str, _, violations = relaxer.process(prot=prot)

    with open(r"$relaxed_pdb", "w") as f:
        f.write(relaxed_pdb_str)

    print("Saved:", r"$relaxed_pdb")
    print("Violations:", float(violations.sum()))

except Exception as e:
    print("ERROR processing:", r"$pdb")
    print(str(e))
EOF

done

echo "[$(date)] ===== ALL RELAXATION DONE ====="