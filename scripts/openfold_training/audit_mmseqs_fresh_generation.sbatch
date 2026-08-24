#!/bin/bash -l
#SBATCH --job-name=of_mmseqs_audit
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail

TRAIN_ROOT="${1:?training root required}"
OPENFOLD="${2:?OpenFold root required}"
OFENV="${3:?OpenFold environment required}"

TRAIN_ROOT="$(readlink -f "$TRAIN_ROOT")"
OPENFOLD="$(readlink -f "$OPENFOLD")"
OFENV="$(readlink -f "$OFENV")"

PYTHON="$OFENV/bin/python"
MMSEQS="$OFENV/bin/mmseqs"

if [[ ! -d "$TRAIN_ROOT" ]]; then
    echo "ERROR: training root missing: $TRAIN_ROOT"
    exit 1
fi

if [[ ! -d "$OPENFOLD/openfold" ]]; then
    echo "ERROR: OpenFold Python package missing: $OPENFOLD/openfold"
    exit 1
fi

if [[ ! -f "$OPENFOLD/openfold/data/tools/hhsearch.py" ]]; then
    echo "ERROR: hhsearch module missing"
    exit 1
fi

if [[ ! -f "$OPENFOLD/scripts/precompute_alignments_mmseqs.py" ]]; then
    echo "ERROR: MMseqs OpenFold driver missing"
    exit 1
fi

if [[ ! -f "$OPENFOLD/scripts/colabfold_search.sh" ]]; then
    echo "ERROR: colabfold_search.sh missing"
    exit 1
fi

if [[ ! -x "$PYTHON" ]]; then
    echo "ERROR: Python missing: $PYTHON"
    exit 1
fi

if [[ ! -x "$MMSEQS" ]]; then
    echo "ERROR: MMseqs2 missing: $MMSEQS"
    exit 1
fi

export PYTHONPATH="$OPENFOLD${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1

echo "============================================================"
echo "FRESH MMSEQS2 MSA GENERATION AUDIT"
echo "Job:          ${SLURM_JOB_ID}"
echo "Node:         $(hostname)"
echo "Training:     $TRAIN_ROOT"
echo "OpenFold:     $OPENFOLD"
echo "Environment:  $OFENV"
echo "Started:      $(date --iso-8601=seconds)"
echo "============================================================"

echo
echo "===== OPENFOLD PYTHON IMPORT ====="

"$PYTHON" - <<'PY'
import openfold
from openfold.data.tools import hhsearch

print("openfold package:", openfold.__file__)
print("hhsearch module:", hhsearch.__file__)
print("STATUS: IMPORT PASS")
PY

echo
echo "===== MMSEQS2 ====="
echo "binary: $MMSEQS"
echo -n "version: "
"$MMSEQS" version

echo
echo "===== CURRENT SNAPSHOT QUERY POPULATION ====="

QUERY_MANIFEST="$TRAIN_ROOT/msa/msa_query_manifest.json"
QUERY_FASTA="$TRAIN_ROOT/msa/msa_queries.fasta"

if [[ ! -f "$QUERY_MANIFEST" ]]; then
    echo "ERROR: query manifest missing: $QUERY_MANIFEST"
    exit 1
fi

if [[ ! -f "$QUERY_FASTA" ]]; then
    echo "ERROR: query FASTA missing: $QUERY_FASTA"
    exit 1
fi

cat "$QUERY_MANIFEST"

echo
echo "FASTA:"
ls -lh "$QUERY_FASTA"

echo
echo "===== OPENFOLD MMSEQS DRIVER HELP ====="

HELP_FILE="$TRAIN_ROOT/msa/audit/mmseqs_driver_help_${SLURM_JOB_ID}.txt"

cd "$OPENFOLD"

"$PYTHON" \
    scripts/precompute_alignments_mmseqs.py \
    --help \
    > "$HELP_FILE"

cat "$HELP_FILE"
rm -f "$HELP_FILE"

echo
echo "===== COLABFOLD SEARCH PROTOCOL ====="
sed -n '1,240p' "$OPENFOLD/scripts/colabfold_search.sh"

echo
echo "===== MMSEQS DATABASE DOWNLOAD DRIVER ====="
sed -n '1,220p' "$OPENFOLD/scripts/download_mmseqs_dbs.sh"

echo
echo "===== UNIREF30 DATABASE DEFINITION ====="
sed -n '1,180p' "$OPENFOLD/scripts/download_uniref30.sh"

echo
echo "===== ENVIRONMENTAL DATABASE DEFINITION ====="
sed -n '1,180p' "$OPENFOLD/scripts/download_colabfold_envdb.sh"

echo
echo "===== DOWNLOAD TOOL AVAILABILITY ====="

if command -v aria2c >/dev/null 2>&1; then
    echo "aria2c: $(command -v aria2c)"
else
    echo "aria2c: NOT ON PATH"
fi

if command -v curl >/dev/null 2>&1; then
    echo "curl: $(command -v curl)"
else
    echo "curl: NOT ON PATH"
fi

if command -v wget >/dev/null 2>&1; then
    echo "wget: $(command -v wget)"
else
    echo "wget: NOT ON PATH"
fi

echo
echo "===== EXISTING MMSEQS DATABASE MATERIAL ====="

for ROOT in "$HOME/scratch" "$HOME/fastscratch"; do
    echo
    echo "--- $ROOT ---"

    if [[ ! -d "$ROOT" ]]; then
        echo "directory absent"
        continue
    fi

    find "$ROOT" \
        -maxdepth 6 \
        \( \
            -iname 'UniRef30_2021_03*' \
            -o -iname '*uniref30*' \
            -o -iname '*colabfold_envdb*' \
            -o -iname '*bfd.mgnify30.metaeuk30.smag30*' \
        \) \
        -print \
        2>/dev/null \
        || true
done

echo
echo "===== STORAGE ====="

for ROOT in "$HOME/scratch" "$HOME/fastscratch"; do
    if [[ -e "$ROOT" ]]; then
        df -h "$ROOT"
    fi
done

echo
echo "============================================================"
echo "STATUS: PASS"
echo "Completed: $(date --iso-8601=seconds)"
echo "============================================================"
