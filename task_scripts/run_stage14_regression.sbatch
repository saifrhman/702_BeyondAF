#!/bin/bash -l
#SBATCH --job-name=stage14_regression
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/stage14_regression_%j.out
#SBATCH --error=logs/stage14_regression_%j.err

# Scientific regression evidence for the configuration refactor.
#
# Re-runs the three Stage-14 entry points with the refactored, config-driven
# code and compares every artefact against the frozen release.
#
# SAFETY:
#   * All output goes to a scratch regression root.
#   * The frozen release and the frozen protocol outputs are read only.
#   * Nothing under outputs/releases or outputs/pdbclean is written.

set -euo pipefail

# Slurm copies the batch script into its spool directory, so the repository
# must never be inferred from BASH_SOURCE. Resolve it from the environment or
# from the submission directory instead.
if [[ -z "${PDBCLEAN_REPO_ROOT:-}" ]]; then
    PDBCLEAN_REPO_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
fi

PDBCLEAN_REPO_ROOT="$(cd "$PDBCLEAN_REPO_ROOT" && pwd)"
export PDBCLEAN_REPO_ROOT

if [[ ! -f "$PDBCLEAN_REPO_ROOT/config/pdbclean/pipeline_env.sh" ]]; then
    echo "ERROR: not a PDBClean repository root: $PDBCLEAN_REPO_ROOT" >&2
    echo "Submit from the repository root, or export PDBCLEAN_REPO_ROOT." >&2
    exit 2
fi

source "$PDBCLEAN_REPO_ROOT/config/pdbclean/pipeline_env.sh"

cd "$PDBCLEAN_REPO_ROOT"

PY="${PDBCLEAN_PYTHON}"

if [[ ! -x "$PY" ]]; then
    PY="$(command -v "$PDBCLEAN_PYTHON")"
fi

REGRESSION_ROOT="${PDBCLEAN_REGRESSION_ROOT:-$HOME/fastscratch/COMP702_stage14_regression/${SLURM_JOB_ID:-manual}}"

PROTOCOL_ROOT="$PDBCLEAN_OUTPUT_ROOT/20260101/protocol3.2-comp702-v1"
FROZEN_RELEASE="$PDBCLEAN_RELEASE_ROOT/PDBClean-20260101-protocol3.2-comp702-v1-dedup-v1"

PROFILE="config/pdbclean/profiles/comp702_frozen_20260101.yaml"

mkdir -p "$REGRESSION_ROOT"

echo "============================================================"
echo "STAGE-14 SCIENTIFIC REGRESSION"
echo "Job:             ${SLURM_JOB_ID:-manual}"
echo "Node:            $(hostname)"
echo "Repository:      $PDBCLEAN_REPO_ROOT"
echo "Git HEAD:        $(git rev-parse HEAD)"
echo "Git dirty:       $(git status --porcelain | wc -l) path(s)"
echo "Python:          $PY"
echo "Protocol root:   $PROTOCOL_ROOT"
echo "Frozen release:  $FROZEN_RELEASE"
echo "Regression root: $REGRESSION_ROOT"
echo "Started:         $(date --iso-8601=seconds)"
echo "============================================================"

for required in \
    "$PROTOCOL_ROOT/full_bri_nn/finalized/candidate_near_duplicates.parquet" \
    "$PROTOCOL_ROOT/stage14_geometric_graph/global_summary.json" \
    "$FROZEN_RELEASE/release_manifest.json"
do
    if [[ ! -f "$required" ]]; then
        echo "ERROR: missing required frozen input: $required"
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# 1. Near-duplicate graph
# ---------------------------------------------------------------------------
GRAPH_OUT="$REGRESSION_ROOT/stage14_geometric_graph"

echo
echo "===== STAGE 14a: NEAR-DUPLICATE GRAPH ====="

rm -rf "$GRAPH_OUT"

"$PY" scripts/build_stage14_geometric_graph.py \
    --edges "$PROTOCOL_ROOT/full_bri_nn/finalized/candidate_near_duplicates.parquet" \
    --m1-edges "$PROTOCOL_ROOT/full_bri_nn/finalized/m1_near_duplicates.parquet" \
    --output-dir "$GRAPH_OUT" \
    --expected-edges 1068256 \
    --expected-m1-edges 4495 \
    --expected-mge2-nodes 577760

# ---------------------------------------------------------------------------
# 2. Representative selection
# ---------------------------------------------------------------------------
REPS_OUT="$REGRESSION_ROOT/stage14_representative_selection_v1"

echo
echo "===== STAGE 14b: REPRESENTATIVE SELECTION ====="

rm -rf "$REPS_OUT"

"$PY" scripts/select_stage14_representatives.py \
    --graph-dir "$GRAPH_OUT" \
    --edges "$PROTOCOL_ROOT/full_bri_nn/finalized/candidate_near_duplicates.parquet" \
    --accepted "$PROTOCOL_ROOT/quality/merged/accepted.parquet" \
    --metadata "$PROTOCOL_ROOT/downstream_metadata/finalized/entry_metadata.parquet" \
    --config "$PDBCLEAN_REPRESENTATIVE_POLICY" \
    --output-dir "$REPS_OUT" \
    --m1-retained-chain-count 764 \
    --expected-canonical-input-chains 578524

# ---------------------------------------------------------------------------
# 3. Comparison against the frozen artefacts
# ---------------------------------------------------------------------------
echo
echo "===== COMPARISON AGAINST FROZEN ARTEFACTS ====="

"$PY" scripts/compare_stage14_regression.py \
    --frozen-graph "$PROTOCOL_ROOT/stage14_geometric_graph" \
    --regression-graph "$GRAPH_OUT" \
    --frozen-representatives "$PROTOCOL_ROOT/stage14_representative_selection_v1" \
    --regression-representatives "$REPS_OUT" \
    --frozen-release "$FROZEN_RELEASE" \
    --report "$REGRESSION_ROOT/regression_report.json"

echo
echo "===== CLI / CONFIGURATION EQUIVALENCE ====="

"$PY" -m pdbclean.cli plan --config "$PROFILE" --json \
    > "$REGRESSION_ROOT/frozen_plan.json"

"$PY" - "$REGRESSION_ROOT/frozen_plan.json" <<'PY'
import json
import sys

plan = json.load(open(sys.argv[1]))

bad = [
    stage
    for stage in plan["stages"]
    if stage["action"] not in {"reuse", "not_applicable"}
]

print("snapshot:", plan["snapshot"])
print("resolved config sha256:", plan["resolved_config_sha256"])
print("scientific config sha256:", plan["scientific_config_sha256"])
print("stages:", len(plan["stages"]))
print("all reusable:", not bad)

# The scientific identity must not depend on which host resolved it. The
# full hash legitimately differs between the login node and a compute node
# because $TMPDIR expands per host; the scientific hash must not.
EXPECTED_SCIENTIFIC_SHA256 = (
    "25b8e62a87cb90797af41cd4149dfd4280e3a7aed99428e70fa97117c5bababa"
)

if plan["scientific_config_sha256"] != EXPECTED_SCIENTIFIC_SHA256:
    print("  SCIENTIFIC HASH MISMATCH")
    print("  expected:", EXPECTED_SCIENTIFIC_SHA256)
    print("  observed:", plan["scientific_config_sha256"])
    raise SystemExit(1)

if bad:
    for stage in bad:
        print("  NOT REUSABLE:", stage["stage_id"], stage["status"], stage["validation"])
    raise SystemExit(1)

print("STATUS: PASS")
PY

echo
echo "============================================================"
echo "STAGE-14 REGRESSION COMPLETE"
echo "Report:    $REGRESSION_ROOT/regression_report.json"
echo "Completed: $(date --iso-8601=seconds)"
echo "============================================================"
