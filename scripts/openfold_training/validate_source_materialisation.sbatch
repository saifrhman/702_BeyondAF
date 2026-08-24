#!/bin/bash -l
#SBATCH --job-name=of26_validate
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/users/%u/fastscratch/COMP702_pdbclean_pipeline/outputs/openfold_training/20260101/pdbclean-dedup-v1/source_materialisation/validation_%j.out
#SBATCH --error=/users/%u/fastscratch/COMP702_pdbclean_pipeline/outputs/openfold_training/20260101/pdbclean-dedup-v1/source_materialisation/validation_%j.err

set -euo pipefail

ROOT="$HOME/fastscratch/COMP702_pdbclean_pipeline"

TRAIN="$ROOT/outputs/openfold_training/20260101/pdbclean-dedup-v1"
MAT="$TRAIN/source_materialisation"

WORKLIST="$MAT/source_worklist.parquet"
STATUS="$MAT/status"

PY="$HOME/fastscratch/envs/bri_env_1.2.2/bin/python"

export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1
export XDG_CACHE_HOME="$HOME/fastscratch/.cache"

cd "$ROOT"

echo "========================================"
echo "OpenFold source materialisation validation"
echo "Job:     ${SLURM_JOB_ID}"
echo "Node:    $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
echo "========================================"

"$PY" - "$WORKLIST" "$STATUS" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

import pyarrow.parquet as pq


EXPECTED_BATCHES = 1182
EXPECTED_PDBS = 118197
EXPECTED_WORKLIST_SHA = (
    "01f20d56079fd1e8c02f99dbe8d1fb0abf59fbd00596f3d5e278ab7d7f579766"
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


worklist_path = Path(sys.argv[1]).resolve()
status_dir = Path(sys.argv[2]).resolve()

print("===== FINAL SOURCE MATERIALISATION AUDIT =====")
print("worklist:", worklist_path)
print("status dir:", status_dir)
print()

worklist_sha = sha256_file(worklist_path)

print("worklist sha256:", worklist_sha)

if worklist_sha != EXPECTED_WORKLIST_SHA:
    raise RuntimeError(
        "Frozen worklist SHA256 mismatch"
    )

worklist = pq.read_table(worklist_path).to_pylist()

if len(worklist) != EXPECTED_PDBS:
    raise RuntimeError(
        f"Expected {EXPECTED_PDBS} worklist rows, got {len(worklist)}"
    )

work_by_batch = {}
work_by_pdb = {}

for row in worklist:
    batch_id = int(row["batch_id"])
    pdb_id = row["pdb_id"]

    work_by_batch.setdefault(batch_id, []).append(row)

    if pdb_id in work_by_pdb:
        raise RuntimeError(
            f"Duplicate PDB in worklist: {pdb_id}"
        )

    work_by_pdb[pdb_id] = row

expected_batch_ids = set(range(EXPECTED_BATCHES))

if set(work_by_batch) != expected_batch_ids:
    missing = sorted(
        expected_batch_ids - set(work_by_batch)
    )
    extra = sorted(
        set(work_by_batch) - expected_batch_ids
    )

    raise RuntimeError(
        f"Worklist batch IDs invalid: "
        f"missing={missing[:20]} extra={extra[:20]}"
    )

reports = sorted(
    status_dir.glob("batch_*.json")
)

print("PASS report files:", len(reports))

if len(reports) != EXPECTED_BATCHES:
    raise RuntimeError(
        f"Expected {EXPECTED_BATCHES} reports, got {len(reports)}"
    )

seen_batches = set()
seen_pdbs = set()

missing_files = []
size_mismatches = []
report_errors = []
provenance_errors = []

for report_path in reports:
    report = json.loads(
        report_path.read_text()
    )

    if report.get("status") != "PASS":
        report_errors.append(
            (str(report_path), "status_not_PASS")
        )
        continue

    batch_id = int(report["batch_id"])

    if batch_id in seen_batches:
        report_errors.append(
            (batch_id, "duplicate_batch_report")
        )
        continue

    seen_batches.add(batch_id)

    if (
        report.get("worklist_sha256")
        != EXPECTED_WORKLIST_SHA
    ):
        report_errors.append(
            (batch_id, "wrong_worklist_sha256")
        )

    expected_rows = work_by_batch.get(batch_id)

    if expected_rows is None:
        report_errors.append(
            (batch_id, "batch_not_in_worklist")
        )
        continue

    expected_by_pdb = {
        r["pdb_id"]: r
        for r in expected_rows
    }

    outputs = report.get("outputs", [])

    if len(outputs) != len(expected_rows):
        report_errors.append(
            (
                batch_id,
                f"row_count:{len(outputs)}"
                f"!=expected:{len(expected_rows)}",
            )
        )

    for output in outputs:
        pdb_id = output["pdb_id"]

        if pdb_id in seen_pdbs:
            report_errors.append(
                (pdb_id, "duplicate_output")
            )
            continue

        seen_pdbs.add(pdb_id)

        expected = expected_by_pdb.get(pdb_id)

        if expected is None:
            report_errors.append(
                (pdb_id, "not_expected_in_batch")
            )
            continue

        if output["s3_key"] != expected["s3_key"]:
            provenance_errors.append(
                (pdb_id, "s3_key")
            )

        if output["source_etag"] != expected["etag"]:
            provenance_errors.append(
                (pdb_id, "etag")
            )

        if (
            int(output["compressed_size_bytes"])
            != int(expected["size_bytes"])
        ):
            provenance_errors.append(
                (pdb_id, "compressed_size")
            )

        expected_path = Path(
            expected["output_cif_path"]
        )

        output_path = Path(
            output["output_cif_path"]
        )

        if output_path != expected_path:
            provenance_errors.append(
                (pdb_id, "output_path")
            )

        if not output_path.is_file():
            missing_files.append(
                str(output_path)
            )
            continue

        actual_size = output_path.stat().st_size
        reported_size = int(
            output["cif_size_bytes"]
        )

        if actual_size != reported_size:
            size_mismatches.append(
                (
                    pdb_id,
                    actual_size,
                    reported_size,
                )
            )

print()
print("===== COUNTS =====")
print("expected batches:", EXPECTED_BATCHES)
print("validated batches:", len(seen_batches))
print("expected PDB entries:", EXPECTED_PDBS)
print("validated unique outputs:", len(seen_pdbs))
print("missing files:", len(missing_files))
print("size mismatches:", len(size_mismatches))
print("report errors:", len(report_errors))
print("provenance errors:", len(provenance_errors))

missing_pdbs = sorted(
    set(work_by_pdb) - seen_pdbs
)

extra_pdbs = sorted(
    seen_pdbs - set(work_by_pdb)
)

print("missing PDB IDs:", len(missing_pdbs))
print("extra PDB IDs:", len(extra_pdbs))

if missing_files:
    print()
    print("FIRST MISSING FILES:")
    for x in missing_files[:20]:
        print(x)

if size_mismatches:
    print()
    print("FIRST SIZE MISMATCHES:")
    for x in size_mismatches[:20]:
        print(x)

if report_errors:
    print()
    print("FIRST REPORT ERRORS:")
    for x in report_errors[:20]:
        print(x)

if provenance_errors:
    print()
    print("FIRST PROVENANCE ERRORS:")
    for x in provenance_errors[:20]:
        print(x)

if missing_pdbs:
    print()
    print(
        "FIRST MISSING PDB IDs:",
        missing_pdbs[:20],
    )

if extra_pdbs:
    print()
    print(
        "FIRST EXTRA PDB IDs:",
        extra_pdbs[:20],
    )

passed = all([
    len(seen_batches) == EXPECTED_BATCHES,
    len(seen_pdbs) == EXPECTED_PDBS,
    not missing_files,
    not size_mismatches,
    not report_errors,
    not provenance_errors,
    not missing_pdbs,
    not extra_pdbs,
])

print()
print("===== DECISION =====")

if not passed:
    print("STATUS: FAIL")
    raise SystemExit(1)

print("STATUS: PASS")
print(
    "All 118197 source mmCIFs are present and "
    "consistent with the frozen 20260101 worklist."
)
print()
print(
    "===== FINAL SOURCE MATERIALISATION AUDIT COMPLETE ====="
)
PY

echo
echo "Completed: $(date --iso-8601=seconds)"
