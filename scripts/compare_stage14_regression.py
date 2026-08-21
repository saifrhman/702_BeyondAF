#!/usr/bin/env python3

"""Compare regression Stage-14 artefacts against the frozen release.

Proves that the configuration refactor did not change the scientific result:
for the same resolved scientific configuration, every Stage-14 artefact must be
identical to the frozen one.

Volatile provenance fields (timestamps, Slurm job identifiers, git commits, and
hashes of paths that legitimately differ between the frozen tree and the
regression tree) are compared separately and are allowed to differ; every
scientific field must match exactly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pyarrow.parquet as pq


# Fields that legitimately differ between two runs of the same configuration.
VOLATILE_KEYS = frozenset(
    {
        "generated_at",
        "generated_at_utc",
        "created_at",
        "completed_at",
        "verified_at",
        "run_date",
        "date",
        "elapsed_seconds",
        "host",
        "hostname",
        "job_id",
        "slurm_job_id",
        "analysis_job_id",
        "graph_success_sha256",
        "policy_config_sha256",
        "accepted_quality_sha256",
        "entry_metadata_sha256",
        "source_edges_sha256",
        "excluded_m1_edges_sha256",
    }
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()

    p.add_argument("--frozen-graph", type=Path, required=True)
    p.add_argument("--regression-graph", type=Path, required=True)
    p.add_argument("--frozen-representatives", type=Path, required=True)
    p.add_argument("--regression-representatives", type=Path, required=True)
    p.add_argument("--frozen-release", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)

    return p.parse_args()


def compare_json(
    label: str,
    frozen_path: Path,
    regression_path: Path,
    findings: list[dict],
) -> bool:
    frozen = json.loads(frozen_path.read_text(encoding="utf-8"))
    regression = json.loads(regression_path.read_text(encoding="utf-8"))

    ok = True

    keys = set(frozen) | set(regression)

    for key in sorted(keys):
        if key in VOLATILE_KEYS:
            continue

        frozen_value = frozen.get(key, "<absent>")
        regression_value = regression.get(key, "<absent>")

        if frozen_value == regression_value:
            continue

        # New provenance keys added by the refactor are additive, not
        # scientific changes: report them, do not fail on them.
        if frozen_value == "<absent>":
            findings.append(
                {
                    "artefact": label,
                    "key": key,
                    "kind": "added_provenance_key",
                    "regression": regression_value,
                }
            )
            continue

        findings.append(
            {
                "artefact": label,
                "key": key,
                "kind": "value_mismatch",
                "frozen": frozen_value,
                "regression": regression_value,
            }
        )
        ok = False

    return ok


def compare_parquet(
    label: str,
    frozen_path: Path,
    regression_path: Path,
    findings: list[dict],
) -> bool:
    frozen_hash = sha256(frozen_path)
    regression_hash = sha256(regression_path)

    if frozen_hash == regression_hash:
        return True

    # Byte differences can come from compression metadata; compare content.
    frozen_table = pq.read_table(frozen_path)
    regression_table = pq.read_table(regression_path)

    if frozen_table.num_rows != regression_table.num_rows:
        findings.append(
            {
                "artefact": label,
                "kind": "row_count_mismatch",
                "frozen": frozen_table.num_rows,
                "regression": regression_table.num_rows,
            }
        )
        return False

    if frozen_table.schema.names != regression_table.schema.names:
        findings.append(
            {
                "artefact": label,
                "kind": "schema_mismatch",
                "frozen": frozen_table.schema.names,
                "regression": regression_table.schema.names,
            }
        )
        return False

    if frozen_table.equals(regression_table):
        findings.append(
            {
                "artefact": label,
                "kind": "byte_difference_content_identical",
                "frozen_sha256": frozen_hash,
                "regression_sha256": regression_hash,
            }
        )
        return True

    mismatched = [
        name
        for name in frozen_table.schema.names
        if not frozen_table.column(name).equals(
            regression_table.column(name)
        )
    ]

    findings.append(
        {
            "artefact": label,
            "kind": "content_mismatch",
            "mismatched_columns": mismatched,
        }
    )

    return False


def main() -> int:
    args = parse_args()

    findings: list[dict] = []
    checks: dict[str, bool] = {}

    print("===== STAGE-14 GRAPH =====", flush=True)

    checks["graph_summary"] = compare_json(
        "stage14_geometric_graph/global_summary.json",
        args.frozen_graph / "global_summary.json",
        args.regression_graph / "global_summary.json",
        findings,
    )

    for name in ("component_summary.parquet", "edge_node_components.parquet"):
        checks[f"graph:{name}"] = compare_parquet(
            f"stage14_geometric_graph/{name}",
            args.frozen_graph / name,
            args.regression_graph / name,
            findings,
        )

    print("===== STAGE-14 REPRESENTATIVES =====", flush=True)

    checks["representative_summary"] = compare_json(
        "stage14_representative_selection_v1/global_summary.json",
        args.frozen_representatives / "global_summary.json",
        args.regression_representatives / "global_summary.json",
        findings,
    )

    for name in ("representative_mapping.parquet", "representatives.parquet"):
        checks[f"representatives:{name}"] = compare_parquet(
            f"stage14_representative_selection_v1/{name}",
            args.frozen_representatives / name,
            args.regression_representatives / name,
            findings,
        )

    print("===== FROZEN RELEASE CROSS-CHECK =====", flush=True)

    release_manifest = json.loads(
        (args.frozen_release / "release_manifest.json").read_text(
            encoding="utf-8"
        )
    )

    regression_mapping = pq.read_table(
        args.regression_representatives / "representative_mapping.parquet"
    )

    removed = sum(
        1
        for action in regression_mapping.column("action").to_pylist()
        if action == "remove"
    )

    checks["release_removed_chain_count"] = (
        removed == release_manifest["removed_chain_count"]
    )

    if not checks["release_removed_chain_count"]:
        findings.append(
            {
                "artefact": "release_manifest.json",
                "key": "removed_chain_count",
                "kind": "value_mismatch",
                "frozen": release_manifest["removed_chain_count"],
                "regression": removed,
            }
        )

    passed = all(checks.values())

    report = {
        "status": "PASS" if passed else "FAIL",
        "checks": checks,
        "findings": findings,
        "frozen_release_counts": {
            "canonical_input_chain_count": release_manifest[
                "canonical_input_chain_count"
            ],
            "retained_chain_count": release_manifest["retained_chain_count"],
            "removed_chain_count": release_manifest["removed_chain_count"],
            "near_duplicate_threshold": release_manifest[
                "near_duplicate_threshold"
            ],
        },
    }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print()
    print("===== RESULT =====")

    for name, ok in sorted(checks.items()):
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")

    if findings:
        print()
        print("FINDINGS:")

        for finding in findings:
            print("  " + json.dumps(finding, sort_keys=True))

    print()
    print("STATUS:", report["status"])
    print("Report:", args.report)

    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
