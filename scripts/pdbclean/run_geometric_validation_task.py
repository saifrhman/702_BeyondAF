#!/usr/bin/env python3
"""Run one post-cleaning geometric-validation manifest partition."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.geometric_validation import (
    GeometricValidationConfig,
)
from pdbclean.geometric_validation_runner import (
    GeometricValidationRunnerError,
    execute_geometric_validation_task,
    validate_upstream_quality_task,
)
from pdbclean.manifest import (
    manifest_partition_count,
    resolve_manifest_snapshot,
    select_manifest_partition,
    validate_manifest_table,
)
from pdbclean.provenance import resolve_clean_git_commit
from pdbclean.quality_runner import quality_stage_output_root


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def _non_negative_integer(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "task ID must be an integer"
        ) from exc

    if parsed < 0:
        raise argparse.ArgumentTypeError(
            "task ID must be non-negative"
        )

    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run one zero-based post-cleaning geometric-validation "
            "partition from completed PDBClean quality outputs."
        )
    )

    parser.add_argument(
        "--config",
        required=True,
        help="Path to the versioned PDBClean YAML configuration.",
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to the immutable source_manifest.parquet.",
    )
    parser.add_argument(
        "--task-id",
        required=True,
        type=_non_negative_integer,
        help="Zero-based manifest partition identifier.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    loaded = load_config(args.config)
    config = loaded.data

    snapshot_config = config["snapshot"]
    execution_config = config["execution"]
    storage_config = config["storage"]
    geometry_config = config[
        "post_cleaning_geometric_validation"
    ]

    if geometry_config["enabled"] is not True:
        raise GeometricValidationRunnerError(
            "Post-cleaning geometric validation is disabled "
            "by configuration"
        )

    protocol_version = config["release"]["protocol_version"]

    manifest_path = Path(args.manifest).resolve()

    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"Source manifest does not exist: {manifest_path}"
        )

    manifest = pq.read_table(manifest_path)

    snapshot = resolve_manifest_snapshot(
        manifest,
        snapshot_config,
    )

    expected_count = None
    expected_total_bytes = None

    if snapshot_config["mode"] == "fixed":
        expected_count = snapshot_config.get(
            "expected_mmcif_count"
        )
        expected_total_bytes = snapshot_config.get(
            "expected_total_bytes"
        )

    manifest_summary = validate_manifest_table(
        manifest,
        expected_snapshot=snapshot,
        expected_count=expected_count,
        expected_total_bytes=expected_total_bytes,
    )

    batch_size = execution_config["batch_size"]

    partition_count = manifest_partition_count(
        manifest_summary.row_count,
        batch_size,
    )

    partition = select_manifest_partition(
        manifest,
        task_id=args.task_id,
        batch_size=batch_size,
    )

    storage_output_root = Path(
        storage_config["output_root"]
    )

    if not storage_output_root.is_absolute():
        storage_output_root = (
            REPOSITORY_ROOT / storage_output_root
        )

    quality_root = quality_stage_output_root(
        storage_output_root,
        snapshot=snapshot,
        protocol_version=protocol_version,
    )

    task_id_text = str(args.task_id)

    quality_summary_path = (
        quality_root
        / "summaries"
        / f"task_{task_id_text}.json"
    )
    accepted_path = (
        quality_root
        / "accepted"
        / f"task_{task_id_text}.parquet"
    )

    if not quality_summary_path.is_file():
        raise FileNotFoundError(
            "Upstream quality-task completion summary does not "
            f"exist: {quality_summary_path}"
        )

    if not accepted_path.is_file():
        raise FileNotFoundError(
            "Upstream accepted-chain shard does not exist: "
            f"{accepted_path}"
        )

    quality_summary = json.loads(
        quality_summary_path.read_text(
            encoding="utf-8"
        )
    )

    accepted = pq.read_table(accepted_path)

    upstream = validate_upstream_quality_task(
        quality_summary,
        accepted.to_pylist(),
        expected_task_id=args.task_id,
        expected_snapshot=snapshot,
        expected_cleaning_protocol=protocol_version,
        expected_manifest_source_object_count=(
            partition.num_rows
        ),
    )

    validation_config = GeometricValidationConfig(
        minimum_backbone_distance_angstrom=(
            config["quality_rules"]
            ["backbone_distance"]
            ["minimum_distance_angstrom"]
        ),
        minimum_triangle_angle_degrees=(
            geometry_config[
                "minimum_triangle_angle_degrees"
            ]
        ),
    )

    geometric_validation_pipeline_git_commit = (
        resolve_clean_git_commit(
            REPOSITORY_ROOT
        )
    )

    # Step 2 is a sibling stage of quality under the same
    # snapshot/protocol release root.
    output_root = (
        quality_root.parent
        / "geometric_validation"
    )

    publication = execute_geometric_validation_task(
        partition.to_pylist(),
        upstream.accepted_rows,
        output_root=output_root,
        task_id=args.task_id,
        snapshot=snapshot,
        bucket_url=snapshot_config["bucket_url"],
        config=validation_config,
        cleaning_protocol=protocol_version,
        quality_pipeline_git_commit=(
            upstream.quality_pipeline_git_commit
        ),
        geometric_validation_pipeline_git_commit=(
            geometric_validation_pipeline_git_commit
        ),
        timeout_seconds=execution_config[
            "connection_timeout_seconds"
        ],
        max_retries=execution_config["max_retries"],
        download_concurrency=execution_config[
            "download_concurrency"
        ],
    )

    print(f"Snapshot: {snapshot}")
    print(f"Protocol: {protocol_version}")
    print(
        "Quality producer Git commit: "
        f"{upstream.quality_pipeline_git_commit}"
    )
    print(
        "Geometric-validation Git commit: "
        f"{geometric_validation_pipeline_git_commit}"
    )
    print(
        f"Manifest rows: {manifest_summary.row_count:,}"
    )
    print(
        f"Partition: {args.task_id} / "
        f"{partition_count - 1}"
    )
    print(
        f"Partition rows: {partition.num_rows:,}"
    )
    print(
        f"Accepted input chains: "
        f"{len(upstream.accepted_rows):,}"
    )
    print(
        f"Geometric-validation output root: {output_root}"
    )
    print(
        f"Task summary: {publication.summary_path}"
    )


if __name__ == "__main__":
    main()
