#!/usr/bin/env python3
"""Run one Stage-3 Definition 3.4 BRI manifest partition."""

from __future__ import annotations

import argparse
from pathlib import Path

import pyarrow.parquet as pq

from pdbclean.bri_runner import (
    BRIRunnerError,
    execute_bri_task,
    validate_upstream_geometric_validation_stage,
)
from pdbclean.config import load_config
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
    """Parse a zero-based manifest task identifier."""

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
            "Run one zero-based Stage-3 Definition 3.4 BRI "
            "partition from the completed geometric-validation "
            "population."
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

    geometric_validation_root = (
        quality_root.parent
        / "geometric_validation"
    )

    upstream = validate_upstream_geometric_validation_stage(
        geometric_validation_root,
        expected_snapshot=snapshot,
        expected_cleaning_protocol=protocol_version,
        expected_task_count=partition_count,
    )

    # Stage-3 task membership is defined only through the immutable
    # manifest partition's source-object keys. Never depend on merged
    # eligible.parquet row-group numbering.
    source_keys = partition["s3_key"].to_pylist()

    eligible = pq.read_table(
        upstream.eligible_path,
        filters=[
            (
                "source_mmcif_key",
                "in",
                source_keys,
            )
        ],
    )

    partition_source_keys = set(source_keys)

    observed_source_keys = {
        value
        for value in eligible[
            "source_mmcif_key"
        ].to_pylist()
    }

    if not observed_source_keys.issubset(
        partition_source_keys
    ):
        raise BRIRunnerError(
            "Canonical eligible selection contains a source "
            "outside the current manifest partition"
        )

    # Production provenance is resolved only after every immutable
    # upstream input and task-selection contract has been validated.
    bri_pipeline_git_commit = resolve_clean_git_commit(
        REPOSITORY_ROOT
    )

    # Stage 3 is a sibling of quality and geometric_validation
    # under the same snapshot/protocol release root.
    output_root = (
        quality_root.parent
        / "bri"
    )

    publication = execute_bri_task(
        partition.to_pylist(),
        eligible.to_pylist(),
        upstream=upstream,
        output_root=output_root,
        task_id=args.task_id,
        bucket_url=snapshot_config["bucket_url"],
        bri_pipeline_git_commit=bri_pipeline_git_commit,
        timeout_seconds=execution_config[
            "connection_timeout_seconds"
        ],
        max_retries=execution_config[
            "max_retries"
        ],
    )

    summary = publication.summary

    print(f"Snapshot: {snapshot}")
    print(f"Protocol: {protocol_version}")
    print(
        "Quality producer Git commit: "
        f"{upstream.quality_pipeline_git_commit}"
    )
    print(
        "Geometric-validation Git commit: "
        f"{upstream.geometric_validation_pipeline_git_commit}"
    )
    print(
        "Geometric-validation finalizer Git commit: "
        f"{upstream.geometric_validation_finalizer_git_commit}"
    )
    print(
        "BRI producer Git commit: "
        f"{bri_pipeline_git_commit}"
    )
    print(
        f"Manifest rows: {manifest_summary.row_count:,}"
    )
    print(
        f"Partition: {args.task_id} / "
        f"{partition_count - 1}"
    )
    print(
        f"Partition source objects: "
        f"{partition.num_rows:,}"
    )
    print(
        f"Eligible input chains: "
        f"{eligible.num_rows:,}"
    )
    print(
        f"BRI output chains: "
        f"{summary['bri_chain_count']:,}"
    )
    print(
        f"Processing errors: "
        f"{summary['processing_error_count']:,}"
    )
    print(
        f"BRI output root: {output_root}"
    )
    print(
        f"Task summary: {publication.summary_path}"
    )


if __name__ == "__main__":
    main()
