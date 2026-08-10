#!/usr/bin/env python3
"""Run one deterministic PDBClean quality-cleaning manifest partition."""

from __future__ import annotations

import argparse
from pathlib import Path

import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.manifest import (
    manifest_partition_count,
    resolve_manifest_snapshot,
    select_manifest_partition,
    validate_manifest_table,
)
from pdbclean.provenance import resolve_clean_git_commit
from pdbclean.quality_runner import (
    execute_quality_task,
    quality_stage_output_root,
)


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
            "Run one zero-based quality-cleaning partition from an "
            "immutable PDBClean source manifest."
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

    # Resolve the actual snapshot represented by the immutable manifest.
    # In fixed mode this must match snapshot.snapshot_id. In
    # latest_complete mode the manifest itself supplies the resolved ID.
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

    # Production provenance is always the full clean HEAD SHA.
    pipeline_git_commit = resolve_clean_git_commit(
        REPOSITORY_ROOT
    )

    storage_output_root = Path(
        storage_config["output_root"]
    )

    # Relative configured output paths are anchored to the repository,
    # never to an arbitrary SLURM working directory.
    if not storage_output_root.is_absolute():
        storage_output_root = (
            REPOSITORY_ROOT / storage_output_root
        )

    output_root = quality_stage_output_root(
        storage_output_root,
        snapshot=snapshot,
        protocol_version=protocol_version,
    )

    publication = execute_quality_task(
        partition.to_pylist(),
        output_root=output_root,
        task_id=args.task_id,
        snapshot=snapshot,
        bucket_url=snapshot_config["bucket_url"],
        selection_config=config["selection"],
        cleaning_protocol=protocol_version,
        pipeline_git_commit=pipeline_git_commit,
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
    print(f"Git commit: {pipeline_git_commit}")
    print(
        f"Manifest rows: {manifest_summary.row_count:,}"
    )
    print(
        f"Partition: {args.task_id} / "
        f"{partition_count - 1}"
    )
    print(f"Partition rows: {partition.num_rows:,}")
    print(f"Quality output root: {output_root}")
    print(f"Task summary: {publication.summary_path}")


if __name__ == "__main__":
    main()
