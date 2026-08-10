#!/usr/bin/env python3
"""Validate and merge a complete distributed PDBClean quality stage."""

from __future__ import annotations

import argparse
from pathlib import Path

import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.manifest import (
    resolve_manifest_snapshot,
    validate_manifest_table,
)
from pdbclean.provenance import resolve_clean_git_commit
from pdbclean.quality_merge import merge_quality_stage
from pdbclean.quality_runner import quality_stage_output_root


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate all distributed PDBClean quality tasks, merge "
            "their Gold outputs, and publish the quality-stage "
            "_SUCCESS marker."
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

    pipeline_git_commit = resolve_clean_git_commit(
        REPOSITORY_ROOT
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

    publication = merge_quality_stage(
        quality_root=quality_root,
        manifest=manifest,
        manifest_row_count=manifest_summary.row_count,
        batch_size=execution_config["batch_size"],
        snapshot=snapshot,
        cleaning_protocol=protocol_version,
        pipeline_git_commit=pipeline_git_commit,
    )

    print(f"Snapshot: {snapshot}")
    print(f"Protocol: {protocol_version}")
    print(f"Git commit: {pipeline_git_commit}")
    print(f"Config SHA256: {loaded.sha256}")
    print(f"Manifest rows: {manifest_summary.row_count:,}")
    print(f"Quality root: {quality_root}")
    print(
        f"Completed tasks: "
        f"{publication.global_summary['task_count']:,}"
    )
    print(
        f"Accepted chains: "
        f"{publication.global_summary['accepted_chain_count']:,}"
    )
    print(
        f"Rejected chains: "
        f"{publication.global_summary['rejected_chain_count']:,}"
    )
    print(
        f"Non-candidate chains: "
        f"{publication.global_summary['non_candidate_chain_count']:,}"
    )
    print(
        f"Processing errors: "
        f"{publication.global_summary['processing_error_count']:,}"
    )
    print(
        f"Global summary: {publication.global_summary_path}"
    )
    print(f"Stage completion: {publication.success_path}")


if __name__ == "__main__":
    main()
