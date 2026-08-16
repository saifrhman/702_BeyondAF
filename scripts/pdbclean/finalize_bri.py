#!/usr/bin/env python3
"""Validate and finalize a complete distributed Stage-3 BRI stage."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pyarrow.parquet as pq

from pdbclean.bri_finalize import (
    BRIFinalizeError,
    finalize_bri_stage,
)
from pdbclean.bri_runner import (
    validate_upstream_geometric_validation_stage,
)
from pdbclean.config import load_config
from pdbclean.manifest import (
    manifest_partition_count,
    resolve_manifest_snapshot,
    validate_manifest_table,
)
from pdbclean.provenance import resolve_clean_git_commit
from pdbclean.quality_runner import quality_stage_output_root


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Globally validate distributed Stage-3 Definition 3.4 "
            "BRI outputs and publish the canonical finalized "
            "BRI population."
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


def _validated_git_commit(
    value: object,
    *,
    field: str,
    path: Path,
) -> str:
    """Validate one exact Git commit identifier."""

    if (
        not isinstance(value, str)
        or len(value) != 40
        or any(
            character not in "0123456789abcdef"
            for character in value.lower()
        )
    ):
        raise BRIFinalizeError(
            f"Invalid {field} in {path}"
        )

    return value


def resolve_bri_producer_commit(
    bri_root: str | Path,
) -> str:
    """Resolve the unique Stage-3 BRI producer from task summaries."""

    summary_dir = (
        Path(bri_root)
        / "summaries"
    )

    if not summary_dir.is_dir():
        raise BRIFinalizeError(
            "Stage-3 BRI summary directory does not exist: "
            f"{summary_dir}"
        )

    summary_paths = sorted(
        summary_dir.glob("task_*.json")
    )

    if not summary_paths:
        raise BRIFinalizeError(
            "No Stage-3 BRI task summaries found in "
            f"{summary_dir}"
        )

    commits: set[str] = set()

    for path in summary_paths:
        try:
            summary = json.loads(
                path.read_text(
                    encoding="utf-8"
                )
            )
        except (
            OSError,
            UnicodeError,
            json.JSONDecodeError,
        ) as exc:
            raise BRIFinalizeError(
                "Cannot read Stage-3 BRI task summary "
                f"{path}: {exc}"
            ) from exc

        if not isinstance(summary, dict):
            raise BRIFinalizeError(
                "Stage-3 BRI task summary must contain "
                f"a JSON object: {path}"
            )

        commits.add(
            _validated_git_commit(
                summary.get(
                    "bri_pipeline_git_commit"
                ),
                field="bri_pipeline_git_commit",
                path=path,
            )
        )

    if len(commits) != 1:
        raise BRIFinalizeError(
            "Stage-3 BRI summaries contain multiple "
            "producer Git commits: "
            + ", ".join(
                sorted(commits)
            )
        )

    return next(iter(commits))


def main() -> None:
    args = parse_args()

    loaded = load_config(
        args.config
    )
    config = loaded.data

    snapshot_config = config["snapshot"]
    execution_config = config["execution"]
    storage_config = config["storage"]

    protocol_version = config[
        "release"
    ]["protocol_version"]

    manifest_path = Path(
        args.manifest
    ).resolve()

    if not manifest_path.is_file():
        raise FileNotFoundError(
            "Source manifest does not exist: "
            f"{manifest_path}"
        )

    manifest = pq.read_table(
        manifest_path
    )

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
        expected_total_bytes=(
            expected_total_bytes
        ),
    )

    batch_size = execution_config[
        "batch_size"
    ]

    task_count = manifest_partition_count(
        manifest_summary.row_count,
        batch_size,
    )

    storage_output_root = Path(
        storage_config["output_root"]
    )

    if not storage_output_root.is_absolute():
        storage_output_root = (
            REPOSITORY_ROOT
            / storage_output_root
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

    upstream = (
        validate_upstream_geometric_validation_stage(
            geometric_validation_root,
            expected_snapshot=snapshot,
            expected_cleaning_protocol=(
                protocol_version
            ),
            expected_task_count=task_count,
        )
    )

    bri_root = (
        quality_root.parent
        / "bri"
    )

    bri_pipeline_git_commit = (
        resolve_bri_producer_commit(
            bri_root
        )
    )

    # Finalizer provenance is resolved only after the immutable
    # manifest, completed Stage-2 publication, and Stage-3 producer
    # provenance have been validated.
    finalizer_pipeline_git_commit = (
        resolve_clean_git_commit(
            REPOSITORY_ROOT
        )
    )

    publication = finalize_bri_stage(
        bri_root=bri_root,
        eligible_path=upstream.eligible_path,
        manifest_row_count=(
            manifest_summary.row_count
        ),
        batch_size=batch_size,
        snapshot=snapshot,
        cleaning_protocol=(
            protocol_version
        ),
        quality_pipeline_git_commit=(
            upstream.quality_pipeline_git_commit
        ),
        geometric_validation_pipeline_git_commit=(
            upstream.geometric_validation_pipeline_git_commit
        ),
        geometric_validation_finalizer_git_commit=(
            upstream.geometric_validation_finalizer_git_commit
        ),
        bri_pipeline_git_commit=(
            bri_pipeline_git_commit
        ),
        finalizer_pipeline_git_commit=(
            finalizer_pipeline_git_commit
        ),
    )

    summary = publication.global_summary

    print(f"Snapshot: {snapshot}")
    print(
        f"Protocol: {protocol_version}"
    )
    print(
        "Quality producer Git commit: "
        f"{upstream.quality_pipeline_git_commit}"
    )
    print(
        "Geometry producer Git commit: "
        f"{upstream.geometric_validation_pipeline_git_commit}"
    )
    print(
        "Geometry finalizer Git commit: "
        f"{upstream.geometric_validation_finalizer_git_commit}"
    )
    print(
        "BRI producer Git commit: "
        f"{bri_pipeline_git_commit}"
    )
    print(
        "BRI finalizer Git commit: "
        f"{finalizer_pipeline_git_commit}"
    )
    print(
        f"Config SHA256: {loaded.sha256}"
    )
    print(
        "Manifest rows: "
        f"{manifest_summary.row_count:,}"
    )
    print(
        "Completed tasks: "
        f"{summary['task_count']:,}"
    )
    print(
        "Input eligible chains: "
        f"{summary['input_eligible_chain_count']:,}"
    )
    print(
        "BRI chains: "
        f"{summary['bri_chain_count']:,}"
    )
    print(
        "Processing errors: "
        f"{summary['processing_error_count']:,}"
    )
    print(
        "Unique BRI identities: "
        f"{summary['unique_bri_identity_count']:,}"
    )
    print(
        "Minimum retained m: "
        f"{summary['minimum_retained_residue_count']}"
    )
    print(
        "Maximum retained m: "
        f"{summary['maximum_retained_residue_count']}"
    )
    print(
        f"Finalized BRI population: "
        f"{publication.bri_path}"
    )
    print(
        f"Global summary: "
        f"{publication.global_summary_path}"
    )
    print(
        f"Stage completion: "
        f"{publication.success_path}"
    )


if __name__ == "__main__":
    main()
