#!/usr/bin/env python3
"""Validate and finalize a complete distributed Step-2 geometry stage."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.geometric_validation_finalize import (
    GeometricValidationFinalizeError,
    finalize_geometric_validation_stage,
)
from pdbclean.manifest import (
    resolve_manifest_snapshot,
    validate_manifest_table,
)
from pdbclean.provenance import resolve_clean_git_commit
from pdbclean.quality_runner import quality_stage_output_root


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Globally validate post-cleaning geometric-validation "
            "outputs and publish the canonical Stage-3 eligible "
            "and quarantined populations."
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
    if (
        not isinstance(value, str)
        or len(value) != 40
        or any(
            character not in "0123456789abcdef"
            for character in value.lower()
        )
    ):
        raise GeometricValidationFinalizeError(
            f"Invalid {field} in {path}"
        )

    return value


def resolve_stage2_producer_commits(
    geometric_validation_root: str | Path,
) -> tuple[str, str]:
    """Resolve unique Step-1 and Step-2 producer commits from summaries."""

    summary_dir = (
        Path(geometric_validation_root)
        / "summaries"
    )

    if not summary_dir.is_dir():
        raise GeometricValidationFinalizeError(
            "Geometric-validation summary directory does not "
            f"exist: {summary_dir}"
        )

    summary_paths = sorted(
        summary_dir.glob("task_*.json")
    )

    if not summary_paths:
        raise GeometricValidationFinalizeError(
            f"No Step-2 task summaries found in {summary_dir}"
        )

    quality_commits: set[str] = set()
    geometry_commits: set[str] = set()

    for path in summary_paths:
        try:
            summary = json.loads(
                path.read_text(encoding="utf-8")
            )
        except (
            OSError,
            UnicodeError,
            json.JSONDecodeError,
        ) as exc:
            raise GeometricValidationFinalizeError(
                f"Cannot read Step-2 task summary {path}: {exc}"
            ) from exc

        if not isinstance(summary, dict):
            raise GeometricValidationFinalizeError(
                f"Step-2 task summary must contain a JSON object: "
                f"{path}"
            )

        quality_commits.add(
            _validated_git_commit(
                summary.get(
                    "quality_pipeline_git_commit"
                ),
                field="quality_pipeline_git_commit",
                path=path,
            )
        )

        geometry_commits.add(
            _validated_git_commit(
                summary.get(
                    "geometric_validation_pipeline_git_commit"
                ),
                field=(
                    "geometric_validation_pipeline_git_commit"
                ),
                path=path,
            )
        )

    if len(quality_commits) != 1:
        raise GeometricValidationFinalizeError(
            "Step-2 summaries contain multiple quality producer "
            "Git commits: "
            + ", ".join(sorted(quality_commits))
        )

    if len(geometry_commits) != 1:
        raise GeometricValidationFinalizeError(
            "Step-2 summaries contain multiple geometry producer "
            "Git commits: "
            + ", ".join(sorted(geometry_commits))
        )

    return (
        next(iter(quality_commits)),
        next(iter(geometry_commits)),
    )


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
        raise GeometricValidationFinalizeError(
            "Post-cleaning geometric validation is disabled "
            "by configuration"
        )

    protocol_version = config[
        "release"
    ]["protocol_version"]

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

    # Stage 2 was derived from a globally completed quality stage.
    quality_success = quality_root / "_SUCCESS"

    if not quality_success.is_file():
        raise GeometricValidationFinalizeError(
            "Upstream quality stage has no _SUCCESS marker: "
            f"{quality_success}"
        )

    (
        quality_pipeline_git_commit,
        geometric_validation_pipeline_git_commit,
    ) = resolve_stage2_producer_commits(
        geometric_validation_root
    )

    # Final publication itself must also have exact clean-tree provenance.
    finalizer_pipeline_git_commit = (
        resolve_clean_git_commit(
            REPOSITORY_ROOT
        )
    )

    minimum_backbone_distance_angstrom = (
        config["quality_rules"]
        ["backbone_distance"]
        ["minimum_distance_angstrom"]
    )

    minimum_triangle_angle_degrees = (
        geometry_config[
            "minimum_triangle_angle_degrees"
        ]
    )

    publication = finalize_geometric_validation_stage(
        quality_root=quality_root,
        geometric_validation_root=(
            geometric_validation_root
        ),
        manifest_row_count=manifest_summary.row_count,
        batch_size=execution_config["batch_size"],
        snapshot=snapshot,
        cleaning_protocol=protocol_version,
        quality_pipeline_git_commit=(
            quality_pipeline_git_commit
        ),
        geometric_validation_pipeline_git_commit=(
            geometric_validation_pipeline_git_commit
        ),
        finalizer_pipeline_git_commit=(
            finalizer_pipeline_git_commit
        ),
        minimum_backbone_distance_angstrom=(
            minimum_backbone_distance_angstrom
        ),
        minimum_triangle_angle_degrees=(
            minimum_triangle_angle_degrees
        ),
    )

    summary = publication.global_summary

    print(f"Snapshot: {snapshot}")
    print(f"Protocol: {protocol_version}")
    print(
        "Quality producer Git commit: "
        f"{quality_pipeline_git_commit}"
    )
    print(
        "Geometry producer Git commit: "
        f"{geometric_validation_pipeline_git_commit}"
    )
    print(
        "Finalizer Git commit: "
        f"{finalizer_pipeline_git_commit}"
    )
    print(f"Config SHA256: {loaded.sha256}")
    print(
        f"Manifest rows: {manifest_summary.row_count:,}"
    )
    print(
        f"Completed tasks: {summary['task_count']:,}"
    )
    print(
        "Input accepted chains: "
        f"{summary['input_accepted_chain_count']:,}"
    )
    print(
        "Eligible chains: "
        f"{summary['eligible_chain_count']:,}"
    )
    print(
        "Quarantined chains: "
        f"{summary['quarantined_chain_count']:,}"
    )
    print(
        "Processing errors: "
        f"{summary['processing_error_count']:,}"
    )
    print(
        "Violation events: "
        f"{summary['violation_event_count']:,}"
    )
    print(
        f"Eligible population: {publication.eligible_path}"
    )
    print(
        "Quarantined population: "
        f"{publication.quarantined_path}"
    )
    print(
        f"Global summary: {publication.global_summary_path}"
    )
    print(
        f"Stage completion: {publication.success_path}"
    )


if __name__ == "__main__":
    main()
