#!/usr/bin/env python3
"""Validate and merge a complete distributed PDBClean quality stage."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.manifest import (
    resolve_manifest_snapshot,
    validate_manifest_table,
)
from pdbclean.provenance import resolve_clean_git_commit
from pdbclean.quality_merge import (
    QualityMergeError,
    merge_quality_stage,
)
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


def resolve_quality_pipeline_git_commit(
    quality_root: str | Path,
) -> str:
    """Resolve the unique Git commit that produced quality tasks."""

    summary_dir = Path(quality_root) / "summaries"

    if not summary_dir.is_dir():
        raise QualityMergeError(
            f"Quality summary directory does not exist: "
            f"{summary_dir}"
        )

    summary_paths = sorted(
        summary_dir.glob("task_*.json")
    )

    if not summary_paths:
        raise QualityMergeError(
            f"No quality-task summaries found in {summary_dir}"
        )

    commits: set[str] = set()

    for summary_path in summary_paths:
        try:
            summary = json.loads(
                summary_path.read_text(encoding="utf-8")
            )
        except (
            OSError,
            UnicodeError,
            json.JSONDecodeError,
        ) as exc:
            raise QualityMergeError(
                f"Cannot read quality-task summary "
                f"{summary_path}: {exc}"
            ) from exc

        if not isinstance(summary, dict):
            raise QualityMergeError(
                f"Quality-task summary must contain a "
                f"JSON object: {summary_path}"
            )

        commit = summary.get("pipeline_git_commit")

        if (
            not isinstance(commit, str)
            or len(commit) != 40
            or any(
                character not in "0123456789abcdef"
                for character in commit.lower()
            )
        ):
            raise QualityMergeError(
                f"Invalid pipeline_git_commit in "
                f"{summary_path}"
            )

        commits.add(commit)

    if len(commits) != 1:
        raise QualityMergeError(
            "Quality-task summaries were produced by "
            "multiple Git commits: "
            + ", ".join(sorted(commits))
        )

    return next(iter(commits))


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

    merge_pipeline_git_commit = resolve_clean_git_commit(
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

    quality_pipeline_git_commit = (
        resolve_quality_pipeline_git_commit(
            quality_root
        )
    )

    publication = merge_quality_stage(
        quality_root=quality_root,
        manifest=manifest,
        manifest_row_count=manifest_summary.row_count,
        batch_size=execution_config["batch_size"],
        snapshot=snapshot,
        cleaning_protocol=protocol_version,
        quality_pipeline_git_commit=(
            quality_pipeline_git_commit
        ),
        merge_pipeline_git_commit=(
            merge_pipeline_git_commit
        ),
    )

    print(f"Snapshot: {snapshot}")
    print(f"Protocol: {protocol_version}")
    print(
        f"Quality producer Git commit: "
        f"{quality_pipeline_git_commit}"
    )
    print(
        f"Merge Git commit: "
        f"{merge_pipeline_git_commit}"
    )
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
