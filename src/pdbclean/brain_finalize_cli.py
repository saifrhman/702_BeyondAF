"""CLI for Stage-5 Brain global finalization."""

from __future__ import annotations

import argparse
from pathlib import Path

import pyarrow.parquet as pq

from pdbclean.brain_finalize import (
    finalize_brain_stage,
)
from pdbclean.brain_runner import (
    validate_upstream_bri_stage,
)
from pdbclean.config import load_config
from pdbclean.manifest import (
    manifest_partition_count,
    resolve_manifest_snapshot,
    validate_manifest_table,
)
from pdbclean.quality_runner import (
    quality_stage_output_root,
)


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--brain-pipeline-git-commit",
        required=True,
    )
    parser.add_argument(
        "--finalizer-pipeline-git-commit",
        required=True,
    )

    args = parser.parse_args()

    repo = Path.cwd()
    loaded = load_config(
        args.config
    )
    config = loaded.data

    protocol = (
        config["release"]["protocol_version"]
    )
    snapshot_config = config["snapshot"]
    batch_size = (
        config["execution"]["batch_size"]
    )

    completed = []

    for manifest_path in sorted(
        (repo / "outputs/pdbclean").glob(
            "*/bronze/source_manifest.parquet"
        )
    ):
        candidate_snapshot = (
            manifest_path.parents[1].name
        )

        candidate_bri = (
            repo
            / "outputs/pdbclean"
            / candidate_snapshot
            / protocol
            / "bri"
        )

        if (
            candidate_bri / "_SUCCESS"
        ).is_file():
            completed.append(
                manifest_path
            )

    if len(completed) != 1:
        raise RuntimeError(
            "Expected exactly one completed Stage-3 "
            f"publication; found {len(completed)}"
        )

    manifest = pq.read_table(
        completed[0]
    )

    snapshot = resolve_manifest_snapshot(
        manifest,
        snapshot_config,
    )

    manifest_summary = validate_manifest_table(
        manifest,
        expected_snapshot=snapshot,
    )

    expected_stage3_tasks = (
        manifest_partition_count(
            manifest_summary.row_count,
            batch_size,
        )
    )

    storage_root = Path(
        config["storage"]["output_root"]
    )

    if not storage_root.is_absolute():
        storage_root = (
            repo / storage_root
        )

    quality_root = quality_stage_output_root(
        storage_root,
        snapshot=snapshot,
        protocol_version=protocol,
    )

    bri_root = (
        quality_root.parent / "bri"
    )

    upstream = validate_upstream_bri_stage(
        bri_root,
        expected_snapshot=snapshot,
        expected_cleaning_protocol=protocol,
        expected_task_count=(
            expected_stage3_tasks
        ),
    )

    publication = finalize_brain_stage(
        brain_root=(
            bri_root.parent / "brain"
        ),
        upstream=upstream,
        brain_pipeline_git_commit=(
            args.brain_pipeline_git_commit
        ),
        finalizer_pipeline_git_commit=(
            args.finalizer_pipeline_git_commit
        ),
    )

    summary = publication.global_summary

    print(
        "Input BRI chains:",
        summary["input_bri_chain_count"],
    )
    print(
        "Brain-defined:",
        summary["brain_chain_count"],
    )
    print(
        "Brain-undefined:",
        summary["undefined_chain_count"],
    )
    print(
        "Processing errors:",
        summary["processing_error_count"],
    )
    print(
        "Unique terminal identities:",
        summary["unique_terminal_identity_count"],
    )
    print(
        "Accounting:",
        summary["chain_accounting_valid"],
    )
    print()
    print(
        "STAGE-5 BRAIN GLOBAL FINALIZATION: PASS"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
