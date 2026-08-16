"""Command-line entrypoint for one Stage-5 Brain production task."""

from __future__ import annotations

import argparse
from pathlib import Path

import pyarrow.parquet as pq

from pdbclean.brain_runner import (
    DEFAULT_BRAIN_ROW_GROUPS_PER_TASK,
    brain_task_partition,
    execute_brain_task,
    validate_upstream_bri_stage,
)
from pdbclean.config import load_config
from pdbclean.manifest import (
    manifest_partition_count,
    resolve_manifest_snapshot,
    validate_manifest_table,
)
from pdbclean.quality_runner import quality_stage_output_root


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--task-id",
        required=True,
        type=int,
    )
    parser.add_argument(
        "--brain-pipeline-git-commit",
        required=True,
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--row-groups-per-task",
        type=int,
        default=DEFAULT_BRAIN_ROW_GROUPS_PER_TASK,
    )

    args = parser.parse_args()

    loaded = load_config(args.config)
    config = loaded.data

    repo = Path.cwd()

    protocol = config["release"]["protocol_version"]
    snapshot_config = config["snapshot"]
    batch_size = config["execution"]["batch_size"]

    manifests = sorted(
        (repo / "outputs/pdbclean").glob(
            "*/bronze/source_manifest.parquet"
        )
    )

    completed = []

    for path in manifests:
        snapshot_candidate = path.parents[1].name

        candidate_root = (
            repo
            / "outputs/pdbclean"
            / snapshot_candidate
            / protocol
            / "bri"
        )

        if (candidate_root / "_SUCCESS").is_file():
            completed.append(path)

    if len(completed) != 1:
        raise RuntimeError(
            "Expected exactly one manifest associated with a "
            "completed Stage-3 BRI publication; found "
            f"{[str(path) for path in completed]!r}"
        )

    manifest_path = completed[0]
    manifest = pq.read_table(manifest_path)

    snapshot = resolve_manifest_snapshot(
        manifest,
        snapshot_config,
    )

    manifest_summary = validate_manifest_table(
        manifest,
        expected_snapshot=snapshot,
    )

    expected_stage3_task_count = manifest_partition_count(
        manifest_summary.row_count,
        batch_size,
    )

    storage_output_root = Path(
        config["storage"]["output_root"]
    )

    if not storage_output_root.is_absolute():
        storage_output_root = repo / storage_output_root

    quality_root = quality_stage_output_root(
        storage_output_root,
        snapshot=snapshot,
        protocol_version=protocol,
    )

    bri_root = quality_root.parent / "bri"

    upstream = validate_upstream_bri_stage(
        bri_root,
        expected_snapshot=snapshot,
        expected_cleaning_protocol=protocol,
        expected_task_count=expected_stage3_task_count,
    )

    parquet = pq.ParquetFile(
        upstream.bri_path
    )

    row_group_counts = tuple(
        parquet.metadata.row_group(i).num_rows
        for i in range(
            parquet.metadata.num_row_groups
        )
    )

    partition = brain_task_partition(
        row_group_counts,
        task_id=args.task_id,
        row_groups_per_task=(
            args.row_groups_per_task
        ),
    )

    output_root = args.output_root

    if output_root is None:
        output_root = bri_root.parent / "brain"

    publication = execute_brain_task(
        upstream=upstream,
        partition=partition,
        output_root=output_root,
        brain_pipeline_git_commit=(
            args.brain_pipeline_git_commit
        ),
    )

    print("Stage-5 Brain task:", partition.task_id)
    print("Task count:", partition.task_count)
    print(
        "Row groups:",
        f"{partition.start_row_group}.."
        f"{partition.stop_row_group - 1}",
    )
    print(
        "Input BRI chains:",
        partition.input_bri_chain_count,
    )
    print(
        "Brain-defined:",
        publication.summary["brain_chain_count"],
    )
    print(
        "Brain-undefined:",
        publication.summary["undefined_chain_count"],
    )
    print(
        "Processing errors:",
        publication.summary["processing_error_count"],
    )
    print(
        "Accounting:",
        publication.summary["chain_accounting_valid"],
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
