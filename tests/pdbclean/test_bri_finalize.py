"""Tests for Stage-3 BRI task-artifact finalization validation."""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from pdbclean.bri_finalize import (
    BRIFinalizeError,
    discover_bri_task_artifacts,
    validate_bri_task_accounting,
)
from pdbclean.schemas import (
    STAGE3_BRI_CHAIN_SCHEMA,
    STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
)


SNAPSHOT = "20310415"
PROTOCOL = "protocol3.2-comp702-v1"

QUALITY_COMMIT = "a" * 40
GEOMETRY_COMMIT = "b" * 40
GEOMETRY_FINALIZER_COMMIT = "c" * 40
BRI_COMMIT = "d" * 40


def _summary(
    task_id: int,
    *,
    input_count: int = 0,
    bri_count: int = 0,
    error_count: int = 0,
) -> dict:
    return {
        "summary_schema_name":
            "pdbclean_stage3_bri_task_summary",
        "summary_schema_version": "1.0",
        "task_id": str(task_id),
        "snapshot": SNAPSHOT,
        "cleaning_protocol": PROTOCOL,
        "quality_pipeline_git_commit":
            QUALITY_COMMIT,
        "geometric_validation_pipeline_git_commit":
            GEOMETRY_COMMIT,
        "geometric_validation_finalizer_git_commit":
            GEOMETRY_FINALIZER_COMMIT,
        "bri_pipeline_git_commit":
            BRI_COMMIT,
        "started_at_utc":
            "2031-04-15T00:00:00.000000Z",
        "completed_at_utc":
            "2031-04-15T00:00:01.000000Z",
        "runtime_seconds": 1.0,
        "slurm_job_id": None,
        "slurm_array_task_id": None,
        "peak_memory_bytes": None,
        "manifest_source_object_count": 0,
        "relevant_source_object_count": 0,
        "downloaded_source_object_count": 0,
        "parsed_source_object_count": 0,
        "input_eligible_chain_count":
            input_count,
        "bri_chain_count":
            bri_count,
        "processing_error_count":
            error_count,
        "processing_errors_by_stage": {},
        "processing_errors_by_type": {},
        "chain_accounting_valid": True,
    }


def _write_task(
    bri_root: Path,
    task_id: int,
    *,
    summary: dict | None = None,
) -> None:
    chains_path = (
        bri_root
        / "chains"
        / f"task_{task_id}.parquet"
    )

    errors_path = (
        bri_root
        / "processing_errors"
        / f"task_{task_id}.parquet"
    )

    summary_path = (
        bri_root
        / "summaries"
        / f"task_{task_id}.json"
    )

    for path in (
        chains_path,
        errors_path,
        summary_path,
    ):
        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    pq.write_table(
        pa.Table.from_pylist(
            [],
            schema=STAGE3_BRI_CHAIN_SCHEMA,
        ),
        chains_path,
    )

    pq.write_table(
        pa.Table.from_pylist(
            [],
            schema=STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
        ),
        errors_path,
    )

    if summary is None:
        summary = _summary(task_id)

    summary_path.write_text(
        json.dumps(summary) + "\n",
        encoding="utf-8",
    )


def _discover(
    bri_root: Path,
    *,
    expected_task_ids: tuple[int, ...],
):
    return discover_bri_task_artifacts(
        bri_root,
        expected_task_ids=expected_task_ids,
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_quality_pipeline_git_commit=(
            QUALITY_COMMIT
        ),
        expected_geometric_validation_pipeline_git_commit=(
            GEOMETRY_COMMIT
        ),
        expected_geometric_validation_finalizer_git_commit=(
            GEOMETRY_FINALIZER_COMMIT
        ),
        expected_bri_pipeline_git_commit=(
            BRI_COMMIT
        ),
    )


def test_discovery_accepts_complete_task_set(
    tmp_path: Path,
) -> None:
    bri_root = tmp_path / "bri"

    _write_task(bri_root, 0)
    _write_task(bri_root, 1)

    artifacts = _discover(
        bri_root,
        expected_task_ids=(0, 1),
    )

    assert tuple(
        artifact.task_id
        for artifact in artifacts
    ) == (0, 1)


def test_discovery_rejects_missing_error_shard(
    tmp_path: Path,
) -> None:
    bri_root = tmp_path / "bri"

    _write_task(bri_root, 0)

    (
        bri_root
        / "processing_errors"
        / "task_0.parquet"
    ).unlink()

    with pytest.raises(
        BRIFinalizeError,
        match=(
            "BRI processing-error shard"
            ".*missing=\\[0\\]"
        ),
    ):
        _discover(
            bri_root,
            expected_task_ids=(0,),
        )


def test_discovery_rejects_temporary_file(
    tmp_path: Path,
) -> None:
    bri_root = tmp_path / "bri"

    _write_task(bri_root, 0)

    temporary = (
        bri_root
        / "chains"
        / "task_0.parquet.tmp"
    )

    temporary.write_text(
        "partial",
        encoding="utf-8",
    )

    with pytest.raises(
        BRIFinalizeError,
        match="Temporary Stage-3 BRI files remain",
    ):
        _discover(
            bri_root,
            expected_task_ids=(0,),
        )


def test_discovery_rejects_bri_provenance_mismatch(
    tmp_path: Path,
) -> None:
    bri_root = tmp_path / "bri"

    summary = _summary(0)
    summary["bri_pipeline_git_commit"] = (
        "e" * 40
    )

    _write_task(
        bri_root,
        0,
        summary=summary,
    )

    with pytest.raises(
        BRIFinalizeError,
        match="BRI producer Git commit mismatch",
    ):
        _discover(
            bri_root,
            expected_task_ids=(0,),
        )


def test_discovery_rejects_wrong_chain_schema(
    tmp_path: Path,
) -> None:
    bri_root = tmp_path / "bri"

    _write_task(bri_root, 0)

    chains_path = (
        bri_root
        / "chains"
        / "task_0.parquet"
    )

    pq.write_table(
        pa.table(
            {
                "unexpected": [1],
            }
        ),
        chains_path,
    )

    with pytest.raises(
        BRIFinalizeError,
        match="Unexpected Parquet",
    ):
        _discover(
            bri_root,
            expected_task_ids=(0,),
        )


def test_task_accounting_recomputes_bri_row_count(
    tmp_path: Path,
) -> None:
    bri_root = tmp_path / "bri"

    summary = _summary(
        0,
        input_count=1,
        bri_count=1,
        error_count=0,
    )

    _write_task(
        bri_root,
        0,
        summary=summary,
    )

    artifacts = _discover(
        bri_root,
        expected_task_ids=(0,),
    )

    with pytest.raises(
        BRIFinalizeError,
        match="BRI chain row-count mismatch",
    ):
        validate_bri_task_accounting(
            artifacts
        )


def test_task_accounting_recomputes_terminal_invariant(
    tmp_path: Path,
) -> None:
    bri_root = tmp_path / "bri"

    summary = _summary(
        0,
        input_count=1,
        bri_count=0,
        error_count=0,
    )

    _write_task(
        bri_root,
        0,
        summary=summary,
    )

    artifacts = _discover(
        bri_root,
        expected_task_ids=(0,),
    )

    with pytest.raises(
        BRIFinalizeError,
        match="Terminal chain accounting mismatch",
    ):
        validate_bri_task_accounting(
            artifacts
        )
