"""Tests for Stage-3 BRI production orchestration."""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from pdbclean.bri_runner import (
    BRIRunnerError,
    validate_upstream_geometric_validation_stage,
)
from pdbclean.geometric_validation_finalize import (
    GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_NAME,
    GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_VERSION,
    GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_NAME,
    GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_VERSION,
)
from pdbclean.schemas import STAGE3_ELIGIBLE_CHAIN_SCHEMA


SNAPSHOT = "20310415"
PROTOCOL = "protocol3.2-comp702-v1"

QUALITY_COMMIT = "1" * 40
GEOMETRY_COMMIT = "2" * 40
FINALIZER_COMMIT = "3" * 40


def _eligible_row() -> dict:
    return {
        "snapshot": SNAPSHOT,
        "pdb_id": "test",
        "model_id": 1,
        "label_chain_id": "A",
        "auth_chain_id": "A",
        "entity_id": "1",
        "original_start_label_seq_id": 1,
        "original_end_label_seq_id": 1,
        "retained_start_label_seq_id": 1,
        "retained_end_label_seq_id": 1,
        "retained_residue_count": 1,
        "retained_label_seq_ids": [1],
        "retained_sequence": "A",
        "terminal_trimmed": False,
        "dirty_residue_count": 0,
        "dirty_rule_ids": [],
        "source_mmcif_key": (
            f"{SNAPSHOT}/pub/pdb/data/structures/"
            "divided/mmCIF/es/test.cif.gz"
        ),
        "source_etag": "etag-1",
        "cleaning_protocol": PROTOCOL,
        "pipeline_git_commit": QUALITY_COMMIT,
    }


def _success() -> dict:
    return {
        "success_schema_name": (
            GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_NAME
        ),
        "success_schema_version": (
            GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_VERSION
        ),
        "snapshot": SNAPSHOT,
        "cleaning_protocol": PROTOCOL,
        "quality_pipeline_git_commit": QUALITY_COMMIT,
        "geometric_validation_pipeline_git_commit": (
            GEOMETRY_COMMIT
        ),
        "finalizer_pipeline_git_commit": FINALIZER_COMMIT,
        "task_count": 4,
        "global_summary": "global_summary.json",
        "finalized_directory": "finalized",
        "eligible_population": "finalized/eligible.parquet",
        "quarantined_population": "finalized/quarantined.parquet",
    }


def _summary() -> dict:
    return {
        "summary_schema_name": (
            GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_NAME
        ),
        "summary_schema_version": (
            GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_VERSION
        ),
        "snapshot": SNAPSHOT,
        "cleaning_protocol": PROTOCOL,
        "quality_pipeline_git_commit": QUALITY_COMMIT,
        "geometric_validation_pipeline_git_commit": (
            GEOMETRY_COMMIT
        ),
        "finalizer_pipeline_git_commit": FINALIZER_COMMIT,
        "task_count": 4,
        "input_accepted_chain_count": 2,
        "eligible_chain_count": 1,
        "quarantined_chain_count": 1,
        "processing_error_count": 0,
        "chain_accounting_valid": True,
    }


def _write_stage(
    tmp_path: Path,
    *,
    success: dict | None = None,
    summary: dict | None = None,
    eligible_rows: list[dict] | None = None,
    eligible_schema: pa.Schema = STAGE3_ELIGIBLE_CHAIN_SCHEMA,
) -> Path:
    root = tmp_path / "geometric_validation"
    finalized = root / "finalized"
    finalized.mkdir(parents=True)

    success = _success() if success is None else success
    summary = _summary() if summary is None else summary

    if eligible_rows is None:
        eligible_rows = [_eligible_row()]

    (root / "_SUCCESS").write_text(
        json.dumps(success) + "\n",
        encoding="utf-8",
    )

    (root / "global_summary.json").write_text(
        json.dumps(summary) + "\n",
        encoding="utf-8",
    )

    pq.write_table(
        pa.Table.from_pylist(
            eligible_rows,
            schema=eligible_schema,
        ),
        finalized / "eligible.parquet",
    )

    return root


def _validate(root: Path):
    return validate_upstream_geometric_validation_stage(
        root,
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_task_count=4,
    )


def test_validate_upstream_geometric_validation_stage(
    tmp_path: Path,
) -> None:
    root = _write_stage(tmp_path)

    observed = _validate(root)

    assert observed.snapshot == SNAPSHOT
    assert observed.cleaning_protocol == PROTOCOL
    assert observed.task_count == 4
    assert observed.eligible_chain_count == 1

    assert observed.quality_pipeline_git_commit == QUALITY_COMMIT
    assert (
        observed.geometric_validation_pipeline_git_commit
        == GEOMETRY_COMMIT
    )
    assert (
        observed.geometric_validation_finalizer_git_commit
        == FINALIZER_COMMIT
    )

    assert observed.eligible_path == (
        root / "finalized" / "eligible.parquet"
    )


def test_upstream_requires_success_marker(tmp_path: Path) -> None:
    root = _write_stage(tmp_path)
    (root / "_SUCCESS").unlink()

    with pytest.raises(
        BRIRunnerError,
        match="Stage-2 _SUCCESS marker does not exist",
    ):
        _validate(root)


def test_upstream_rejects_success_schema_mismatch(
    tmp_path: Path,
) -> None:
    success = _success()
    success["success_schema_version"] = "999"

    root = _write_stage(
        tmp_path,
        success=success,
    )

    with pytest.raises(
        BRIRunnerError,
        match="Unexpected Stage-2 _SUCCESS schema",
    ):
        _validate(root)


def test_upstream_rejects_provenance_mismatch(
    tmp_path: Path,
) -> None:
    summary = _summary()
    summary["geometric_validation_pipeline_git_commit"] = "4" * 40

    root = _write_stage(
        tmp_path,
        summary=summary,
    )

    with pytest.raises(
        BRIRunnerError,
        match=(
            "Stage-2 provenance mismatch: "
            "geometric_validation_pipeline_git_commit"
        ),
    ):
        _validate(root)


def test_upstream_rejects_processing_errors(
    tmp_path: Path,
) -> None:
    summary = _summary()
    summary["processing_error_count"] = 1
    summary["quarantined_chain_count"] = 0

    root = _write_stage(
        tmp_path,
        summary=summary,
    )

    with pytest.raises(
        BRIRunnerError,
        match="contains processing errors",
    ):
        _validate(root)


def test_upstream_rejects_eligible_row_count_mismatch(
    tmp_path: Path,
) -> None:
    summary = _summary()
    summary["input_accepted_chain_count"] = 3
    summary["eligible_chain_count"] = 2

    root = _write_stage(
        tmp_path,
        summary=summary,
    )

    with pytest.raises(
        BRIRunnerError,
        match="eligible row count does not match",
    ):
        _validate(root)


def test_upstream_rejects_eligible_schema_mismatch(
    tmp_path: Path,
) -> None:
    incompatible_schema = pa.schema(
        [
            pa.field(
                "pdb_id",
                pa.string(),
                nullable=False,
            )
        ],
        metadata={
            b"schema_name": b"wrong",
            b"schema_version": b"1.0",
        },
    )

    root = _write_stage(
        tmp_path,
        eligible_rows=[{"pdb_id": "test"}],
        eligible_schema=incompatible_schema,
    )

    with pytest.raises(
        BRIRunnerError,
        match="eligible population schema mismatch",
    ):
        _validate(root)
