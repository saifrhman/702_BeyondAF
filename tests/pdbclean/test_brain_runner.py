"""Tests for Stage-5 Brain upstream validation."""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from pdbclean.brain_runner import (
    BrainRunnerError,
    validate_upstream_bri_stage,
)
from pdbclean.schemas import STAGE3_BRI_CHAIN_SCHEMA


SNAPSHOT = "20310415"
PROTOCOL = "protocol3.2-comp702-v1"

QUALITY_COMMIT = "a" * 40
GEOMETRY_COMMIT = "b" * 40
GEOMETRY_FINALIZER_COMMIT = "c" * 40
BRI_COMMIT = "d" * 40
BRI_FINALIZER_COMMIT = "e" * 40


def _write_completed_bri(
    root: Path,
    *,
    success_updates: dict | None = None,
    summary_updates: dict | None = None,
) -> Path:
    finalized = root / "finalized"
    finalized.mkdir(parents=True)

    bri_path = finalized / "bri.parquet"

    # One canonical-schema row is unnecessary for these contract tests;
    # an empty Parquet file is written and its metadata row count is
    # overridden through the published summary below only when desired.
    pq.write_table(
        STAGE3_BRI_CHAIN_SCHEMA.empty_table(),
        bri_path,
    )

    success = {
        "success_schema_name": "pdbclean_stage3_bri_success",
        "success_schema_version": "1.0",
        "snapshot": SNAPSHOT,
        "cleaning_protocol": PROTOCOL,
        "quality_pipeline_git_commit": QUALITY_COMMIT,
        "geometric_validation_pipeline_git_commit": GEOMETRY_COMMIT,
        "geometric_validation_finalizer_git_commit": (
            GEOMETRY_FINALIZER_COMMIT
        ),
        "bri_pipeline_git_commit": BRI_COMMIT,
        "finalizer_pipeline_git_commit": BRI_FINALIZER_COMMIT,
        "task_count": 3,
        "global_summary": "global_summary.json",
        "finalized_directory": "finalized",
        "bri_population": "finalized/bri.parquet",
    }

    summary = {
        "summary_schema_name": "pdbclean_stage3_bri_global_summary",
        "summary_schema_version": "1.0",
        "snapshot": SNAPSHOT,
        "cleaning_protocol": PROTOCOL,
        "quality_pipeline_git_commit": QUALITY_COMMIT,
        "geometric_validation_pipeline_git_commit": GEOMETRY_COMMIT,
        "geometric_validation_finalizer_git_commit": (
            GEOMETRY_FINALIZER_COMMIT
        ),
        "bri_pipeline_git_commit": BRI_COMMIT,
        "finalizer_pipeline_git_commit": BRI_FINALIZER_COMMIT,
        "task_count": 3,
        "manifest_source_object_count": 5,
        "relevant_source_object_count": 2,
        "downloaded_source_object_count": 2,
        "parsed_source_object_count": 2,
        "input_eligible_chain_count": 1,
        "bri_chain_count": 1,
        "processing_error_count": 0,
        "unique_eligible_identity_count": 1,
        "unique_bri_identity_count": 1,
        "minimum_retained_residue_count": 1,
        "maximum_retained_residue_count": 20,
        "chain_accounting_valid": True,
    }

    if success_updates:
        success.update(success_updates)

    if summary_updates:
        summary.update(summary_updates)

    (root / "_SUCCESS").write_text(
        json.dumps(success) + "\n",
        encoding="utf-8",
    )

    (root / "global_summary.json").write_text(
        json.dumps(summary) + "\n",
        encoding="utf-8",
    )

    return bri_path


def _patch_row_count(
    monkeypatch: pytest.MonkeyPatch,
    *,
    row_count: int,
) -> None:
    real_parquet_file = pq.ParquetFile

    class _Metadata:
        num_rows = row_count

    class _ParquetFile:
        metadata = _Metadata()

    def fake_parquet_file(path):
        # Ensure the validator still receives the expected path.
        assert Path(path).name == "bri.parquet"
        return _ParquetFile()

    monkeypatch.setattr(
        "pdbclean.brain_runner.pq.ParquetFile",
        fake_parquet_file,
    )


def test_upstream_bri_accepts_completed_canonical_publication(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "bri"

    _write_completed_bri(root)
    _patch_row_count(
        monkeypatch,
        row_count=1,
    )

    result = validate_upstream_bri_stage(
        root,
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_task_count=3,
    )

    assert result.bri_path == root / "finalized" / "bri.parquet"
    assert result.task_count == 3
    assert result.bri_chain_count == 1
    assert result.minimum_retained_residue_count == 1
    assert result.maximum_retained_residue_count == 20
    assert result.quality_pipeline_git_commit == QUALITY_COMMIT
    assert result.bri_pipeline_git_commit == BRI_COMMIT
    assert result.bri_finalizer_git_commit == BRI_FINALIZER_COMMIT


def test_upstream_bri_rejects_noncanonical_completion_pointer(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "bri"

    _write_completed_bri(
        root,
        success_updates={
            "bri_population": "other.parquet",
        },
    )
    _patch_row_count(monkeypatch, row_count=1)

    with pytest.raises(
        BrainRunnerError,
        match="completion pointer",
    ):
        validate_upstream_bri_stage(
            root,
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
        )


def test_upstream_bri_rejects_provenance_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "bri"

    _write_completed_bri(
        root,
        summary_updates={
            "bri_pipeline_git_commit": "f" * 40,
        },
    )
    _patch_row_count(monkeypatch, row_count=1)

    with pytest.raises(
        BrainRunnerError,
        match="provenance mismatch",
    ):
        validate_upstream_bri_stage(
            root,
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
        )


def test_upstream_bri_rejects_processing_errors(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "bri"

    _write_completed_bri(
        root,
        summary_updates={
            "input_eligible_chain_count": 2,
            "bri_chain_count": 1,
            "processing_error_count": 1,
            "unique_eligible_identity_count": 2,
            "unique_bri_identity_count": 1,
        },
    )
    _patch_row_count(monkeypatch, row_count=1)

    with pytest.raises(
        BrainRunnerError,
        match="processing errors",
    ):
        validate_upstream_bri_stage(
            root,
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
        )


def test_upstream_bri_rejects_task_count_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "bri"

    _write_completed_bri(root)
    _patch_row_count(monkeypatch, row_count=1)

    with pytest.raises(
        BrainRunnerError,
        match="manifest partition contract",
    ):
        validate_upstream_bri_stage(
            root,
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
            expected_task_count=4,
        )


def test_upstream_bri_rejects_schema_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "bri"

    _write_completed_bri(root)
    _patch_row_count(monkeypatch, row_count=1)

    monkeypatch.setattr(
        "pdbclean.brain_runner.pq.read_schema",
        lambda path: STAGE3_BRI_CHAIN_SCHEMA.remove(
            STAGE3_BRI_CHAIN_SCHEMA.get_field_index("bri")
        ),
    )

    with pytest.raises(
        BrainRunnerError,
        match="schema mismatch",
    ):
        validate_upstream_bri_stage(
            root,
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
        )


def test_upstream_bri_rejects_finalized_row_count_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "bri"

    _write_completed_bri(root)
    _patch_row_count(
        monkeypatch,
        row_count=2,
    )

    with pytest.raises(
        BrainRunnerError,
        match="row count",
    ):
        validate_upstream_bri_stage(
            root,
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
        )
