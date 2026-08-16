"""Tests for Stage-5 Brain upstream validation."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from pdbclean.brain_runner import (
    UpstreamBRI,
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


def _upstream_for_records(
    tmp_path: Path,
) -> UpstreamBRI:
    from pdbclean.brain_runner import UpstreamBRI

    return UpstreamBRI(
        bri_root=tmp_path / "bri",
        bri_path=tmp_path / "bri" / "finalized" / "bri.parquet",
        snapshot=SNAPSHOT,
        cleaning_protocol=PROTOCOL,
        task_count=3,
        bri_chain_count=2,
        minimum_retained_residue_count=1,
        maximum_retained_residue_count=20,
        quality_pipeline_git_commit=QUALITY_COMMIT,
        geometric_validation_pipeline_git_commit=GEOMETRY_COMMIT,
        geometric_validation_finalizer_git_commit=(
            GEOMETRY_FINALIZER_COMMIT
        ),
        bri_pipeline_git_commit=BRI_COMMIT,
        bri_finalizer_git_commit=BRI_FINALIZER_COMMIT,
    )


def _stage3_bri_row(
    *,
    m: int,
    bri: list[list[float]],
) -> dict:
    row = {
        field.name: None
        for field in STAGE3_BRI_CHAIN_SCHEMA
    }

    row.update(
        {
            "snapshot": SNAPSHOT,
            "pdb_id": "1abc",
            "model_id": 1,
            "label_chain_id": "A",
            "auth_chain_id": None,
            "entity_id": None,
            "retained_residue_count": m,
            "retained_label_seq_ids": list(
                range(1, m + 1)
            ),
            "source_mmcif_key": (
                f"{SNAPSHOT}/coordinates/1abc.cif.gz"
            ),
            "source_etag": "etag",
            "cleaning_protocol": PROTOCOL,
            "quality_pipeline_git_commit": QUALITY_COMMIT,
            "geometric_validation_pipeline_git_commit": (
                GEOMETRY_COMMIT
            ),
            "geometric_validation_finalizer_git_commit": (
                GEOMETRY_FINALIZER_COMMIT
            ),
            "bri_pipeline_git_commit": BRI_COMMIT,
            "bri": bri,
        }
    )

    # Fill any additional non-null Stage-3 lineage fields with values
    # acceptable to the schema-oriented processor tests.
    for field in STAGE3_BRI_CHAIN_SCHEMA:
        if row[field.name] is not None:
            continue

        if field.nullable:
            continue

        if pa.types.is_string(field.type):
            row[field.name] = "test"
        elif pa.types.is_integer(field.type):
            row[field.name] = 1
        elif pa.types.is_boolean(field.type):
            row[field.name] = False
        elif pa.types.is_list(field.type):
            row[field.name] = [1]

    return row


def test_process_brain_record_defined_preserves_lineage_and_unrounded_mean(
    tmp_path: Path,
) -> None:
    from pdbclean.brain_runner import process_brain_record

    upstream = _upstream_for_records(
        tmp_path
    )

    bri = [
        [0.0] * 9,
        [0.001] * 9,
        [0.002] * 9,
        [0.002] * 9,
    ]

    row = _stage3_bri_row(
        m=4,
        bri=bri,
    )

    result = process_brain_record(
        row,
        upstream=upstream,
        brain_pipeline_git_commit="f" * 40,
    )

    assert result.chain_accounting_valid
    assert len(result.brain_records) == 1
    assert result.undefined_records == ()
    assert result.processing_errors == ()

    record = result.brain_records[0]

    expected = np.full(
        9,
        0.005 / 3.0,
        dtype=np.float64,
    )

    assert np.array_equal(
        np.asarray(
            record["brain"],
            dtype=np.float64,
        ),
        expected,
    )

    assert not np.array_equal(
        expected,
        np.around(expected, 3),
    )

    assert "bri" not in record
    assert record["bri_pipeline_git_commit"] == BRI_COMMIT
    assert record["bri_finalizer_git_commit"] == (
        BRI_FINALIZER_COMMIT
    )
    assert record["brain_pipeline_git_commit"] == "f" * 40


def test_process_brain_record_m1_is_explicit_undefined_not_error(
    tmp_path: Path,
) -> None:
    from pdbclean.brain_runner import (
        BRAIN_UNDEFINED_M1_REASON,
        process_brain_record,
    )

    upstream = _upstream_for_records(
        tmp_path
    )

    row = _stage3_bri_row(
        m=1,
        bri=[[0.0] * 9],
    )

    result = process_brain_record(
        row,
        upstream=upstream,
        brain_pipeline_git_commit="f" * 40,
    )

    assert result.chain_accounting_valid
    assert result.brain_records == ()
    assert result.processing_errors == ()
    assert len(result.undefined_records) == 1

    record = result.undefined_records[0]

    assert record["undefined_reason"] == (
        BRAIN_UNDEFINED_M1_REASON
    )
    assert "brain" not in record
    assert "bri" not in record
    assert record["bri_finalizer_git_commit"] == (
        BRI_FINALIZER_COMMIT
    )


def test_process_brain_record_two_residues_uses_second_row(
    tmp_path: Path,
) -> None:
    from pdbclean.brain_runner import process_brain_record

    upstream = _upstream_for_records(
        tmp_path
    )

    second = [
        0.101,
        -0.202,
        0.303,
        1.404,
        -1.505,
        2.606,
        -2.707,
        3.808,
        -3.909,
    ]

    row = _stage3_bri_row(
        m=2,
        bri=[
            [0.0] * 9,
            second,
        ],
    )

    result = process_brain_record(
        row,
        upstream=upstream,
        brain_pipeline_git_commit="f" * 40,
    )

    assert np.array_equal(
        np.asarray(
            result.brain_records[0]["brain"],
            dtype=np.float64,
        ),
        np.asarray(
            second,
            dtype=np.float64,
        ),
    )


def test_process_brain_record_rejects_lineage_mismatch_as_terminal_error(
    tmp_path: Path,
) -> None:
    from pdbclean.brain_runner import process_brain_record

    upstream = _upstream_for_records(
        tmp_path
    )

    row = _stage3_bri_row(
        m=2,
        bri=[
            [0.0] * 9,
            [1.0] * 9,
        ],
    )

    row["bri_pipeline_git_commit"] = "0" * 40

    result = process_brain_record(
        row,
        upstream=upstream,
        brain_pipeline_git_commit="f" * 40,
    )

    assert result.chain_accounting_valid
    assert result.brain_records == ()
    assert result.undefined_records == ()
    assert len(result.processing_errors) == 1

    error = result.processing_errors[0]

    assert error["processing_stage"] == (
        "bri_lineage_validation"
    )
    assert error["error_type"] == "BRILineageError"


def test_process_brain_record_invalid_bri_is_terminal_error(
    tmp_path: Path,
) -> None:
    from pdbclean.brain_runner import process_brain_record

    upstream = _upstream_for_records(
        tmp_path
    )

    row = _stage3_bri_row(
        m=2,
        bri=[
            [0.0] * 9,
            [1.0] * 9,
        ],
    )

    row["bri"][1][0] = 0.0005

    result = process_brain_record(
        row,
        upstream=upstream,
        brain_pipeline_git_commit="f" * 40,
    )

    assert result.chain_accounting_valid
    assert len(result.processing_errors) == 1
    assert result.processing_errors[0][
        "processing_stage"
    ] == "bri_input_validation"


def test_process_brain_record_computation_failure_is_terminal_error(
    tmp_path: Path,
) -> None:
    from pdbclean.brain_runner import process_brain_record

    upstream = _upstream_for_records(
        tmp_path
    )

    row = _stage3_bri_row(
        m=2,
        bri=[
            [0.0] * 9,
            [1.0] * 9,
        ],
    )

    def failing_brain(_bri):
        raise ValueError("synthetic Brain failure")

    result = process_brain_record(
        row,
        upstream=upstream,
        brain_pipeline_git_commit="f" * 40,
        brain_computer=failing_brain,
    )

    assert result.chain_accounting_valid
    assert len(result.processing_errors) == 1

    error = result.processing_errors[0]

    assert error["processing_stage"] == "brain_computation"
    assert error["error_type"] == "ValueError"
    assert error["error_message"] == "synthetic Brain failure"


def test_process_brain_record_rejects_invalid_brain_producer_commit(
    tmp_path: Path,
) -> None:
    from pdbclean.brain_runner import process_brain_record

    upstream = _upstream_for_records(
        tmp_path
    )

    row = _stage3_bri_row(
        m=1,
        bri=[[0.0] * 9],
    )

    with pytest.raises(
        BrainRunnerError,
        match="brain_pipeline_git_commit",
    ):
        process_brain_record(
            row,
            upstream=upstream,
            brain_pipeline_git_commit="not-a-commit",
        )


def test_brain_task_count_uses_ceiling_partition() -> None:
    from pdbclean.brain_runner import brain_task_count

    assert brain_task_count(
        9286,
        row_groups_per_task=64,
    ) == 146


def test_brain_task_partition_is_contiguous() -> None:
    from pdbclean.brain_runner import brain_task_partition

    counts = tuple([64] * 130)

    first = brain_task_partition(
        counts,
        task_id=0,
        row_groups_per_task=64,
    )
    second = brain_task_partition(
        counts,
        task_id=1,
        row_groups_per_task=64,
    )

    assert first.start_row_group == 0
    assert first.stop_row_group == 64
    assert first.row_group_count == 64
    assert first.input_bri_chain_count == 4096

    assert second.start_row_group == 64
    assert second.stop_row_group == 128


def test_brain_task_partition_final_task_is_short() -> None:
    from pdbclean.brain_runner import brain_task_partition

    counts = tuple([64] * 130)

    final = brain_task_partition(
        counts,
        task_id=2,
        row_groups_per_task=64,
    )

    assert final.task_count == 3
    assert final.start_row_group == 128
    assert final.stop_row_group == 130
    assert final.row_group_count == 2
    assert final.input_bri_chain_count == 128


def test_brain_task_partition_uses_actual_row_counts() -> None:
    from pdbclean.brain_runner import brain_task_partition

    counts = (
        64,
        64,
        17,
        3,
    )

    partition = brain_task_partition(
        counts,
        task_id=1,
        row_groups_per_task=2,
    )

    assert partition.start_row_group == 2
    assert partition.stop_row_group == 4
    assert partition.input_bri_chain_count == 20


def test_brain_task_partition_rejects_out_of_range_task() -> None:
    from pdbclean.brain_runner import (
        BrainRunnerError,
        brain_task_partition,
    )

    with pytest.raises(
        BrainRunnerError,
        match="outside the physical partition range",
    ):
        brain_task_partition(
            (64, 64),
            task_id=1,
            row_groups_per_task=64,
        )


def test_brain_task_partition_rejects_invalid_row_group_counts() -> None:
    from pdbclean.brain_runner import (
        BrainRunnerError,
        brain_task_partition,
    )

    with pytest.raises(
        BrainRunnerError,
        match="row-group row count",
    ):
        brain_task_partition(
            (64, 0, 64),
            task_id=0,
            row_groups_per_task=64,
        )
