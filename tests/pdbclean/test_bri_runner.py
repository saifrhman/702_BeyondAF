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


def _source_upstream() -> UpstreamGeometricValidation:
    from pdbclean.bri_runner import UpstreamGeometricValidation

    return UpstreamGeometricValidation(
        geometric_validation_root=Path("/unused/geometry"),
        eligible_path=Path("/unused/eligible.parquet"),
        snapshot=SNAPSHOT,
        cleaning_protocol=PROTOCOL,
        task_count=4,
        eligible_chain_count=1,
        quality_pipeline_git_commit=QUALITY_COMMIT,
        geometric_validation_pipeline_git_commit=GEOMETRY_COMMIT,
        geometric_validation_finalizer_git_commit=FINALIZER_COMMIT,
    )


def _source_manifest() -> dict:
    return {
        "snapshot": SNAPSHOT,
        "pdb_id": "TEST",
        "s3_key": (
            f"{SNAPSHOT}/pub/pdb/data/structures/"
            "divided/mmCIF/es/test.cif.gz"
        ),
        "size_bytes": 123,
        "etag": "etag-1",
    }


def _source_chain(
    *,
    residue_id: int = 1,
) -> ChainObservation:
    from pdbclean.mmcif_parser import (
        AtomObservation,
        ChainObservation,
    )

    def atom(
        name: str,
        xyz: tuple[float, float, float],
    ) -> AtomObservation:
        return AtomObservation(
            model_id=1,
            label_chain_id="A",
            auth_chain_id="A",
            entity_id="1",
            label_seq_id=residue_id,
            auth_seq_id=str(residue_id),
            residue_name="ALA",
            atom_name=name,
            alt_id=None,
            occupancy=1.0,
            x=xyz[0],
            y=xyz[1],
            z=xyz[2],
            group_pdb="ATOM",
            occupancy_raw="1.00",
        )

    return ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        auth_chain_id="A",
        entity_id="1",
        polymer_type="polypeptide(L)",
        entry_has_polypeptide=True,
        atoms=[
            atom("N", (0.0, 0.0, 0.0)),
            atom("CA", (1.0, 0.0, 0.0)),
            atom("C", (1.0, 1.0, 0.0)),
        ],
    )


def _source_eligible(
    *,
    retained_ids: list[int] | None = None,
    source_etag: str = "etag-1",
) -> dict:
    row = _eligible_row().copy()

    if retained_ids is not None:
        row["retained_label_seq_ids"] = retained_ids
        row["retained_residue_count"] = len(retained_ids)
        row["retained_start_label_seq_id"] = retained_ids[0]
        row["retained_end_label_seq_id"] = retained_ids[-1]

    row["source_etag"] = source_etag
    return row


def _process_source(
    rows,
    *,
    downloader=None,
    parser=None,
    bri_computer=None,
    max_retries: int = 0,
):
    from pdbclean.bri_runner import process_bri_source
    from pdbclean.snapshot import download_verified_s3_object_bytes
    from pdbclean.mmcif_parser import parse_coordinate_mmcif_bytes
    from pdbclean.bri import compute_bri

    return process_bri_source(
        _source_manifest(),
        rows,
        upstream=_source_upstream(),
        bucket_url="https://example.invalid",
        bri_pipeline_git_commit="4" * 40,
        max_retries=max_retries,
        downloader=(
            download_verified_s3_object_bytes
            if downloader is None
            else downloader
        ),
        parser=(
            parse_coordinate_mmcif_bytes
            if parser is None
            else parser
        ),
        bri_computer=(
            compute_bri
            if bri_computer is None
            else bri_computer
        ),
    )


def test_process_bri_source_success() -> None:
    result = _process_source(
        [_source_eligible()],
        downloader=lambda **kwargs: b"source",
        parser=lambda *args, **kwargs: [_source_chain()],
    )

    assert result.chain_accounting_valid is True
    assert result.input_eligible_chain_count == 1
    assert len(result.bri_records) == 1
    assert len(result.processing_errors) == 0
    assert result.source_downloaded is True
    assert result.source_parsed is True

    record = result.bri_records[0]

    assert record["pdb_id"] == "test"
    assert record["retained_residue_count"] == 1
    assert record["quality_pipeline_git_commit"] == QUALITY_COMMIT
    assert (
        record["geometric_validation_pipeline_git_commit"]
        == GEOMETRY_COMMIT
    )
    assert (
        record["geometric_validation_finalizer_git_commit"]
        == FINALIZER_COMMIT
    )
    assert record["bri_pipeline_git_commit"] == "4" * 40

    assert record["bri"] == [
        [
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
        ]
    ]


def test_process_bri_source_with_no_eligible_rows_does_not_download() -> None:
    def forbidden_download(**kwargs):
        raise AssertionError("Downloader must not be called")

    result = _process_source(
        [],
        downloader=forbidden_download,
    )

    assert result.chain_accounting_valid is True
    assert result.input_eligible_chain_count == 0
    assert result.bri_records == ()
    assert result.processing_errors == ()
    assert result.source_downloaded is False
    assert result.source_parsed is False


def test_process_bri_source_rejects_lineage_before_download() -> None:
    def forbidden_download(**kwargs):
        raise AssertionError("Downloader must not be called")

    result = _process_source(
        [
            _source_eligible(
                source_etag="wrong-etag",
            )
        ],
        downloader=forbidden_download,
    )

    assert result.chain_accounting_valid is True
    assert len(result.bri_records) == 0
    assert len(result.processing_errors) == 1

    error = result.processing_errors[0]

    assert (
        error["processing_stage"]
        == "eligible_lineage_validation"
    )
    assert error["error_type"] == "EligibleLineageError"


def test_process_bri_source_retries_transport_failure() -> None:
    from pdbclean.snapshot import SnapshotTransportError

    attempts = []

    def downloader(**kwargs):
        attempts.append(1)

        if len(attempts) == 1:
            raise SnapshotTransportError("temporary")

        return b"source"

    result = _process_source(
        [_source_eligible()],
        downloader=downloader,
        parser=lambda *args, **kwargs: [_source_chain()],
        max_retries=1,
    )

    assert len(attempts) == 2
    assert result.chain_accounting_valid is True
    assert len(result.bri_records) == 1


def test_process_bri_source_parse_error_is_terminal_chain_error() -> None:
    from pdbclean.mmcif_parser import MMCIFParseError

    def parser(*args, **kwargs):
        raise MMCIFParseError("broken source")

    result = _process_source(
        [_source_eligible()],
        downloader=lambda **kwargs: b"source",
        parser=parser,
    )

    assert result.chain_accounting_valid is True
    assert len(result.bri_records) == 0
    assert len(result.processing_errors) == 1
    assert (
        result.processing_errors[0]["processing_stage"]
        == "mmcif_parse"
    )


def test_process_bri_source_missing_chain_is_terminal_error() -> None:
    result = _process_source(
        [_source_eligible()],
        downloader=lambda **kwargs: b"source",
        parser=lambda *args, **kwargs: [],
    )

    assert result.chain_accounting_valid is True
    assert len(result.bri_records) == 0
    assert len(result.processing_errors) == 1
    assert (
        result.processing_errors[0]["processing_stage"]
        == "source_chain_lookup"
    )


def test_process_bri_source_reconstruction_failure_is_terminal_error() -> None:
    result = _process_source(
        [_source_eligible(retained_ids=[2])],
        downloader=lambda **kwargs: b"source",
        parser=lambda *args, **kwargs: [_source_chain()],
    )

    assert result.chain_accounting_valid is True
    assert len(result.bri_records) == 0
    assert len(result.processing_errors) == 1
    assert (
        result.processing_errors[0]["processing_stage"]
        == "retained_chain_reconstruction"
    )


def test_process_bri_source_bri_failure_is_terminal_error() -> None:
    def bri_failure(chain):
        raise ValueError("undefined BRI")

    result = _process_source(
        [_source_eligible()],
        downloader=lambda **kwargs: b"source",
        parser=lambda *args, **kwargs: [_source_chain()],
        bri_computer=bri_failure,
    )

    assert result.chain_accounting_valid is True
    assert len(result.bri_records) == 0
    assert len(result.processing_errors) == 1
    assert (
        result.processing_errors[0]["processing_stage"]
        == "bri_computation"
    )


def _successful_source_result(
    rows,
) -> SourceBRIResult:
    from pdbclean.bri_runner import SourceBRIResult

    single = _process_source(
        rows,
        downloader=lambda **kwargs: b"source",
        parser=lambda *args, **kwargs: [
            _source_chain()
        ],
    )

    return single


def test_process_bri_batch_materializes_frozen_schemas() -> None:
    from pdbclean.bri_runner import process_bri_batch
    from pdbclean.schemas import (
        STAGE3_BRI_CHAIN_SCHEMA,
        STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
    )

    row = _source_eligible()

    def processor(
        manifest_row,
        rows,
        **kwargs,
    ):
        return _successful_source_result(rows)

    batch = process_bri_batch(
        [_source_manifest()],
        [row],
        upstream=_source_upstream(),
        bucket_url="https://example.invalid",
        bri_pipeline_git_commit="4" * 40,
        source_processor=processor,
    )

    assert batch.chain_accounting_valid is True
    assert batch.manifest_source_object_count == 1
    assert batch.relevant_source_object_count == 1
    assert batch.input_eligible_chain_count == 1
    assert batch.downloaded_source_object_count == 1
    assert batch.parsed_source_object_count == 1

    assert batch.tables.chains.num_rows == 1
    assert batch.tables.processing_errors.num_rows == 0

    assert batch.tables.chains.schema.equals(
        STAGE3_BRI_CHAIN_SCHEMA,
        check_metadata=True,
    )

    assert batch.tables.processing_errors.schema.equals(
        STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
        check_metadata=True,
    )


def test_process_bri_batch_skips_irrelevant_manifest_sources() -> None:
    from pdbclean.bri_runner import process_bri_batch

    second = _source_manifest().copy()
    second["pdb_id"] = "NONE"
    second["s3_key"] = (
        f"{SNAPSHOT}/none.cif.gz"
    )
    second["etag"] = "etag-2"

    calls = []

    def processor(
        manifest_row,
        rows,
        **kwargs,
    ):
        calls.append(manifest_row["s3_key"])
        return _successful_source_result(rows)

    batch = process_bri_batch(
        [_source_manifest(), second],
        [_source_eligible()],
        upstream=_source_upstream(),
        bucket_url="https://example.invalid",
        bri_pipeline_git_commit="4" * 40,
        source_processor=processor,
    )

    assert len(calls) == 1
    assert batch.manifest_source_object_count == 2
    assert batch.relevant_source_object_count == 1


def test_process_bri_batch_rejects_source_outside_partition() -> None:
    from pdbclean.bri_runner import (
        BRIRunnerError,
        process_bri_batch,
    )

    row = _source_eligible()
    row["source_mmcif_key"] = (
        f"{SNAPSHOT}/outside.cif.gz"
    )

    with pytest.raises(
        BRIRunnerError,
        match="outside the current manifest partition",
    ):
        process_bri_batch(
            [_source_manifest()],
            [row],
            upstream=_source_upstream(),
            bucket_url="https://example.invalid",
            bri_pipeline_git_commit="4" * 40,
        )


def test_process_bri_batch_rejects_duplicate_chain_identity() -> None:
    from pdbclean.bri_runner import (
        BRIRunnerError,
        process_bri_batch,
    )

    row = _source_eligible()

    with pytest.raises(
        BRIRunnerError,
        match="duplicate eligible chain identity",
    ):
        process_bri_batch(
            [_source_manifest()],
            [row, row.copy()],
            upstream=_source_upstream(),
            bucket_url="https://example.invalid",
            bri_pipeline_git_commit="4" * 40,
        )


def test_write_bri_shards_round_trip(
    tmp_path: Path,
) -> None:
    import pyarrow.parquet as pq

    from pdbclean.bri_runner import (
        process_bri_batch,
        write_bri_shards,
    )
    from pdbclean.schemas import (
        STAGE3_BRI_CHAIN_SCHEMA,
        STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
    )

    def processor(
        manifest_row,
        rows,
        **kwargs,
    ):
        return _successful_source_result(rows)

    batch = process_bri_batch(
        [_source_manifest()],
        [_source_eligible()],
        upstream=_source_upstream(),
        bucket_url="https://example.invalid",
        bri_pipeline_git_commit="4" * 40,
        source_processor=processor,
    )

    paths = write_bri_shards(
        batch.tables,
        tmp_path,
        "7",
    )

    assert paths["chains"] == (
        tmp_path / "chains" / "task_7.parquet"
    )

    assert paths["processing_errors"] == (
        tmp_path
        / "processing_errors"
        / "task_7.parquet"
    )

    assert pq.read_table(
        paths["chains"]
    ).num_rows == 1

    assert pq.read_table(
        paths["processing_errors"]
    ).num_rows == 0

    # Parquet may rename list child fields, so compare the
    # in-memory canonical tables rather than requiring strict
    # nested-child-name equality after Parquet round-trip.
    assert batch.tables.chains.schema.equals(
        STAGE3_BRI_CHAIN_SCHEMA,
        check_metadata=True,
    )

    assert batch.tables.processing_errors.schema.equals(
        STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
        check_metadata=True,
    )


def test_publish_bri_batch_writes_summary_last(
    tmp_path: Path,
) -> None:
    from pdbclean.bri_runner import (
        BRITaskContext,
        process_bri_batch,
        publish_bri_batch,
    )

    def processor(
        manifest_row,
        rows,
        **kwargs,
    ):
        return _successful_source_result(rows)

    batch = process_bri_batch(
        [_source_manifest()],
        [_source_eligible()],
        upstream=_source_upstream(),
        bucket_url="https://example.invalid",
        bri_pipeline_git_commit="4" * 40,
        source_processor=processor,
    )

    events = []

    def shard_writer(
        tables,
        output_root,
        task_id,
    ):
        events.append("shards")

        return {
            "chains": tmp_path / "chains.parquet",
            "processing_errors": (
                tmp_path / "errors.parquet"
            ),
        }

    def summary_writer(
        summary,
        output_root,
    ):
        events.append("summary")
        return tmp_path / "summary.json"

    context = BRITaskContext(
        task_id="7",
        snapshot=SNAPSHOT,
        cleaning_protocol=PROTOCOL,
        quality_pipeline_git_commit=QUALITY_COMMIT,
        geometric_validation_pipeline_git_commit=(
            GEOMETRY_COMMIT
        ),
        geometric_validation_finalizer_git_commit=(
            FINALIZER_COMMIT
        ),
        bri_pipeline_git_commit="4" * 40,
        started_at_utc="2031-04-15T00:00:00Z",
    )

    publication = publish_bri_batch(
        batch,
        output_root=tmp_path,
        context=context,
        started_perf_counter=10.0,
        shard_writer=shard_writer,
        summary_writer=summary_writer,
        utc_now=lambda: "2031-04-15T00:00:01Z",
        perf_counter=lambda: 11.5,
        peak_memory_reader=lambda: 1234,
    )

    assert events == ["shards", "summary"]

    summary = publication.summary

    assert summary["task_id"] == "7"
    assert summary["input_eligible_chain_count"] == 1
    assert summary["bri_chain_count"] == 1
    assert summary["processing_error_count"] == 0
    assert summary["chain_accounting_valid"] is True
    assert summary["runtime_seconds"] == 1.5
    assert summary["peak_memory_bytes"] == 1234


def test_publish_bri_batch_invalidates_stale_summary(
    tmp_path: Path,
) -> None:
    from pdbclean.bri_runner import (
        BRIRunnerError,
        BRITaskContext,
        process_bri_batch,
        publish_bri_batch,
    )

    def processor(
        manifest_row,
        rows,
        **kwargs,
    ):
        return _successful_source_result(rows)

    batch = process_bri_batch(
        [_source_manifest()],
        [_source_eligible()],
        upstream=_source_upstream(),
        bucket_url="https://example.invalid",
        bri_pipeline_git_commit="4" * 40,
        source_processor=processor,
    )

    stale = (
        tmp_path
        / "summaries"
        / "task_7.json"
    )

    stale.parent.mkdir(parents=True)
    stale.write_text(
        "stale\n",
        encoding="utf-8",
    )

    context = BRITaskContext(
        task_id="7",
        snapshot=SNAPSHOT,
        cleaning_protocol=PROTOCOL,
        quality_pipeline_git_commit=QUALITY_COMMIT,
        geometric_validation_pipeline_git_commit=(
            GEOMETRY_COMMIT
        ),
        geometric_validation_finalizer_git_commit=(
            FINALIZER_COMMIT
        ),
        bri_pipeline_git_commit="4" * 40,
        started_at_utc="2031-04-15T00:00:00Z",
    )

    def failing_shard_writer(*args, **kwargs):
        raise BRIRunnerError("publication failed")

    with pytest.raises(
        BRIRunnerError,
        match="publication failed",
    ):
        publish_bri_batch(
            batch,
            output_root=tmp_path,
            context=context,
            started_perf_counter=1.0,
            shard_writer=failing_shard_writer,
        )

    assert not stale.exists()


def test_execute_bri_task_propagates_validated_provenance(
    tmp_path: Path,
) -> None:
    from pdbclean.bri_runner import (
        BRIBatchResult,
        BRITables,
        BRITaskPublication,
        execute_bri_task,
    )
    from pdbclean.schemas import (
        STAGE3_BRI_CHAIN_SCHEMA,
        STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
    )

    empty_tables = BRITables(
        chains=pa.Table.from_pylist(
            [],
            schema=STAGE3_BRI_CHAIN_SCHEMA,
        ),
        processing_errors=pa.Table.from_pylist(
            [],
            schema=STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
        ),
    )

    batch = BRIBatchResult(
        manifest_source_object_count=1,
        relevant_source_object_count=0,
        input_eligible_chain_count=0,
        downloaded_source_object_count=0,
        parsed_source_object_count=0,
        tables=empty_tables,
    )

    captured = {}

    def batch_processor(
        manifest_rows,
        eligible_rows,
        **kwargs,
    ):
        return batch

    def publisher(
        observed_batch,
        *,
        output_root,
        context,
        started_perf_counter,
        **kwargs,
    ):
        captured["context"] = context

        return BRITaskPublication(
            shard_paths={},
            summary_path=tmp_path / "summary.json",
            summary={},
        )

    execute_bri_task(
        [_source_manifest()],
        [],
        upstream=_source_upstream(),
        output_root=tmp_path,
        task_id=7,
        bucket_url="https://example.invalid",
        bri_pipeline_git_commit="4" * 40,
        environ={
            "SLURM_JOB_ID": "100",
            "SLURM_ARRAY_TASK_ID": "7",
        },
        batch_processor=batch_processor,
        publisher=publisher,
        utc_now=lambda: "2031-04-15T00:00:00Z",
        perf_counter=lambda: 10.0,
    )

    context = captured["context"]

    assert context.task_id == "7"
    assert context.snapshot == SNAPSHOT
    assert context.quality_pipeline_git_commit == QUALITY_COMMIT
    assert (
        context.geometric_validation_pipeline_git_commit
        == GEOMETRY_COMMIT
    )
    assert (
        context.geometric_validation_finalizer_git_commit
        == FINALIZER_COMMIT
    )
    assert context.bri_pipeline_git_commit == "4" * 40
    assert context.slurm_job_id == "100"
    assert context.slurm_array_task_id == "7"
