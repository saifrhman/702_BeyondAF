"""Tests for source-level post-cleaning geometric validation."""

from __future__ import annotations

import json

import pyarrow.parquet as pq
import pytest

from pdbclean.geometric_validation import GeometricValidationConfig
from pdbclean.geometric_validation_runner import (
    GeometricValidationBatchResult,
    GeometricValidationRunnerError,
    GeometricValidationTaskContext,
    SourceGeometricValidationResult,
    build_geometric_validation_task_summary,
    execute_geometric_validation_task,
    geometric_validation_records_to_tables,
    process_geometric_validation_batch,
    process_geometric_validation_source,
    publish_geometric_validation_batch,
    validate_upstream_quality_task,
    write_geometric_validation_shards,
    write_geometric_validation_task_summary_atomic,
)
from pdbclean.mmcif_parser import (
    AtomObservation,
    ChainObservation,
    MMCIFParseError,
)
from pdbclean.schemas import (
    GEOMETRIC_VALIDATION_AUDIT_SCHEMA,
    GEOMETRIC_VALIDATION_PROCESSING_ERROR_SCHEMA,
)


def _atom(
    chain_id: str,
    residue_id: int,
    atom_name: str,
    xyz: tuple[float, float, float],
) -> AtomObservation:
    return AtomObservation(
        model_id=1,
        label_chain_id=chain_id,
        auth_chain_id=chain_id,
        entity_id="1",
        label_seq_id=residue_id,
        auth_seq_id=str(residue_id),
        residue_name="ALA",
        atom_name=atom_name,
        alt_id=None,
        occupancy=1.0,
        x=xyz[0],
        y=xyz[1],
        z=xyz[2],
        group_pdb="ATOM",
        occupancy_raw="1.00",
    )


def _chain(
    chain_id: str,
    *,
    collinear: bool = False,
) -> ChainObservation:
    if collinear:
        n = (0.0, 0.0, 0.0)
        ca = (1.0, 0.0, 0.0)
        c = (2.0, 0.0, 0.0)
    else:
        n = (0.0, 0.0, 0.0)
        ca = (1.0, 0.0, 0.0)
        c = (1.0, 1.0, 0.0)

    return ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id=chain_id,
        auth_chain_id=chain_id,
        entity_id="1",
        polymer_type="polypeptide(L)",
        entry_has_polypeptide=True,
        atoms=[
            _atom(chain_id, 1, "N", n),
            _atom(chain_id, 1, "CA", ca),
            _atom(chain_id, 1, "C", c),
        ],
    )


def _manifest() -> dict:
    return {
        "snapshot": "20310415",
        "pdb_id": "TEST",
        "s3_key": "20310415/coordinates/test.cif.gz",
        "size_bytes": 123,
        "etag": "etag-1",
    }


def _accepted(
    chain_id: str,
    *,
    etag: str = "etag-1",
) -> dict:
    return {
        "snapshot": "20310415",
        "pdb_id": "test",
        "model_id": 1,
        "label_chain_id": chain_id,
        "retained_residue_count": 1,
        "retained_label_seq_ids": [1],
        "source_mmcif_key": (
            "20310415/coordinates/test.cif.gz"
        ),
        "source_etag": etag,
        "cleaning_protocol": "protocol3.2-comp702-v1",
        "pipeline_git_commit": "q" * 40,
    }


def _process(
    accepted_rows,
    *,
    parser,
    downloader=None,
):
    if downloader is None:
        downloader = lambda **kwargs: b"compressed"

    return process_geometric_validation_source(
        _manifest(),
        accepted_rows,
        bucket_url="https://example.invalid",
        config=GeometricValidationConfig(),
        geometric_validation_pipeline_git_commit="g" * 40,
        downloader=downloader,
        parser=parser,
    )


def test_source_is_downloaded_and_parsed_once_for_multiple_chains() -> None:
    calls = {
        "download": 0,
        "parse": 0,
    }

    def downloader(**kwargs):
        calls["download"] += 1
        return b"compressed"

    def parser(compressed_bytes, *, pdb_id):
        calls["parse"] += 1
        assert compressed_bytes == b"compressed"
        assert pdb_id == "test"
        return [
            _chain("A"),
            _chain("B"),
        ]

    result = _process(
        [
            _accepted("A"),
            _accepted("B"),
        ],
        downloader=downloader,
        parser=parser,
    )

    assert calls == {
        "download": 1,
        "parse": 1,
    }

    assert result.input_accepted_chain_count == 2
    assert len(result.audit_records) == 2
    assert result.processing_errors == ()
    assert result.chain_accounting_valid is True
    assert result.source_downloaded is True
    assert result.source_parsed is True

    assert [
        row["label_chain_id"]
        for row in result.audit_records
    ] == ["A", "B"]

    assert all(
        row["passed"] is True
        for row in result.audit_records
    )


def test_lineage_mismatch_becomes_chain_error_without_download() -> None:
    calls = {"download": 0}

    def downloader(**kwargs):
        calls["download"] += 1
        return b"compressed"

    result = _process(
        [
            _accepted(
                "A",
                etag="wrong-etag",
            )
        ],
        downloader=downloader,
        parser=lambda *args, **kwargs: [_chain("A")],
    )

    assert calls["download"] == 0
    assert result.audit_records == ()
    assert len(result.processing_errors) == 1
    assert result.chain_accounting_valid is True

    error = result.processing_errors[0]

    assert (
        error["processing_stage"]
        == "accepted_lineage_validation"
    )
    assert error["error_type"] == "AcceptedLineageError"


def test_missing_source_chain_is_a_processing_error() -> None:
    result = _process(
        [_accepted("A")],
        parser=lambda *args, **kwargs: [_chain("B")],
    )

    assert result.audit_records == ()
    assert len(result.processing_errors) == 1
    assert result.chain_accounting_valid is True

    error = result.processing_errors[0]

    assert error["processing_stage"] == "source_chain_lookup"
    assert error["error_type"] == "SourceChainNotFoundError"


def test_gold_residue_mismatch_is_a_reconstruction_error() -> None:
    accepted = _accepted("A")
    accepted["retained_residue_count"] = 2
    accepted["retained_label_seq_ids"] = [1, 2]

    result = _process(
        [accepted],
        parser=lambda *args, **kwargs: [_chain("A")],
    )

    assert result.audit_records == ()
    assert len(result.processing_errors) == 1
    assert result.chain_accounting_valid is True

    error = result.processing_errors[0]

    assert (
        error["processing_stage"]
        == "retained_chain_reconstruction"
    )
    assert error["error_type"] == "ValueError"


def test_real_geometric_violation_is_an_audit_row_not_error() -> None:
    result = _process(
        [_accepted("A")],
        parser=lambda *args, **kwargs: [
            _chain(
                "A",
                collinear=True,
            )
        ],
    )

    assert len(result.audit_records) == 1
    assert result.processing_errors == ()
    assert result.chain_accounting_valid is True

    audit = result.audit_records[0]

    assert audit["passed"] is False
    assert audit["violation_count"] >= 1
    assert "definition_3_4_undefined_h" in (
        audit["violation_types"]
    )


def test_parse_failure_fans_out_to_chain_level_errors() -> None:
    calls = {"parse": 0}

    def parser(compressed_bytes, *, pdb_id):
        calls["parse"] += 1
        raise MMCIFParseError("synthetic parse failure")

    result = _process(
        [
            _accepted("A"),
            _accepted("B"),
        ],
        parser=parser,
    )

    assert calls["parse"] == 1
    assert result.audit_records == ()
    assert len(result.processing_errors) == 2
    assert result.chain_accounting_valid is True
    assert result.source_downloaded is True
    assert result.source_parsed is False

    assert {
        error["label_chain_id"]
        for error in result.processing_errors
    } == {"A", "B"}

    assert all(
        error["processing_stage"] == "mmcif_parse"
        for error in result.processing_errors
    )


def test_records_materialize_with_explicit_step2_schemas() -> None:
    result = _process(
        [_accepted("A")],
        parser=lambda *args, **kwargs: [_chain("A")],
    )

    tables = geometric_validation_records_to_tables(
        result.audit_records,
        result.processing_errors,
    )

    assert tables.audit.num_rows == 1
    assert tables.processing_errors.num_rows == 0

    assert tables.audit.schema.equals(
        GEOMETRIC_VALIDATION_AUDIT_SCHEMA,
        check_metadata=True,
    )
    assert tables.processing_errors.schema.equals(
        GEOMETRIC_VALIDATION_PROCESSING_ERROR_SCHEMA,
        check_metadata=True,
    )


def test_geometric_validation_shards_are_written_atomically(
    tmp_path,
) -> None:
    result = _process(
        [_accepted("A")],
        parser=lambda *args, **kwargs: [_chain("A")],
    )

    tables = geometric_validation_records_to_tables(
        result.audit_records,
        result.processing_errors,
    )

    paths = write_geometric_validation_shards(
        tables,
        tmp_path,
        task_id=7,
    )

    assert paths == {
        "audit": tmp_path / "audit" / "task_7.parquet",
        "errors": tmp_path / "errors" / "task_7.parquet",
    }

    assert paths["audit"].is_file()
    assert paths["errors"].is_file()

    assert not paths["audit"].with_suffix(
        ".parquet.tmp"
    ).exists()
    assert not paths["errors"].with_suffix(
        ".parquet.tmp"
    ).exists()

    assert pq.read_schema(paths["audit"]).equals(
        GEOMETRIC_VALIDATION_AUDIT_SCHEMA,
        check_metadata=True,
    )

    assert pq.read_schema(paths["errors"]).equals(
        GEOMETRIC_VALIDATION_PROCESSING_ERROR_SCHEMA,
        check_metadata=True,
    )


def _manifest_for(pdb_id: str) -> dict:
    """Build one synthetic manifest row with matching source lineage."""

    normalized = pdb_id.lower()

    return {
        "snapshot": "20310415",
        "pdb_id": normalized,
        "s3_key": (
            f"20310415/coordinates/{normalized}.cif.gz"
        ),
        "size_bytes": 123,
        "etag": f"etag-{normalized}",
    }


def _accepted_for(
    pdb_id: str,
    chain_id: str,
) -> dict:
    """Build one accepted Gold row matching _manifest_for()."""

    normalized = pdb_id.lower()

    row = _accepted(chain_id)

    row["pdb_id"] = normalized
    row["source_mmcif_key"] = (
        f"20310415/coordinates/{normalized}.cif.gz"
    )
    row["source_etag"] = f"etag-{normalized}"

    return row


def test_batch_processes_only_relevant_manifest_sources() -> None:
    manifests = [
        _manifest_for("aaaa"),
        _manifest_for("bbbb"),
        _manifest_for("cccc"),
    ]

    # Deliberately reverse accepted-source order. The batch should still
    # aggregate in validated manifest order.
    accepted = [
        _accepted_for("cccc", "C"),
        _accepted_for("bbbb", "B"),
    ]

    calls = []

    def source_processor(
        manifest_row,
        accepted_rows,
        **kwargs,
    ):
        rows = list(accepted_rows)
        pdb_id = manifest_row["pdb_id"].lower()

        calls.append(pdb_id)

        return process_geometric_validation_source(
            manifest_row,
            rows,
            bucket_url=kwargs["bucket_url"],
            config=kwargs["config"],
            geometric_validation_pipeline_git_commit=(
                kwargs[
                    "geometric_validation_pipeline_git_commit"
                ]
            ),
            timeout_seconds=kwargs["timeout_seconds"],
            max_retries=kwargs["max_retries"],
            downloader=lambda **download_kwargs: b"compressed",
            parser=lambda *args, **parser_kwargs: [
                _chain(row["label_chain_id"])
                for row in rows
            ],
        )

    batch = process_geometric_validation_batch(
        manifests,
        accepted,
        bucket_url="https://example.invalid",
        config=GeometricValidationConfig(),
        geometric_validation_pipeline_git_commit="g" * 40,
        download_concurrency=1,
        source_processor=source_processor,
    )

    # "aaaa" has no accepted chain and therefore must never be processed.
    assert calls == ["bbbb", "cccc"]

    assert batch.manifest_source_object_count == 3
    assert batch.relevant_source_object_count == 2
    assert batch.input_accepted_chain_count == 2
    assert batch.downloaded_source_object_count == 2
    assert batch.parsed_source_object_count == 2

    assert batch.audit_chain_count == 2
    assert batch.processing_error_count == 0
    assert batch.chain_accounting_valid is True

    # Aggregation follows manifest order, not accepted-shard input order.
    assert [
        value
        for value in batch.tables.audit["pdb_id"].to_pylist()
    ] == ["bbbb", "cccc"]


def test_batch_allows_empty_accepted_shard_without_processing_sources() -> None:
    calls = []

    def source_processor(*args, **kwargs):
        calls.append("called")
        raise AssertionError(
            "No source should be processed for an empty accepted shard"
        )

    batch = process_geometric_validation_batch(
        [
            _manifest_for("aaaa"),
            _manifest_for("bbbb"),
        ],
        [],
        bucket_url="https://example.invalid",
        config=GeometricValidationConfig(),
        geometric_validation_pipeline_git_commit="g" * 40,
        source_processor=source_processor,
    )

    assert calls == []

    assert batch.manifest_source_object_count == 2
    assert batch.relevant_source_object_count == 0
    assert batch.input_accepted_chain_count == 0
    assert batch.downloaded_source_object_count == 0
    assert batch.parsed_source_object_count == 0

    assert batch.audit_chain_count == 0
    assert batch.processing_error_count == 0
    assert batch.chain_accounting_valid is True


def test_batch_rejects_accepted_source_outside_manifest_partition() -> None:
    with pytest.raises(
        GeometricValidationRunnerError,
        match="outside the matching manifest partition",
    ):
        process_geometric_validation_batch(
            [_manifest_for("aaaa")],
            [_accepted_for("bbbb", "B")],
            bucket_url="https://example.invalid",
            config=GeometricValidationConfig(),
            geometric_validation_pipeline_git_commit="g" * 40,
        )


def test_batch_rejects_duplicate_accepted_chain_identity() -> None:
    row = _accepted_for("aaaa", "A")

    with pytest.raises(
        GeometricValidationRunnerError,
        match="duplicate chain identity",
    ):
        process_geometric_validation_batch(
            [_manifest_for("aaaa")],
            [
                dict(row),
                dict(row),
            ],
            bucket_url="https://example.invalid",
            config=GeometricValidationConfig(),
            geometric_validation_pipeline_git_commit="g" * 40,
        )


def test_batch_rejects_invalid_source_level_chain_accounting() -> None:
    def source_processor(
        manifest_row,
        accepted_rows,
        **kwargs,
    ):
        return SourceGeometricValidationResult(
            pdb_id=manifest_row["pdb_id"],
            input_accepted_chain_count=1,
            audit_records=(),
            processing_errors=(),
            source_downloaded=True,
            source_parsed=True,
        )

    with pytest.raises(
        GeometricValidationRunnerError,
        match="Source-level Step-2 chain accounting failed",
    ):
        process_geometric_validation_batch(
            [_manifest_for("aaaa")],
            [_accepted_for("aaaa", "A")],
            bucket_url="https://example.invalid",
            config=GeometricValidationConfig(),
            geometric_validation_pipeline_git_commit="g" * 40,
            source_processor=source_processor,
        )


def _mixed_terminal_batch() -> GeometricValidationBatchResult:
    """One pass, one real geometry violation, and one processing error."""

    result = _process(
        [
            _accepted("A"),
            _accepted("B"),
            _accepted("C"),
        ],
        parser=lambda *args, **kwargs: [
            _chain("A"),
            _chain("B", collinear=True),
            # C deliberately absent -> source-chain lookup error.
        ],
    )

    assert result.chain_accounting_valid is True

    tables = geometric_validation_records_to_tables(
        result.audit_records,
        result.processing_errors,
    )

    return GeometricValidationBatchResult(
        manifest_source_object_count=1,
        relevant_source_object_count=1,
        input_accepted_chain_count=3,
        downloaded_source_object_count=1,
        parsed_source_object_count=1,
        tables=tables,
    )


def _task_context(
    *,
    task_id: str = "7",
) -> GeometricValidationTaskContext:
    return GeometricValidationTaskContext(
        task_id=task_id,
        snapshot="20310415",
        cleaning_protocol="protocol3.2-comp702-v1",
        quality_pipeline_git_commit="q" * 40,
        geometric_validation_pipeline_git_commit="g" * 40,
        configured_minimum_backbone_distance_angstrom=0.010,
        configured_minimum_triangle_angle_degrees=3.0,
        started_at_utc="2031-04-15T12:00:00Z",
        completed_at_utc="2031-04-15T12:00:05Z",
        runtime_seconds=5.0,
        slurm_job_id="12345",
        slurm_array_task_id="7",
        peak_memory_bytes=123456,
    )


def test_task_summary_accounts_all_terminal_chain_outcomes() -> None:
    batch = _mixed_terminal_batch()

    summary = build_geometric_validation_task_summary(
        batch,
        _task_context(),
    )

    assert summary["manifest_source_object_count"] == 1
    assert summary["relevant_source_object_count"] == 1
    assert summary["downloaded_source_object_count"] == 1
    assert summary["parsed_source_object_count"] == 1

    assert summary["input_accepted_chain_count"] == 3
    assert summary["audit_chain_count"] == 2
    assert summary["geometric_passed_chain_count"] == 1
    assert summary["geometric_violated_chain_count"] == 1
    assert summary["processing_error_count"] == 1

    assert summary["chain_accounting_valid"] is True

    assert (
        summary["processing_errors_by_stage"]
        ["source_chain_lookup"]
        == 1
    )
    assert (
        summary["processing_errors_by_type"]
        ["SourceChainNotFoundError"]
        == 1
    )

    assert (
        summary["violations_by_type"]
        ["definition_3_4_undefined_h"]
        == 1
    )


def test_task_summary_is_written_atomically(
    tmp_path,
) -> None:
    batch = _mixed_terminal_batch()

    summary = build_geometric_validation_task_summary(
        batch,
        _task_context(),
    )

    path = write_geometric_validation_task_summary_atomic(
        summary,
        tmp_path,
    )

    assert path == (
        tmp_path
        / "summaries"
        / "task_7.json"
    )
    assert path.is_file()

    temporary = path.with_suffix(
        path.suffix + ".tmp"
    )
    assert not temporary.exists()

    written = json.loads(
        path.read_text(encoding="utf-8")
    )

    assert written == summary


def test_publish_writes_parquet_before_summary(
    tmp_path,
) -> None:
    batch = _mixed_terminal_batch()
    context = _task_context()

    calls = []

    expected_shards = {
        "audit": tmp_path / "audit" / "task_7.parquet",
        "errors": tmp_path / "errors" / "task_7.parquet",
    }
    expected_summary = (
        tmp_path
        / "summaries"
        / "task_7.json"
    )

    def shard_writer(
        tables,
        output_root,
        task_id,
    ):
        calls.append("shards")

        assert tables is batch.tables
        assert output_root == tmp_path
        assert task_id == "7"

        return expected_shards

    def summary_writer(
        summary,
        output_root,
    ):
        calls.append("summary")

        assert output_root == tmp_path
        assert summary["chain_accounting_valid"] is True

        return expected_summary

    def utc_now():
        calls.append("utc")
        return "2031-04-15T12:00:05Z"

    def perf_counter():
        calls.append("perf")
        return 105.0

    def peak_memory_reader():
        calls.append("peak")
        return 654321

    publication = publish_geometric_validation_batch(
        batch,
        output_root=tmp_path,
        context=context,
        started_perf_counter=100.0,
        shard_writer=shard_writer,
        summary_writer=summary_writer,
        utc_now=utc_now,
        perf_counter=perf_counter,
        peak_memory_reader=peak_memory_reader,
    )

    assert calls == [
        "shards",
        "utc",
        "perf",
        "peak",
        "summary",
    ]

    assert publication.shard_paths == expected_shards
    assert publication.summary_path == expected_summary
    assert publication.summary["chain_accounting_valid"] is True
    assert publication.summary["completed_at_utc"] == (
        "2031-04-15T12:00:05Z"
    )
    assert publication.summary["runtime_seconds"] == 5.0
    assert publication.summary["peak_memory_bytes"] == 654321


def test_publish_rejects_bad_accounting_before_any_write(
    tmp_path,
) -> None:
    valid = _mixed_terminal_batch()

    invalid = GeometricValidationBatchResult(
        manifest_source_object_count=(
            valid.manifest_source_object_count
        ),
        relevant_source_object_count=(
            valid.relevant_source_object_count
        ),
        # There are only three terminal rows in the tables.
        input_accepted_chain_count=4,
        downloaded_source_object_count=(
            valid.downloaded_source_object_count
        ),
        parsed_source_object_count=(
            valid.parsed_source_object_count
        ),
        tables=valid.tables,
    )

    calls = []

    def shard_writer(*args, **kwargs):
        calls.append("shards")
        raise AssertionError("Shard writer must not be called")

    def summary_writer(*args, **kwargs):
        calls.append("summary")
        raise AssertionError("Summary writer must not be called")

    with pytest.raises(
        GeometricValidationRunnerError,
        match="chain accounting failed",
    ):
        publish_geometric_validation_batch(
            invalid,
            output_root=tmp_path,
            context=_task_context(),
            started_perf_counter=100.0,
            shard_writer=shard_writer,
            summary_writer=summary_writer,
        )

    assert calls == []


def _quality_summary(
    *,
    accepted_chain_count: int = 1,
    input_source_object_count: int = 1,
) -> dict:
    return {
        "summary_schema_name": "pdbclean_quality_task_summary",
        "summary_schema_version": "1.0",
        "task_id": "7",
        "snapshot": "20310415",
        "cleaning_protocol": "protocol3.2-comp702-v1",
        "pipeline_git_commit": "q" * 40,
        "input_source_object_count": input_source_object_count,
        "accepted_chain_count": accepted_chain_count,
        "source_object_accounting_valid": True,
        "selected_chain_accounting_valid": True,
    }


def test_upstream_quality_task_validation_accepts_completed_task() -> None:
    row = _accepted("A")

    validated = validate_upstream_quality_task(
        _quality_summary(),
        [row],
        expected_task_id=7,
        expected_snapshot="20310415",
        expected_cleaning_protocol="protocol3.2-comp702-v1",
        expected_manifest_source_object_count=1,
    )

    assert validated.accepted_rows == (row,)
    assert validated.quality_pipeline_git_commit == "q" * 40


def test_upstream_quality_task_validation_supports_empty_accepted_shard() -> None:
    validated = validate_upstream_quality_task(
        _quality_summary(
            accepted_chain_count=0,
            input_source_object_count=500,
        ),
        [],
        expected_task_id=7,
        expected_snapshot="20310415",
        expected_cleaning_protocol="protocol3.2-comp702-v1",
        expected_manifest_source_object_count=500,
    )

    assert validated.accepted_rows == ()
    assert validated.quality_pipeline_git_commit == "q" * 40


def test_upstream_quality_task_validation_rejects_bad_completion_marker() -> None:
    summary = _quality_summary()
    summary["selected_chain_accounting_valid"] = False

    with pytest.raises(
        GeometricValidationRunnerError,
        match="selected-chain accounting",
    ):
        validate_upstream_quality_task(
            summary,
            [_accepted("A")],
            expected_task_id=7,
            expected_snapshot="20310415",
            expected_cleaning_protocol="protocol3.2-comp702-v1",
            expected_manifest_source_object_count=1,
        )


def test_upstream_quality_task_validation_rejects_count_mismatch() -> None:
    with pytest.raises(
        GeometricValidationRunnerError,
        match="accepted-chain count",
    ):
        validate_upstream_quality_task(
            _quality_summary(
                accepted_chain_count=2,
            ),
            [_accepted("A")],
            expected_task_id=7,
            expected_snapshot="20310415",
            expected_cleaning_protocol="protocol3.2-comp702-v1",
            expected_manifest_source_object_count=1,
        )


def test_upstream_quality_task_validation_rejects_producer_mismatch() -> None:
    row = _accepted("A")
    row["pipeline_git_commit"] = "x" * 40

    with pytest.raises(
        GeometricValidationRunnerError,
        match="producer Git commit",
    ):
        validate_upstream_quality_task(
            _quality_summary(),
            [row],
            expected_task_id=7,
            expected_snapshot="20310415",
            expected_cleaning_protocol="protocol3.2-comp702-v1",
            expected_manifest_source_object_count=1,
        )


def test_execute_geometric_validation_task_orchestrates_batch_then_publication(
    tmp_path,
) -> None:
    calls = []
    batch = _mixed_terminal_batch()

    clock_values = iter(
        [
            100.0,
            105.0,
        ]
    )

    def utc_now():
        calls.append("utc")
        return "2031-04-15T12:00:00Z"

    def perf_counter():
        calls.append("perf")
        return next(clock_values)

    def batch_processor(
        manifest_rows,
        accepted_rows,
        **kwargs,
    ):
        calls.append("batch")

        assert list(manifest_rows) == [
            {"pdb_id": "aaaa"},
        ]
        assert list(accepted_rows) == [
            {"pdb_id": "aaaa", "model_id": 1},
        ]

        assert kwargs == {
            "bucket_url": "https://example.invalid",
            "config": GeometricValidationConfig(
                minimum_backbone_distance_angstrom=0.02,
                minimum_triangle_angle_degrees=5.0,
            ),
            "geometric_validation_pipeline_git_commit": "g" * 40,
            "timeout_seconds": 37,
            "max_retries": 4,
            "download_concurrency": 3,
        }

        return batch

    expected = object()

    def publisher(
        received_batch,
        **kwargs,
    ):
        calls.append("publisher")

        assert received_batch is batch
        assert kwargs["output_root"] == tmp_path
        assert kwargs["started_perf_counter"] == 100.0
        assert kwargs["utc_now"] is utc_now
        assert kwargs["perf_counter"] is perf_counter

        context = kwargs["context"]

        assert context.task_id == "7"
        assert context.snapshot == "20310415"
        assert context.cleaning_protocol == (
            "protocol3.2-comp702-v1"
        )
        assert context.quality_pipeline_git_commit == "q" * 40
        assert (
            context.geometric_validation_pipeline_git_commit
            == "g" * 40
        )

        assert (
            context.configured_minimum_backbone_distance_angstrom
            == 0.02
        )
        assert (
            context.configured_minimum_triangle_angle_degrees
            == 5.0
        )

        assert context.started_at_utc == (
            "2031-04-15T12:00:00Z"
        )
        assert context.completed_at_utc is None
        assert context.runtime_seconds is None

        assert context.slurm_job_id == "12345"
        assert context.slurm_array_task_id == "7"

        return expected

    result = execute_geometric_validation_task(
        [{"pdb_id": "aaaa"}],
        [{"pdb_id": "aaaa", "model_id": 1}],
        output_root=tmp_path,
        task_id=7,
        snapshot="20310415",
        bucket_url="https://example.invalid",
        config=GeometricValidationConfig(
            minimum_backbone_distance_angstrom=0.02,
            minimum_triangle_angle_degrees=5.0,
        ),
        cleaning_protocol="protocol3.2-comp702-v1",
        quality_pipeline_git_commit="q" * 40,
        geometric_validation_pipeline_git_commit="g" * 40,
        timeout_seconds=37,
        max_retries=4,
        download_concurrency=3,
        environ={
            "SLURM_JOB_ID": "12345",
            "SLURM_ARRAY_TASK_ID": "7",
        },
        batch_processor=batch_processor,
        publisher=publisher,
        utc_now=utc_now,
        perf_counter=perf_counter,
    )

    assert result is expected

    # Start timestamp and monotonic clock are captured before processing.
    assert calls == [
        "utc",
        "perf",
        "batch",
        "publisher",
    ]


def test_execute_geometric_validation_task_rejects_unsafe_task_id(
    tmp_path,
) -> None:
    with pytest.raises(
        GeometricValidationRunnerError,
        match="Unsafe geometric-validation task_id",
    ):
        execute_geometric_validation_task(
            [],
            [],
            output_root=tmp_path,
            task_id="../7",
            snapshot="20310415",
            bucket_url="https://example.invalid",
            config=GeometricValidationConfig(),
            cleaning_protocol="protocol3.2-comp702-v1",
            quality_pipeline_git_commit="q" * 40,
            geometric_validation_pipeline_git_commit="g" * 40,
        )
