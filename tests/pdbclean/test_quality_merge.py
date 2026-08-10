"""Tests for distributed quality-stage merge and validation."""

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from pdbclean.gold import (
    QUALITY_TASK_SUMMARY_SCHEMA_NAME,
    QUALITY_TASK_SUMMARY_SCHEMA_VERSION,
)
from pdbclean.quality_merge import (
    QUALITY_SHARD_SCHEMAS,
    QualityMergeError,
    discover_quality_task_artifacts,
    expected_quality_task_ids,
)


SNAPSHOT = "20260101"
PROTOCOL = "protocol3.2-comp702-v1"
GIT_COMMIT = "a" * 40


def _summary(task_id: int) -> dict:
    return {
        "summary_schema_name": QUALITY_TASK_SUMMARY_SCHEMA_NAME,
        "summary_schema_version": QUALITY_TASK_SUMMARY_SCHEMA_VERSION,
        "task_id": str(task_id),
        "snapshot": SNAPSHOT,
        "cleaning_protocol": PROTOCOL,
        "pipeline_git_commit": GIT_COMMIT,
        "source_object_accounting_valid": True,
        "selected_chain_accounting_valid": True,
    }


def _write_empty_shard(
    path: Path,
    schema: pa.Schema,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pylist([], schema=schema)
    pq.write_table(table, path)


def _write_task(
    root: Path,
    task_id: int,
) -> None:
    summary_path = (
        root
        / "summaries"
        / f"task_{task_id}.json"
    )
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps(_summary(task_id)) + "\n",
        encoding="utf-8",
    )

    for shard_name, schema in QUALITY_SHARD_SCHEMAS.items():
        _write_empty_shard(
            root
            / shard_name
            / f"task_{task_id}.parquet",
            schema,
        )


def test_expected_quality_task_ids_are_dynamic() -> None:
    assert expected_quality_task_ids(1, 500) == (0,)
    assert expected_quality_task_ids(500, 500) == (0,)
    assert expected_quality_task_ids(501, 500) == (0, 1)
    assert expected_quality_task_ids(246905, 500) == tuple(
        range(494)
    )


def test_discover_complete_quality_task_set(
    tmp_path: Path,
) -> None:
    _write_task(tmp_path, 0)
    _write_task(tmp_path, 1)

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0, 1),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    assert tuple(
        artifact.task_id
        for artifact in artifacts
    ) == (0, 1)

    for artifact in artifacts:
        assert artifact.summary_path.is_file()
        assert set(artifact.shard_paths) == set(
            QUALITY_SHARD_SCHEMAS
        )


def test_missing_summary_is_rejected(
    tmp_path: Path,
) -> None:
    _write_task(tmp_path, 0)
    _write_task(tmp_path, 1)

    (
        tmp_path
        / "summaries"
        / "task_1.json"
    ).unlink()

    with pytest.raises(
        QualityMergeError,
        match="missing=\\[1\\]",
    ):
        discover_quality_task_artifacts(
            tmp_path,
            expected_task_ids=(0, 1),
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
            expected_pipeline_git_commit=GIT_COMMIT,
        )


def test_missing_shard_is_rejected(
    tmp_path: Path,
) -> None:
    _write_task(tmp_path, 0)

    (
        tmp_path
        / "accepted"
        / "task_0.parquet"
    ).unlink()

    with pytest.raises(
        QualityMergeError,
        match="accepted shard.*missing=\\[0\\]",
    ):
        discover_quality_task_artifacts(
            tmp_path,
            expected_task_ids=(0,),
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
            expected_pipeline_git_commit=GIT_COMMIT,
        )


def test_unexpected_task_is_rejected(
    tmp_path: Path,
) -> None:
    _write_task(tmp_path, 0)
    _write_task(tmp_path, 1)

    with pytest.raises(
        QualityMergeError,
        match="unexpected=\\[1\\]",
    ):
        discover_quality_task_artifacts(
            tmp_path,
            expected_task_ids=(0,),
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
            expected_pipeline_git_commit=GIT_COMMIT,
        )


def test_duplicate_numeric_task_id_is_rejected(
    tmp_path: Path,
) -> None:
    _write_task(tmp_path, 1)

    duplicate = (
        tmp_path
        / "summaries"
        / "task_01.json"
    )
    duplicate.write_text(
        json.dumps(_summary(1)) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(
        QualityMergeError,
        match="Duplicate summary task ID: 1",
    ):
        discover_quality_task_artifacts(
            tmp_path,
            expected_task_ids=(0, 1),
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
            expected_pipeline_git_commit=GIT_COMMIT,
        )


def test_temporary_file_is_rejected(
    tmp_path: Path,
) -> None:
    _write_task(tmp_path, 0)

    temporary = (
        tmp_path
        / "accepted"
        / "task_0.parquet.tmp"
    )
    temporary.write_text("partial", encoding="utf-8")

    with pytest.raises(
        QualityMergeError,
        match="Temporary quality-stage files remain",
    ):
        discover_quality_task_artifacts(
            tmp_path,
            expected_task_ids=(0,),
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
            expected_pipeline_git_commit=GIT_COMMIT,
        )


def test_summary_provenance_mismatch_is_rejected(
    tmp_path: Path,
) -> None:
    _write_task(tmp_path, 0)

    summary_path = (
        tmp_path
        / "summaries"
        / "task_0.json"
    )
    summary = _summary(0)
    summary["snapshot"] = "wrong-snapshot"
    summary_path.write_text(
        json.dumps(summary) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(
        QualityMergeError,
        match="Snapshot mismatch",
    ):
        discover_quality_task_artifacts(
            tmp_path,
            expected_task_ids=(0,),
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
            expected_pipeline_git_commit=GIT_COMMIT,
        )


def test_shard_schema_mismatch_is_rejected(
    tmp_path: Path,
) -> None:
    _write_task(tmp_path, 0)

    bad_schema = pa.schema(
        [pa.field("wrong", pa.string())],
        metadata={
            b"schema_name": b"wrong",
            b"schema_version": b"1.0",
        },
    )

    pq.write_table(
        pa.Table.from_pylist([], schema=bad_schema),
        tmp_path / "accepted" / "task_0.parquet",
    )

    with pytest.raises(
        QualityMergeError,
        match="Unexpected Parquet schema",
    ):
        discover_quality_task_artifacts(
            tmp_path,
            expected_task_ids=(0,),
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
            expected_pipeline_git_commit=GIT_COMMIT,
        )


def _write_accounting_summary(
    root: Path,
    task_id: int,
    *,
    input_sources: int,
    successful_sources: int | None = None,
    failed_sources: int = 0,
    selected_chains: int = 0,
    accepted: int = 0,
    accepted_trimmed: int = 0,
    rejected: int = 0,
    non_candidates: int = 0,
    dirty: int = 0,
    errors: int = 0,
    chain_errors: int = 0,
    source_errors: int = 0,
) -> None:
    if successful_sources is None:
        successful_sources = input_sources - failed_sources

    summary = _summary(task_id)
    summary.update(
        {
            "input_source_object_count": input_sources,
            "successful_source_object_count": successful_sources,
            "failed_source_object_count": failed_sources,
            "parsed_silver_chain_count": selected_chains,
            "selected_silver_chain_count": selected_chains,
            "candidate_entry_count": 0,
            "candidate_chain_count": 0,
            "non_candidate_chain_count": non_candidates,
            "accepted_chain_count": accepted,
            "accepted_trimmed_chain_count": accepted_trimmed,
            "rejected_chain_count": rejected,
            "dirty_residue_count": dirty,
            "processing_error_count": errors,
            "chain_level_processing_error_count": chain_errors,
            "source_entry_processing_error_count": source_errors,
            "total_gold_record_count": (
                accepted
                + rejected
                + non_candidates
                + dirty
                + errors
            ),
        }
    )

    path = root / "summaries" / f"task_{task_id}.json"
    path.write_text(
        json.dumps(summary) + "\n",
        encoding="utf-8",
    )


def test_validate_quality_task_accounting_dynamic_partition_sizes(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import (
        validate_quality_task_accounting,
    )

    _write_task(tmp_path, 0)
    _write_task(tmp_path, 1)

    _write_accounting_summary(
        tmp_path,
        0,
        input_sources=500,
    )
    _write_accounting_summary(
        tmp_path,
        1,
        input_sources=1,
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0, 1),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    validated = validate_quality_task_accounting(
        artifacts,
        manifest_row_count=501,
        batch_size=500,
    )

    assert tuple(
        item.expected_input_source_object_count
        for item in validated
    ) == (500, 1)


def test_task_input_source_count_must_match_partition(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import (
        validate_quality_task_accounting,
    )

    _write_task(tmp_path, 0)
    _write_accounting_summary(
        tmp_path,
        0,
        input_sources=499,
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0,),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    with pytest.raises(
        QualityMergeError,
        match="input source count mismatch",
    ):
        validate_quality_task_accounting(
            artifacts,
            manifest_row_count=500,
            batch_size=500,
        )


def test_task_shard_row_count_must_match_summary(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import (
        validate_quality_task_accounting,
    )

    _write_task(tmp_path, 0)
    _write_accounting_summary(
        tmp_path,
        0,
        input_sources=1,
        selected_chains=1,
        accepted=1,
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0,),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    with pytest.raises(
        QualityMergeError,
        match="accepted row-count mismatch",
    ):
        validate_quality_task_accounting(
            artifacts,
            manifest_row_count=1,
            batch_size=500,
        )


def test_task_source_accounting_is_recomputed(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import (
        validate_quality_task_accounting,
    )

    _write_task(tmp_path, 0)
    _write_accounting_summary(
        tmp_path,
        0,
        input_sources=1,
        successful_sources=1,
        failed_sources=1,
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0,),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    with pytest.raises(
        QualityMergeError,
        match="source-object accounting mismatch",
    ):
        validate_quality_task_accounting(
            artifacts,
            manifest_row_count=1,
            batch_size=500,
        )


def test_task_selected_chain_accounting_is_recomputed(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import (
        validate_quality_task_accounting,
    )

    _write_task(tmp_path, 0)
    _write_accounting_summary(
        tmp_path,
        0,
        input_sources=1,
        selected_chains=1,
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0,),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    with pytest.raises(
        QualityMergeError,
        match="selected-chain accounting mismatch",
    ):
        validate_quality_task_accounting(
            artifacts,
            manifest_row_count=1,
            batch_size=500,
        )


def test_total_gold_count_is_recomputed(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import (
        validate_quality_task_accounting,
    )

    _write_task(tmp_path, 0)
    _write_accounting_summary(
        tmp_path,
        0,
        input_sources=1,
    )

    summary_path = (
        tmp_path
        / "summaries"
        / "task_0.json"
    )
    summary = json.loads(
        summary_path.read_text(encoding="utf-8")
    )
    summary["total_gold_record_count"] = 1
    summary_path.write_text(
        json.dumps(summary) + "\n",
        encoding="utf-8",
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0,),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    with pytest.raises(
        QualityMergeError,
        match="total Gold record count mismatch",
    ):
        validate_quality_task_accounting(
            artifacts,
            manifest_row_count=1,
            batch_size=500,
        )


def _manifest_table(
    rows: list[tuple[str, str, str]],
) -> pa.Table:
    return pa.table(
        {
            "snapshot": [SNAPSHOT for _ in rows],
            "source_layout": [
                "canonical_divided_mmcif"
                for _ in rows
            ],
            "pdb_id": [row[0] for row in rows],
            "s3_key": [row[1] for row in rows],
            "size_bytes": [123 for _ in rows],
            "etag": [row[2] for row in rows],
        }
    )


def _accepted_row(
    *,
    pdb_id: str,
    chain_id: str = "A",
    source_key: str,
    etag: str,
) -> dict:
    return {
        "snapshot": SNAPSHOT,
        "pdb_id": pdb_id,
        "model_id": 1,
        "label_chain_id": chain_id,
        "auth_chain_id": chain_id,
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
        "source_mmcif_key": source_key,
        "source_etag": etag,
        "cleaning_protocol": PROTOCOL,
        "pipeline_git_commit": GIT_COMMIT,
    }


def _non_candidate_row(
    *,
    pdb_id: str,
    chain_id: str = "A",
    source_key: str,
    etag: str,
) -> dict:
    return {
        "snapshot": SNAPSHOT,
        "pdb_id": pdb_id,
        "model_id": 1,
        "label_chain_id": chain_id,
        "auth_chain_id": chain_id,
        "entity_id": "1",
        "terminal_status": "non_candidate",
        "terminal_reason": "no_protocol32_backbone_atoms",
        "terminal_stage": "candidate_selection",
        "source_mmcif_key": source_key,
        "source_etag": etag,
        "cleaning_protocol": PROTOCOL,
        "pipeline_git_commit": GIT_COMMIT,
    }


def _processing_error_row(
    *,
    pdb_id: str,
    chain_id: str | None,
    source_key: str,
    etag: str,
) -> dict:
    return {
        "snapshot": SNAPSHOT,
        "pdb_id": pdb_id,
        "model_id": 1 if chain_id is not None else None,
        "label_chain_id": chain_id,
        "processing_stage": "quality_cleaning",
        "error_type": "RuntimeError",
        "error_message": "synthetic failure",
        "source_mmcif_key": source_key,
        "source_etag": etag,
        "pipeline_git_commit": GIT_COMMIT,
    }


def _dirty_row(
    *,
    pdb_id: str,
    chain_id: str = "A",
    residue_id: int = 1,
    source_key: str,
    etag: str,
) -> dict:
    return {
        "snapshot": SNAPSHOT,
        "pdb_id": pdb_id,
        "model_id": 1,
        "label_chain_id": chain_id,
        "auth_chain_id": chain_id,
        "entity_id": "1",
        "label_seq_id": residue_id,
        "deposited_residue_name": "ALA",
        "mapped_residue_code": "A",
        "rule_id": "Q002",
        "dirty_type": "disorder",
        "cleaning_stage": "Q002",
        "details_json": "{}",
        "source_mmcif_key": source_key,
        "source_etag": etag,
    }


def _overwrite_shard(
    root: Path,
    task_id: int,
    shard_name: str,
    rows: list[dict],
) -> None:
    schema = QUALITY_SHARD_SCHEMAS[shard_name]
    path = root / shard_name / f"task_{task_id}.parquet"
    pq.write_table(
        pa.Table.from_pylist(rows, schema=schema),
        path,
    )


def test_validate_quality_global_state_accepts_valid_provenance(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import validate_quality_global_state

    key = "20260101/pub/pdb/data/structures/divided/mmCIF/aa/1aaa.cif.gz"
    etag = '"etag-1"'

    _write_task(tmp_path, 0)
    _overwrite_shard(
        tmp_path,
        0,
        "accepted",
        [
            _accepted_row(
                pdb_id="1aaa",
                source_key=key,
                etag=etag,
            )
        ],
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0,),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    result = validate_quality_global_state(
        artifacts,
        manifest=_manifest_table(
            [("1aaa", key, etag)]
        ),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    assert result.accepted_chain_count == 1
    assert result.unique_outcome_chain_count == 1


def test_conflicting_gold_chain_identity_is_rejected(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import validate_quality_global_state

    key = "20260101/pub/pdb/data/structures/divided/mmCIF/aa/1aaa.cif.gz"
    etag = '"etag-1"'

    _write_task(tmp_path, 0)
    _write_task(tmp_path, 1)

    _overwrite_shard(
        tmp_path,
        0,
        "accepted",
        [
            _accepted_row(
                pdb_id="1aaa",
                source_key=key,
                etag=etag,
            )
        ],
    )
    _overwrite_shard(
        tmp_path,
        1,
        "non_candidates",
        [
            _non_candidate_row(
                pdb_id="1aaa",
                source_key=key,
                etag=etag,
            )
        ],
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0, 1),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    with pytest.raises(
        QualityMergeError,
        match="Duplicate or conflicting Gold chain identity",
    ):
        validate_quality_global_state(
            artifacts,
            manifest=_manifest_table(
                [("1aaa", key, etag)]
            ),
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
            expected_pipeline_git_commit=GIT_COMMIT,
        )


def test_gold_source_etag_must_match_manifest(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import validate_quality_global_state

    key = "20260101/pub/pdb/data/structures/divided/mmCIF/aa/1aaa.cif.gz"

    _write_task(tmp_path, 0)
    _overwrite_shard(
        tmp_path,
        0,
        "accepted",
        [
            _accepted_row(
                pdb_id="1aaa",
                source_key=key,
                etag='"wrong-etag"',
            )
        ],
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0,),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    with pytest.raises(
        QualityMergeError,
        match="source provenance mismatch",
    ):
        validate_quality_global_state(
            artifacts,
            manifest=_manifest_table(
                [("1aaa", key, '"correct-etag"')]
            ),
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
            expected_pipeline_git_commit=GIT_COMMIT,
        )


def test_chain_cannot_have_outcome_and_processing_error(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import validate_quality_global_state

    key = "20260101/pub/pdb/data/structures/divided/mmCIF/aa/1aaa.cif.gz"
    etag = '"etag-1"'

    _write_task(tmp_path, 0)

    _overwrite_shard(
        tmp_path,
        0,
        "accepted",
        [
            _accepted_row(
                pdb_id="1aaa",
                source_key=key,
                etag=etag,
            )
        ],
    )
    _overwrite_shard(
        tmp_path,
        0,
        "errors",
        [
            _processing_error_row(
                pdb_id="1aaa",
                chain_id="A",
                source_key=key,
                etag=etag,
            )
        ],
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0,),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    with pytest.raises(
        QualityMergeError,
        match="both a Gold outcome and a processing error",
    ):
        validate_quality_global_state(
            artifacts,
            manifest=_manifest_table(
                [("1aaa", key, etag)]
            ),
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
            expected_pipeline_git_commit=GIT_COMMIT,
        )


def test_duplicate_dirty_residue_identity_is_rejected(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import validate_quality_global_state

    key = "20260101/pub/pdb/data/structures/divided/mmCIF/aa/1aaa.cif.gz"
    etag = '"etag-1"'

    _write_task(tmp_path, 0)
    _write_task(tmp_path, 1)

    dirty = _dirty_row(
        pdb_id="1aaa",
        residue_id=7,
        source_key=key,
        etag=etag,
    )

    _overwrite_shard(
        tmp_path,
        0,
        "dirty_residues",
        [dirty],
    )
    _overwrite_shard(
        tmp_path,
        1,
        "dirty_residues",
        [dirty],
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0, 1),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    with pytest.raises(
        QualityMergeError,
        match="Duplicate dirty-residue identity",
    ):
        validate_quality_global_state(
            artifacts,
            manifest=_manifest_table(
                [("1aaa", key, etag)]
            ),
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
            expected_pipeline_git_commit=GIT_COMMIT,
        )


def test_global_validation_rejects_wrong_snapshot_manifest(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import validate_quality_global_state

    key = "wrong/pub/pdb/data/structures/divided/mmCIF/aa/1aaa.cif.gz"
    etag = '"etag-1"'

    _write_task(tmp_path, 0)

    manifest = pa.table(
        {
            "snapshot": ["wrong"],
            "source_layout": ["canonical_divided_mmcif"],
            "pdb_id": ["1aaa"],
            "s3_key": [key],
            "size_bytes": [123],
            "etag": [etag],
        }
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0,),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    with pytest.raises(
        QualityMergeError,
        match="Bronze manifest validation failed",
    ):
        validate_quality_global_state(
            artifacts,
            manifest=manifest,
            expected_snapshot=SNAPSHOT,
            expected_cleaning_protocol=PROTOCOL,
            expected_pipeline_git_commit=GIT_COMMIT,
        )


def test_publish_quality_merge_writes_merged_outputs_in_task_order(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import (
        QualityGlobalValidation,
        publish_quality_merge,
    )

    key_a = (
        "20260101/pub/pdb/data/structures/divided/"
        "mmCIF/aa/1aaa.cif.gz"
    )
    key_b = (
        "20260101/pub/pdb/data/structures/divided/"
        "mmCIF/aa/1aab.cif.gz"
    )

    _write_task(tmp_path, 0)
    _write_task(tmp_path, 1)

    _overwrite_shard(
        tmp_path,
        0,
        "accepted",
        [
            _accepted_row(
                pdb_id="1aaa",
                source_key=key_a,
                etag='"etag-a"',
            )
        ],
    )
    _overwrite_shard(
        tmp_path,
        1,
        "accepted",
        [
            _accepted_row(
                pdb_id="1aab",
                source_key=key_b,
                etag='"etag-b"',
            )
        ],
    )

    _write_accounting_summary(
        tmp_path,
        0,
        input_sources=1,
        selected_chains=1,
        accepted=1,
    )
    _write_accounting_summary(
        tmp_path,
        1,
        input_sources=1,
        selected_chains=1,
        accepted=1,
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0, 1),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    publication = publish_quality_merge(
        artifacts,
        quality_root=tmp_path,
        manifest_row_count=2,
        batch_size=1,
        snapshot=SNAPSHOT,
        cleaning_protocol=PROTOCOL,
        pipeline_git_commit=GIT_COMMIT,
        global_validation=QualityGlobalValidation(
            accepted_chain_count=2,
            rejected_chain_count=0,
            non_candidate_chain_count=0,
            dirty_residue_count=0,
            processing_error_count=0,
            unique_outcome_chain_count=2,
        ),
    )

    accepted = pq.read_table(
        publication.merged_paths["accepted"]
    )

    assert accepted["pdb_id"].to_pylist() == [
        "1aaa",
        "1aab",
    ]

    for shard_name, schema in QUALITY_SHARD_SCHEMAS.items():
        merged_path = publication.merged_paths[shard_name]
        observed = pq.read_schema(merged_path)

        assert observed.metadata == schema.metadata
        assert observed.equals(
            schema,
            check_metadata=False,
        )

    assert publication.global_summary_path.is_file()
    assert publication.success_path.is_file()

    assert publication.global_summary[
        "input_source_object_count"
    ] == 2
    assert publication.global_summary[
        "accepted_chain_count"
    ] == 2
    assert publication.global_summary[
        "source_object_accounting_valid"
    ] is True
    assert publication.global_summary[
        "selected_chain_accounting_valid"
    ] is True

    assert not list(tmp_path.rglob("*.tmp"))


def test_publish_quality_merge_rerun_does_not_duplicate_rows(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import (
        QualityGlobalValidation,
        publish_quality_merge,
    )

    key = (
        "20260101/pub/pdb/data/structures/divided/"
        "mmCIF/aa/1aaa.cif.gz"
    )

    _write_task(tmp_path, 0)
    _overwrite_shard(
        tmp_path,
        0,
        "accepted",
        [
            _accepted_row(
                pdb_id="1aaa",
                source_key=key,
                etag='"etag-a"',
            )
        ],
    )
    _write_accounting_summary(
        tmp_path,
        0,
        input_sources=1,
        selected_chains=1,
        accepted=1,
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0,),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    validation = QualityGlobalValidation(
        accepted_chain_count=1,
        rejected_chain_count=0,
        non_candidate_chain_count=0,
        dirty_residue_count=0,
        processing_error_count=0,
        unique_outcome_chain_count=1,
    )

    for _ in range(2):
        publish_quality_merge(
            artifacts,
            quality_root=tmp_path,
            manifest_row_count=1,
            batch_size=500,
            snapshot=SNAPSHOT,
            cleaning_protocol=PROTOCOL,
            pipeline_git_commit=GIT_COMMIT,
            global_validation=validation,
        )

    accepted = pq.read_table(
        tmp_path / "merged" / "accepted.parquet"
    )

    assert accepted.num_rows == 1
    assert accepted["pdb_id"].to_pylist() == ["1aaa"]
    assert (tmp_path / "_SUCCESS").is_file()
    assert not list(tmp_path.rglob("*.tmp"))


def test_publish_quality_merge_does_not_mark_invalid_global_accounting(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import (
        QualityGlobalValidation,
        publish_quality_merge,
    )

    _write_task(tmp_path, 0)
    _write_accounting_summary(
        tmp_path,
        0,
        input_sources=1,
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0,),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    with pytest.raises(
        QualityMergeError,
        match="complete Bronze manifest",
    ):
        publish_quality_merge(
            artifacts,
            quality_root=tmp_path,
            manifest_row_count=2,
            batch_size=500,
            snapshot=SNAPSHOT,
            cleaning_protocol=PROTOCOL,
            pipeline_git_commit=GIT_COMMIT,
            global_validation=QualityGlobalValidation(
                accepted_chain_count=0,
                rejected_chain_count=0,
                non_candidate_chain_count=0,
                dirty_residue_count=0,
                processing_error_count=0,
                unique_outcome_chain_count=0,
            ),
        )

    assert not (tmp_path / "_SUCCESS").exists()


def test_failed_republish_removes_previous_success_marker(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import pdbclean.quality_merge as quality_merge

    _write_task(tmp_path, 0)
    _write_accounting_summary(
        tmp_path,
        0,
        input_sources=1,
    )

    artifacts = discover_quality_task_artifacts(
        tmp_path,
        expected_task_ids=(0,),
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_pipeline_git_commit=GIT_COMMIT,
    )

    validation = quality_merge.QualityGlobalValidation(
        accepted_chain_count=0,
        rejected_chain_count=0,
        non_candidate_chain_count=0,
        dirty_residue_count=0,
        processing_error_count=0,
        unique_outcome_chain_count=0,
    )

    quality_merge.publish_quality_merge(
        artifacts,
        quality_root=tmp_path,
        manifest_row_count=1,
        batch_size=500,
        snapshot=SNAPSHOT,
        cleaning_protocol=PROTOCOL,
        pipeline_git_commit=GIT_COMMIT,
        global_validation=validation,
    )

    assert (tmp_path / "_SUCCESS").is_file()

    def fail_writer(*args, **kwargs):
        raise RuntimeError("synthetic merge failure")

    monkeypatch.setattr(
        quality_merge,
        "_write_merged_parquet_atomic",
        fail_writer,
    )

    with pytest.raises(
        RuntimeError,
        match="synthetic merge failure",
    ):
        quality_merge.publish_quality_merge(
            artifacts,
            quality_root=tmp_path,
            manifest_row_count=1,
            batch_size=500,
            snapshot=SNAPSHOT,
            cleaning_protocol=PROTOCOL,
            pipeline_git_commit=GIT_COMMIT,
            global_validation=validation,
        )

    assert not (tmp_path / "_SUCCESS").exists()


def test_merge_quality_stage_end_to_end_success(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import merge_quality_stage

    key = (
        "20260101/pub/pdb/data/structures/divided/"
        "mmCIF/aa/1aaa.cif.gz"
    )
    etag = '"etag-a"'

    _write_task(tmp_path, 0)
    _write_accounting_summary(
        tmp_path,
        0,
        input_sources=1,
    )

    publication = merge_quality_stage(
        quality_root=tmp_path,
        manifest=_manifest_table(
            [("1aaa", key, etag)]
        ),
        manifest_row_count=1,
        batch_size=500,
        snapshot=SNAPSHOT,
        cleaning_protocol=PROTOCOL,
        pipeline_git_commit=GIT_COMMIT,
    )

    assert publication.success_path == tmp_path / "_SUCCESS"
    assert publication.success_path.is_file()
    assert publication.global_summary_path.is_file()

    for path in publication.merged_paths.values():
        assert path.is_file()

    assert publication.global_summary[
        "input_source_object_count"
    ] == 1
    assert publication.global_summary["task_count"] == 1


def test_merge_quality_stage_removes_stale_success_before_discovery(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import merge_quality_stage

    key = (
        "20260101/pub/pdb/data/structures/divided/"
        "mmCIF/aa/1aaa.cif.gz"
    )

    success = tmp_path / "_SUCCESS"
    success.write_text("stale\n", encoding="utf-8")

    with pytest.raises(
        QualityMergeError,
        match="Invalid summary task set",
    ):
        merge_quality_stage(
            quality_root=tmp_path,
            manifest=_manifest_table(
                [("1aaa", key, '"etag-a"')]
            ),
            manifest_row_count=1,
            batch_size=500,
            snapshot=SNAPSHOT,
            cleaning_protocol=PROTOCOL,
            pipeline_git_commit=GIT_COMMIT,
        )

    assert not success.exists()


def test_merge_quality_stage_removes_stale_success_before_accounting(
    tmp_path: Path,
) -> None:
    from pdbclean.quality_merge import merge_quality_stage

    key = (
        "20260101/pub/pdb/data/structures/divided/"
        "mmCIF/aa/1aaa.cif.gz"
    )

    _write_task(tmp_path, 0)
    _write_accounting_summary(
        tmp_path,
        0,
        input_sources=0,
    )

    success = tmp_path / "_SUCCESS"
    success.write_text("stale\n", encoding="utf-8")

    with pytest.raises(
        QualityMergeError,
        match="input source count mismatch",
    ):
        merge_quality_stage(
            quality_root=tmp_path,
            manifest=_manifest_table(
                [("1aaa", key, '"etag-a"')]
            ),
            manifest_row_count=1,
            batch_size=500,
            snapshot=SNAPSHOT,
            cleaning_protocol=PROTOCOL,
            pipeline_git_commit=GIT_COMMIT,
        )

    assert not success.exists()
