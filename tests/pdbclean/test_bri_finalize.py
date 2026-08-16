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


def _eligible_row(
    *,
    pdb_id: str = "1abc",
    chain_id: str = "A",
    retained_ids: list[int] | None = None,
) -> dict:
    if retained_ids is None:
        retained_ids = [1, 2]

    return {
        "snapshot": SNAPSHOT,
        "pdb_id": pdb_id,
        "model_id": 1,
        "label_chain_id": chain_id,
        "auth_chain_id": chain_id,
        "entity_id": "1",
        "original_start_label_seq_id": retained_ids[0],
        "original_end_label_seq_id": retained_ids[-1],
        "retained_start_label_seq_id": retained_ids[0],
        "retained_end_label_seq_id": retained_ids[-1],
        "retained_residue_count": len(retained_ids),
        "retained_label_seq_ids": retained_ids,
        "retained_sequence": "A" * len(retained_ids),
        "terminal_trimmed": False,
        "dirty_residue_count": 0,
        "dirty_rule_ids": [],
        "source_mmcif_key": (
            f"{SNAPSHOT}/{pdb_id}.cif.gz"
        ),
        "source_etag": f"etag-{pdb_id}",
        "cleaning_protocol": PROTOCOL,
        "pipeline_git_commit": QUALITY_COMMIT,
    }


def _bri_row(
    eligible: dict,
    *,
    bri: list[list[float]] | None = None,
) -> dict:
    if bri is None:
        m = eligible["retained_residue_count"]

        bri = [
            [
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                0.0,
            ]
        ]

        for _ in range(1, m):
            bri.append(
                [
                    1.001,
                    2.002,
                    3.003,
                    4.004,
                    5.005,
                    6.006,
                    7.007,
                    8.008,
                    9.009,
                ]
            )

    return {
        **{
            key: value
            for key, value in eligible.items()
            if key != "pipeline_git_commit"
        },
        "quality_pipeline_git_commit":
            QUALITY_COMMIT,
        "geometric_validation_pipeline_git_commit":
            GEOMETRY_COMMIT,
        "geometric_validation_finalizer_git_commit":
            GEOMETRY_FINALIZER_COMMIT,
        "bri_pipeline_git_commit":
            BRI_COMMIT,
        "bri": bri,
    }


def _write_global_task(
    bri_root: Path,
    *,
    task_id: int,
    chain_rows: list[dict],
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
            chain_rows,
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

    summary_path.write_text(
        json.dumps(
            _summary(
                task_id,
                input_count=len(chain_rows),
                bri_count=len(chain_rows),
                error_count=0,
            )
        )
        + "\n",
        encoding="utf-8",
    )


def _write_eligible_population(
    path: Path,
    rows: list[dict],
) -> None:
    from pdbclean.schemas import (
        STAGE3_ELIGIBLE_CHAIN_SCHEMA,
    )

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    pq.write_table(
        pa.Table.from_pylist(
            rows,
            schema=STAGE3_ELIGIBLE_CHAIN_SCHEMA,
        ),
        path,
    )


def _validate_global(
    artifacts,
    eligible_path: Path,
):
    from pdbclean.bri_finalize import (
        validate_bri_global_state,
    )

    return validate_bri_global_state(
        artifacts,
        eligible_path=eligible_path,
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


def test_global_validation_accepts_exact_population(
    tmp_path: Path,
) -> None:
    bri_root = tmp_path / "bri"
    eligible_path = tmp_path / "eligible.parquet"

    eligible_a = _eligible_row(
        pdb_id="1aaa",
        retained_ids=[1],
    )

    eligible_b = _eligible_row(
        pdb_id="1bbb",
        retained_ids=[4, 5],
    )

    _write_eligible_population(
        eligible_path,
        [eligible_a, eligible_b],
    )

    _write_global_task(
        bri_root,
        task_id=0,
        chain_rows=[_bri_row(eligible_a)],
    )

    _write_global_task(
        bri_root,
        task_id=1,
        chain_rows=[_bri_row(eligible_b)],
    )

    artifacts = _discover(
        bri_root,
        expected_task_ids=(0, 1),
    )

    result = _validate_global(
        artifacts,
        eligible_path,
    )

    assert result.eligible_chain_count == 2
    assert result.bri_chain_count == 2
    assert result.processing_error_count == 0
    assert result.unique_eligible_identity_count == 2
    assert result.unique_bri_identity_count == 2
    assert result.minimum_retained_residue_count == 1
    assert result.maximum_retained_residue_count == 2


def test_global_validation_rejects_missing_identity(
    tmp_path: Path,
) -> None:
    bri_root = tmp_path / "bri"
    eligible_path = tmp_path / "eligible.parquet"

    eligible_a = _eligible_row(
        pdb_id="1aaa",
    )

    eligible_b = _eligible_row(
        pdb_id="1bbb",
    )

    _write_eligible_population(
        eligible_path,
        [eligible_a, eligible_b],
    )

    _write_global_task(
        bri_root,
        task_id=0,
        chain_rows=[_bri_row(eligible_a)],
    )

    artifacts = _discover(
        bri_root,
        expected_task_ids=(0,),
    )

    with pytest.raises(
        BRIFinalizeError,
        match="identity population mismatch",
    ):
        _validate_global(
            artifacts,
            eligible_path,
        )


def test_global_validation_rejects_retained_lineage_change(
    tmp_path: Path,
) -> None:
    bri_root = tmp_path / "bri"
    eligible_path = tmp_path / "eligible.parquet"

    eligible = _eligible_row()

    _write_eligible_population(
        eligible_path,
        [eligible],
    )

    bri_row = _bri_row(eligible)
    bri_row["retained_label_seq_ids"] = [1, 3]

    _write_global_task(
        bri_root,
        task_id=0,
        chain_rows=[bri_row],
    )

    artifacts = _discover(
        bri_root,
        expected_task_ids=(0,),
    )

    with pytest.raises(
        BRIFinalizeError,
        match="Retained residue-ID lineage mismatch",
    ):
        _validate_global(
            artifacts,
            eligible_path,
        )


def test_global_validation_rejects_wrong_bri_shape(
    tmp_path: Path,
) -> None:
    bri_root = tmp_path / "bri"
    eligible_path = tmp_path / "eligible.parquet"

    eligible = _eligible_row(
        retained_ids=[1, 2],
    )

    _write_eligible_population(
        eligible_path,
        [eligible],
    )

    bri_row = _bri_row(
        eligible,
        bri=[
            [
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                0.0,
            ]
        ],
    )

    _write_global_task(
        bri_root,
        task_id=0,
        chain_rows=[bri_row],
    )

    artifacts = _discover(
        bri_root,
        expected_task_ids=(0,),
    )

    with pytest.raises(
        BRIFinalizeError,
        match="BRI shape mismatch",
    ):
        _validate_global(
            artifacts,
            eligible_path,
        )


def test_global_validation_rejects_noncanonical_precision(
    tmp_path: Path,
) -> None:
    bri_root = tmp_path / "bri"
    eligible_path = tmp_path / "eligible.parquet"

    eligible = _eligible_row(
        retained_ids=[1],
    )

    _write_eligible_population(
        eligible_path,
        [eligible],
    )

    bri_row = _bri_row(
        eligible,
        bri=[
            [
                1.0004,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                0.0,
            ]
        ],
    )

    _write_global_task(
        bri_root,
        task_id=0,
        chain_rows=[bri_row],
    )

    artifacts = _discover(
        bri_root,
        expected_task_ids=(0,),
    )

    with pytest.raises(
        BRIFinalizeError,
        match="Non-canonical Definition 3.4 BRI precision",
    ):
        _validate_global(
            artifacts,
            eligible_path,
        )


def test_global_validation_rejects_first_row_xa_nonzero(
    tmp_path: Path,
) -> None:
    bri_root = tmp_path / "bri"
    eligible_path = tmp_path / "eligible.parquet"

    eligible = _eligible_row(
        retained_ids=[1],
    )

    _write_eligible_population(
        eligible_path,
        [eligible],
    )

    bri_row = _bri_row(
        eligible,
        bri=[
            [
                1.0,
                0.0,
                0.0,
                0.001,
                0.0,
                0.0,
                1.0,
                1.0,
                0.0,
            ]
        ],
    )

    _write_global_task(
        bri_root,
        task_id=0,
        chain_rows=[bri_row],
    )

    artifacts = _discover(
        bri_root,
        expected_task_ids=(0,),
    )

    with pytest.raises(
        BRIFinalizeError,
        match="first-row zero structure mismatch",
    ):
        _validate_global(
            artifacts,
            eligible_path,
        )


def test_finalize_publishes_canonical_bri_population_and_success(
    tmp_path: Path,
) -> None:
    from pdbclean.bri_finalize import (
        finalize_bri_stage,
    )

    bri_root = tmp_path / "bri"
    eligible_path = tmp_path / "eligible.parquet"

    eligible_a = _eligible_row(
        pdb_id="1aaa",
        retained_ids=[1],
    )
    eligible_b = _eligible_row(
        pdb_id="1bbb",
        retained_ids=[4, 5],
    )

    _write_eligible_population(
        eligible_path,
        [eligible_a, eligible_b],
    )

    _write_global_task(
        bri_root,
        task_id=0,
        chain_rows=[_bri_row(eligible_a)],
    )
    _write_global_task(
        bri_root,
        task_id=1,
        chain_rows=[_bri_row(eligible_b)],
    )

    publication = finalize_bri_stage(
        bri_root=bri_root,
        eligible_path=eligible_path,
        manifest_row_count=2,
        batch_size=1,
        snapshot=SNAPSHOT,
        cleaning_protocol=PROTOCOL,
        quality_pipeline_git_commit=QUALITY_COMMIT,
        geometric_validation_pipeline_git_commit=(
            GEOMETRY_COMMIT
        ),
        geometric_validation_finalizer_git_commit=(
            GEOMETRY_FINALIZER_COMMIT
        ),
        bri_pipeline_git_commit=BRI_COMMIT,
        finalizer_pipeline_git_commit="f" * 40,
    )

    assert publication.bri_path.is_file()
    assert publication.global_summary_path.is_file()
    assert publication.success_path.is_file()

    finalized = pq.read_table(
        publication.bri_path
    )

    assert finalized.num_rows == 2
    assert finalized.schema.metadata == (
        STAGE3_BRI_CHAIN_SCHEMA.metadata
    )

    identities = [
        (
            row["pdb_id"],
            row["label_chain_id"],
        )
        for row in finalized.to_pylist()
    ]

    assert identities == [
        ("1aaa", "A"),
        ("1bbb", "A"),
    ]

    summary = publication.global_summary

    assert summary["task_count"] == 2
    assert summary["input_eligible_chain_count"] == 2
    assert summary["bri_chain_count"] == 2
    assert summary["processing_error_count"] == 0
    assert summary["unique_eligible_identity_count"] == 2
    assert summary["unique_bri_identity_count"] == 2
    assert summary["minimum_retained_residue_count"] == 1
    assert summary["maximum_retained_residue_count"] == 2
    assert summary["chain_accounting_valid"] is True

    success = json.loads(
        publication.success_path.read_text(
            encoding="utf-8"
        )
    )

    assert success["success_schema_name"] == (
        "pdbclean_stage3_bri_success"
    )
    assert success["bri_population"] == (
        "finalized/bri.parquet"
    )
    assert success[
        "finalizer_pipeline_git_commit"
    ] == "f" * 40


def test_finalize_removes_stale_success_before_failure(
    tmp_path: Path,
) -> None:
    from pdbclean.bri_finalize import (
        finalize_bri_stage,
    )

    bri_root = tmp_path / "bri"
    eligible_path = tmp_path / "eligible.parquet"

    eligible = _eligible_row(
        pdb_id="1aaa",
    )

    _write_eligible_population(
        eligible_path,
        [eligible],
    )

    bri_row = _bri_row(eligible)
    bri_row["source_etag"] = "wrong-etag"

    _write_global_task(
        bri_root,
        task_id=0,
        chain_rows=[bri_row],
    )

    bri_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    success_path = bri_root / "_SUCCESS"
    success_path.write_text(
        "stale\n",
        encoding="utf-8",
    )

    with pytest.raises(
        BRIFinalizeError,
        match="source_etag",
    ):
        finalize_bri_stage(
            bri_root=bri_root,
            eligible_path=eligible_path,
            manifest_row_count=1,
            batch_size=500,
            snapshot=SNAPSHOT,
            cleaning_protocol=PROTOCOL,
            quality_pipeline_git_commit=(
                QUALITY_COMMIT
            ),
            geometric_validation_pipeline_git_commit=(
                GEOMETRY_COMMIT
            ),
            geometric_validation_finalizer_git_commit=(
                GEOMETRY_FINALIZER_COMMIT
            ),
            bri_pipeline_git_commit=(
                BRI_COMMIT
            ),
            finalizer_pipeline_git_commit="f" * 40,
        )

    assert not success_path.exists()


def test_global_summary_rejects_source_accounting_mismatch(
    tmp_path: Path,
) -> None:
    from pdbclean.bri_finalize import (
        build_bri_global_summary,
    )

    bri_root = tmp_path / "bri"
    eligible_path = tmp_path / "eligible.parquet"

    eligible = _eligible_row(
        pdb_id="1aaa",
    )

    _write_eligible_population(
        eligible_path,
        [eligible],
    )

    _write_global_task(
        bri_root,
        task_id=0,
        chain_rows=[_bri_row(eligible)],
    )

    artifacts = _discover(
        bri_root,
        expected_task_ids=(0,),
    )

    result = _validate_global(
        artifacts,
        eligible_path,
    )

    artifacts[0].summary[
        "manifest_source_object_count"
    ] = 1
    artifacts[0].summary[
        "relevant_source_object_count"
    ] = 1

    with pytest.raises(
        BRIFinalizeError,
        match="relevant/downloaded/parsed",
    ):
        build_bri_global_summary(
            artifacts,
            snapshot=SNAPSHOT,
            cleaning_protocol=PROTOCOL,
            quality_pipeline_git_commit=(
                QUALITY_COMMIT
            ),
            geometric_validation_pipeline_git_commit=(
                GEOMETRY_COMMIT
            ),
            geometric_validation_finalizer_git_commit=(
                GEOMETRY_FINALIZER_COMMIT
            ),
            bri_pipeline_git_commit=(
                BRI_COMMIT
            ),
            finalizer_pipeline_git_commit="f" * 40,
            global_validation=result,
        )
