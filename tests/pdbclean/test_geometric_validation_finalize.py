"""Tests for Step-2 geometric-validation finalization."""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from pdbclean.geometric_validation_finalize import (
    GeometricValidationFinalizeError,
    discover_geometric_validation_task_artifacts,
    validate_geometric_validation_task_accounting,
)
from pdbclean.schemas import (
    GEOMETRIC_VALIDATION_AUDIT_SCHEMA,
    GEOMETRIC_VALIDATION_PROCESSING_ERROR_SCHEMA,
    GOLD_ACCEPTED_CHAIN_SCHEMA,
)


SNAPSHOT = "20310415"
PROTOCOL = "protocol3.2-comp702-v1"
QUALITY_COMMIT = "q" * 40
GEOMETRY_COMMIT = "g" * 40


def _summary(
    task_id: int,
    *,
    input_count: int = 0,
    audit_count: int = 0,
    passed_count: int = 0,
    violated_count: int = 0,
    error_count: int = 0,
) -> dict:
    return {
        "summary_schema_name": (
            "pdbclean_geometric_validation_task_summary"
        ),
        "summary_schema_version": "1.0",
        "task_id": str(task_id),
        "snapshot": SNAPSHOT,
        "cleaning_protocol": PROTOCOL,
        "quality_pipeline_git_commit": QUALITY_COMMIT,
        "geometric_validation_pipeline_git_commit": (
            GEOMETRY_COMMIT
        ),
        "configured_minimum_backbone_distance_angstrom": 0.01,
        "configured_minimum_triangle_angle_degrees": 3.0,
        "input_accepted_chain_count": input_count,
        "audit_chain_count": audit_count,
        "geometric_passed_chain_count": passed_count,
        "geometric_violated_chain_count": violated_count,
        "processing_error_count": error_count,
        "chain_accounting_valid": True,
    }


def _write_task(
    quality_root: Path,
    geometry_root: Path,
    task_id: int,
    *,
    accepted_rows: list[dict] | None = None,
    audit_rows: list[dict] | None = None,
    error_rows: list[dict] | None = None,
    summary: dict | None = None,
) -> None:
    accepted_rows = [] if accepted_rows is None else accepted_rows
    audit_rows = [] if audit_rows is None else audit_rows
    error_rows = [] if error_rows is None else error_rows

    accepted_path = (
        quality_root / "accepted" / f"task_{task_id}.parquet"
    )
    audit_path = (
        geometry_root / "audit" / f"task_{task_id}.parquet"
    )
    error_path = (
        geometry_root / "errors" / f"task_{task_id}.parquet"
    )
    summary_path = (
        geometry_root / "summaries" / f"task_{task_id}.json"
    )

    for path in (
        accepted_path,
        audit_path,
        error_path,
        summary_path,
    ):
        path.parent.mkdir(parents=True, exist_ok=True)

    pq.write_table(
        pa.Table.from_pylist(
            accepted_rows,
            schema=GOLD_ACCEPTED_CHAIN_SCHEMA,
        ),
        accepted_path,
    )

    pq.write_table(
        pa.Table.from_pylist(
            audit_rows,
            schema=GEOMETRIC_VALIDATION_AUDIT_SCHEMA,
        ),
        audit_path,
    )

    pq.write_table(
        pa.Table.from_pylist(
            error_rows,
            schema=GEOMETRIC_VALIDATION_PROCESSING_ERROR_SCHEMA,
        ),
        error_path,
    )

    if summary is None:
        summary = _summary(task_id)

    summary_path.write_text(
        json.dumps(summary) + "\n",
        encoding="utf-8",
    )


def _discover(
    quality_root: Path,
    geometry_root: Path,
    *,
    expected_task_ids: tuple[int, ...],
):
    return discover_geometric_validation_task_artifacts(
        quality_root,
        geometry_root,
        expected_task_ids=expected_task_ids,
        expected_snapshot=SNAPSHOT,
        expected_cleaning_protocol=PROTOCOL,
        expected_quality_pipeline_git_commit=QUALITY_COMMIT,
        expected_geometric_validation_pipeline_git_commit=(
            GEOMETRY_COMMIT
        ),
        expected_minimum_backbone_distance_angstrom=0.01,
        expected_minimum_triangle_angle_degrees=3.0,
    )


def test_discovery_accepts_complete_task_set(
    tmp_path: Path,
) -> None:
    quality = tmp_path / "quality"
    geometry = tmp_path / "geometry"

    _write_task(quality, geometry, 0)
    _write_task(quality, geometry, 1)

    artifacts = _discover(
        quality,
        geometry,
        expected_task_ids=(0, 1),
    )

    assert tuple(
        artifact.task_id for artifact in artifacts
    ) == (0, 1)


def test_discovery_rejects_missing_audit_shard(
    tmp_path: Path,
) -> None:
    quality = tmp_path / "quality"
    geometry = tmp_path / "geometry"

    _write_task(quality, geometry, 0)
    (geometry / "audit" / "task_0.parquet").unlink()

    with pytest.raises(
        GeometricValidationFinalizeError,
        match="audit shard.*missing=\\[0\\]",
    ):
        _discover(
            quality,
            geometry,
            expected_task_ids=(0,),
        )


def test_discovery_rejects_temporary_file(
    tmp_path: Path,
) -> None:
    quality = tmp_path / "quality"
    geometry = tmp_path / "geometry"

    _write_task(quality, geometry, 0)

    temporary = geometry / "audit" / "task_0.parquet.tmp"
    temporary.write_text("partial", encoding="utf-8")

    with pytest.raises(
        GeometricValidationFinalizeError,
        match="Temporary upstream/finalization input files remain",
    ):
        _discover(
            quality,
            geometry,
            expected_task_ids=(0,),
        )


def test_discovery_rejects_summary_provenance_mismatch(
    tmp_path: Path,
) -> None:
    quality = tmp_path / "quality"
    geometry = tmp_path / "geometry"

    summary = _summary(0)
    summary["geometric_validation_pipeline_git_commit"] = "x" * 40

    _write_task(
        quality,
        geometry,
        0,
        summary=summary,
    )

    with pytest.raises(
        GeometricValidationFinalizeError,
        match="Geometry producer Git commit mismatch",
    ):
        _discover(
            quality,
            geometry,
            expected_task_ids=(0,),
        )


def test_task_accounting_recomputes_parquet_counts(
    tmp_path: Path,
) -> None:
    quality = tmp_path / "quality"
    geometry = tmp_path / "geometry"

    summary = _summary(
        0,
        input_count=1,
        audit_count=0,
        passed_count=0,
        violated_count=0,
        error_count=1,
    )

    _write_task(
        quality,
        geometry,
        0,
        summary=summary,
    )

    artifacts = _discover(
        quality,
        geometry,
        expected_task_ids=(0,),
    )

    with pytest.raises(
        GeometricValidationFinalizeError,
        match="Accepted row-count mismatch",
    ):
        validate_geometric_validation_task_accounting(
            artifacts
        )


def test_task_accounting_recomputes_terminal_invariant(
    tmp_path: Path,
) -> None:
    quality = tmp_path / "quality"
    geometry = tmp_path / "geometry"

    summary = _summary(
        0,
        input_count=0,
        audit_count=0,
        passed_count=1,
        violated_count=0,
        error_count=0,
    )

    _write_task(
        quality,
        geometry,
        0,
        summary=summary,
    )

    artifacts = _discover(
        quality,
        geometry,
        expected_task_ids=(0,),
    )

    with pytest.raises(
        GeometricValidationFinalizeError,
        match="Audit pass/violation accounting mismatch",
    ):
        validate_geometric_validation_task_accounting(
            artifacts
        )


def _accepted_row(
    *,
    pdb_id: str,
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
        "source_mmcif_key": f"{SNAPSHOT}/{pdb_id}.cif.gz",
        "source_etag": f"etag-{pdb_id}",
        "cleaning_protocol": PROTOCOL,
        "pipeline_git_commit": QUALITY_COMMIT,
    }


def _audit_row(
    accepted: dict,
    *,
    passed: bool,
) -> dict:
    if passed:
        violation_count = 0
        violation_types = []
        violation_residue_ids = []
        violation_details = []
        minimum_angle = 90.0
    else:
        violation_count = 1
        violation_types = [
            "triangle_angle_below_minimum"
        ]
        violation_residue_ids = [
            accepted["retained_label_seq_ids"][0]
        ]
        violation_details = [
            "vertex=N: 2.5 < 3"
        ]
        minimum_angle = 2.5

    return {
        "snapshot": accepted["snapshot"],
        "pdb_id": accepted["pdb_id"],
        "model_id": accepted["model_id"],
        "label_chain_id": accepted["label_chain_id"],
        "retained_residue_count": (
            accepted["retained_residue_count"]
        ),
        "retained_label_seq_ids": (
            accepted["retained_label_seq_ids"]
        ),
        "source_mmcif_key": accepted["source_mmcif_key"],
        "source_etag": accepted["source_etag"],
        "cleaning_protocol": accepted["cleaning_protocol"],
        "quality_pipeline_git_commit": QUALITY_COMMIT,
        "geometric_validation_pipeline_git_commit": (
            GEOMETRY_COMMIT
        ),
        "configured_minimum_backbone_distance_angstrom": 0.01,
        "configured_minimum_triangle_angle_degrees": 3.0,
        "passed": passed,
        "minimum_observed_backbone_distance_angstrom": 1.3,
        "minimum_observed_triangle_angle_degrees": minimum_angle,
        "minimum_observed_basis_h_norm_angstrom": 0.5,
        "violation_count": violation_count,
        "violation_types": violation_types,
        "violation_residue_ids": violation_residue_ids,
        "violation_details": violation_details,
    }


def _complete_two_chain_task(
    quality: Path,
    geometry: Path,
) -> tuple[dict, dict]:
    passed = _accepted_row(pdb_id="1aaa")
    failed = _accepted_row(pdb_id="1bbb")

    _write_task(
        quality,
        geometry,
        0,
        accepted_rows=[passed, failed],
        audit_rows=[
            _audit_row(passed, passed=True),
            _audit_row(failed, passed=False),
        ],
        summary=_summary(
            0,
            input_count=2,
            audit_count=2,
            passed_count=1,
            violated_count=1,
            error_count=0,
        )
        | {
            "relevant_source_object_count": 2,
            "downloaded_source_object_count": 2,
            "parsed_source_object_count": 2,
        },
    )

    return passed, failed


def test_global_validation_partitions_eligible_and_quarantined(
    tmp_path: Path,
) -> None:
    from pdbclean.geometric_validation_finalize import (
        validate_geometric_validation_global_state,
    )

    quality = tmp_path / "quality"
    geometry = tmp_path / "geometry"

    _complete_two_chain_task(
        quality,
        geometry,
    )

    artifacts = _discover(
        quality,
        geometry,
        expected_task_ids=(0,),
    )

    validate_geometric_validation_task_accounting(
        artifacts
    )

    result = validate_geometric_validation_global_state(
        artifacts
    )

    assert result.input_accepted_chain_count == 2
    assert result.audit_chain_count == 2
    assert result.eligible_chain_count == 1
    assert result.quarantined_chain_count == 1
    assert result.processing_error_count == 0
    assert result.unique_chain_identity_count == 2
    assert result.violation_event_count == 1
    assert result.violations_by_type == {
        "triangle_angle_below_minimum": 1,
    }


def test_global_validation_rejects_exact_lineage_mismatch(
    tmp_path: Path,
) -> None:
    from pdbclean.geometric_validation_finalize import (
        validate_geometric_validation_global_state,
    )

    quality = tmp_path / "quality"
    geometry = tmp_path / "geometry"

    passed = _accepted_row(pdb_id="1aaa")
    audit = _audit_row(passed, passed=True)
    audit["retained_label_seq_ids"] = [1, 3]

    _write_task(
        quality,
        geometry,
        0,
        accepted_rows=[passed],
        audit_rows=[audit],
        summary=_summary(
            0,
            input_count=1,
            audit_count=1,
            passed_count=1,
        ),
    )

    artifacts = _discover(
        quality,
        geometry,
        expected_task_ids=(0,),
    )

    with pytest.raises(
        GeometricValidationFinalizeError,
        match="retained_label_seq_ids",
    ):
        validate_geometric_validation_global_state(
            artifacts
        )


def test_global_validation_rejects_exact_etag_mismatch(
    tmp_path: Path,
) -> None:
    from pdbclean.geometric_validation_finalize import (
        validate_geometric_validation_global_state,
    )

    quality = tmp_path / "quality"
    geometry = tmp_path / "geometry"

    passed = _accepted_row(pdb_id="1aaa")
    audit = _audit_row(passed, passed=True)
    audit["source_etag"] = "different-etag"

    _write_task(
        quality,
        geometry,
        0,
        accepted_rows=[passed],
        audit_rows=[audit],
        summary=_summary(
            0,
            input_count=1,
            audit_count=1,
            passed_count=1,
        ),
    )

    artifacts = _discover(
        quality,
        geometry,
        expected_task_ids=(0,),
    )

    with pytest.raises(
        GeometricValidationFinalizeError,
        match="source_etag",
    ):
        validate_geometric_validation_global_state(
            artifacts
        )


def test_finalize_publishes_stage3_population_and_success_last(
    tmp_path: Path,
) -> None:
    from pdbclean.geometric_validation_finalize import (
        finalize_geometric_validation_stage,
    )

    quality = tmp_path / "quality"
    geometry = tmp_path / "geometry"

    passed, failed = _complete_two_chain_task(
        quality,
        geometry,
    )

    publication = finalize_geometric_validation_stage(
        quality_root=quality,
        geometric_validation_root=geometry,
        manifest_row_count=1,
        batch_size=500,
        snapshot=SNAPSHOT,
        cleaning_protocol=PROTOCOL,
        quality_pipeline_git_commit=QUALITY_COMMIT,
        geometric_validation_pipeline_git_commit=(
            GEOMETRY_COMMIT
        ),
        finalizer_pipeline_git_commit="f" * 40,
        minimum_backbone_distance_angstrom=0.01,
        minimum_triangle_angle_degrees=3.0,
    )

    assert publication.eligible_path.is_file()
    assert publication.quarantined_path.is_file()
    assert publication.global_summary_path.is_file()
    assert publication.success_path.is_file()

    eligible = pq.read_table(
        publication.eligible_path
    ).to_pylist()

    quarantined = pq.read_table(
        publication.quarantined_path
    ).to_pylist()

    assert eligible == [passed]
    assert quarantined == [failed]

    assert publication.global_summary[
        "input_accepted_chain_count"
    ] == 2
    assert publication.global_summary[
        "eligible_chain_count"
    ] == 1
    assert publication.global_summary[
        "quarantined_chain_count"
    ] == 1

    success = json.loads(
        publication.success_path.read_text(
            encoding="utf-8"
        )
    )

    assert success["success_schema_name"] == (
        "pdbclean_geometric_validation_success"
    )
    assert success["eligible_population"] == (
        "finalized/eligible.parquet"
    )


def test_finalize_removes_stale_success_before_failure(
    tmp_path: Path,
) -> None:
    from pdbclean.geometric_validation_finalize import (
        finalize_geometric_validation_stage,
    )

    quality = tmp_path / "quality"
    geometry = tmp_path / "geometry"

    passed = _accepted_row(pdb_id="1aaa")
    audit = _audit_row(passed, passed=True)
    audit["source_etag"] = "wrong"

    _write_task(
        quality,
        geometry,
        0,
        accepted_rows=[passed],
        audit_rows=[audit],
        summary=_summary(
            0,
            input_count=1,
            audit_count=1,
            passed_count=1,
        ),
    )

    geometry.mkdir(parents=True, exist_ok=True)

    success = geometry / "_SUCCESS"
    success.write_text(
        "stale\n",
        encoding="utf-8",
    )

    with pytest.raises(
        GeometricValidationFinalizeError,
        match="source_etag",
    ):
        finalize_geometric_validation_stage(
            quality_root=quality,
            geometric_validation_root=geometry,
            manifest_row_count=1,
            batch_size=500,
            snapshot=SNAPSHOT,
            cleaning_protocol=PROTOCOL,
            quality_pipeline_git_commit=QUALITY_COMMIT,
            geometric_validation_pipeline_git_commit=(
                GEOMETRY_COMMIT
            ),
            finalizer_pipeline_git_commit="f" * 40,
            minimum_backbone_distance_angstrom=0.01,
            minimum_triangle_angle_degrees=3.0,
        )

    assert not success.exists()
