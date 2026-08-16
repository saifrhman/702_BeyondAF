"""Validate and publish the canonical finalized Stage-5 Brain population."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.brain_runner import (
    BRAIN_TASK_SUMMARY_SCHEMA_NAME,
    BRAIN_TASK_SUMMARY_SCHEMA_VERSION,
    BRAIN_UNDEFINED_M1_REASON,
    DEFAULT_BRAIN_ROW_GROUPS_PER_TASK,
    UpstreamBRI,
    _parquet_schema_matches,
    _validated_git_commit,
    brain_task_count,
    brain_task_partition,
)
from pdbclean.schemas import (
    STAGE5_BRAIN_CHAIN_SCHEMA,
    STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA,
    STAGE5_BRAIN_UNDEFINED_CHAIN_SCHEMA,
)


class BrainFinalizeError(RuntimeError):
    """Raised when Stage-5 Brain finalization cannot proceed safely."""


BRAIN_GLOBAL_SUMMARY_SCHEMA_NAME = (
    "pdbclean_stage5_brain_global_summary"
)
BRAIN_GLOBAL_SUMMARY_SCHEMA_VERSION = "1.0"

BRAIN_SUCCESS_SCHEMA_NAME = (
    "pdbclean_stage5_brain_success"
)
BRAIN_SUCCESS_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class BrainTaskArtifacts:
    task_id: int
    chains_path: Path
    undefined_path: Path
    processing_errors_path: Path
    summary_path: Path
    summary: dict[str, Any]


@dataclass(frozen=True)
class BrainGlobalValidation:
    input_bri_chain_count: int
    brain_chain_count: int
    undefined_chain_count: int
    processing_error_count: int
    unique_input_identity_count: int
    unique_terminal_identity_count: int


@dataclass(frozen=True)
class BrainFinalizePublication:
    brain_path: Path
    undefined_path: Path
    processing_errors_path: Path
    global_summary_path: Path
    success_path: Path
    global_summary: dict[str, Any]


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise BrainFinalizeError(
            f"Required Stage-5 artifact does not exist: {path}"
        )

    try:
        value = json.loads(
            path.read_text(encoding="utf-8")
        )
    except Exception as exc:
        raise BrainFinalizeError(
            f"Cannot read Stage-5 JSON artifact {path}: {exc}"
        ) from exc

    if not isinstance(value, dict):
        raise BrainFinalizeError(
            f"Stage-5 JSON artifact is not an object: {path}"
        )

    return value


def _nonnegative_int(
    payload: dict[str, Any],
    field: str,
    *,
    task_id: int,
) -> int:
    value = payload.get(field)

    if (
        not isinstance(value, int)
        or isinstance(value, bool)
        or value < 0
    ):
        raise BrainFinalizeError(
            f"Invalid {field} in Stage-5 task {task_id}"
        )

    return value


def _validate_schema(
    path: Path,
    expected: pa.Schema,
) -> None:
    try:
        observed = pq.read_schema(path)
    except Exception as exc:
        raise BrainFinalizeError(
            f"Cannot read Stage-5 Parquet schema {path}: {exc}"
        ) from exc

    if not _parquet_schema_matches(
        observed,
        expected,
    ):
        raise BrainFinalizeError(
            f"Unexpected Stage-5 Parquet schema: {path}"
        )


def discover_brain_task_artifacts(
    brain_root: str | Path,
    *,
    upstream: UpstreamBRI,
    brain_pipeline_git_commit: str,
    row_groups_per_task: int = (
        DEFAULT_BRAIN_ROW_GROUPS_PER_TASK
    ),
) -> tuple[BrainTaskArtifacts, ...]:
    """Validate every expected Stage-5 task publication."""

    root = Path(brain_root)

    brain_pipeline_git_commit = _validated_git_commit(
        brain_pipeline_git_commit,
        field="brain_pipeline_git_commit",
    )

    parquet = pq.ParquetFile(
        upstream.bri_path
    )

    row_counts = tuple(
        parquet.metadata.row_group(i).num_rows
        for i in range(parquet.metadata.num_row_groups)
    )

    expected_task_count = brain_task_count(
        len(row_counts),
        row_groups_per_task=row_groups_per_task,
    )

    provenance = {
        "quality_pipeline_git_commit": (
            upstream.quality_pipeline_git_commit
        ),
        "geometric_validation_pipeline_git_commit": (
            upstream.geometric_validation_pipeline_git_commit
        ),
        "geometric_validation_finalizer_git_commit": (
            upstream.geometric_validation_finalizer_git_commit
        ),
        "bri_pipeline_git_commit": (
            upstream.bri_pipeline_git_commit
        ),
        "bri_finalizer_git_commit": (
            upstream.bri_finalizer_git_commit
        ),
        "brain_pipeline_git_commit": (
            brain_pipeline_git_commit
        ),
    }

    artifacts = []

    for task_id in range(expected_task_count):
        partition = brain_task_partition(
            row_counts,
            task_id=task_id,
            row_groups_per_task=row_groups_per_task,
        )

        chains = (
            root / "chains"
            / f"task_{task_id}.parquet"
        )
        undefined = (
            root / "undefined"
            / f"task_{task_id}.parquet"
        )
        errors = (
            root / "processing_errors"
            / f"task_{task_id}.parquet"
        )
        summary_path = (
            root / "summaries"
            / f"task_{task_id}.json"
        )

        for path in (
            chains,
            undefined,
            errors,
            summary_path,
        ):
            if not path.is_file():
                raise BrainFinalizeError(
                    f"Missing Stage-5 task artifact: {path}"
                )

        summary = _read_json(
            summary_path
        )

        if (
            summary.get("summary_schema_name")
            != BRAIN_TASK_SUMMARY_SCHEMA_NAME
            or summary.get("summary_schema_version")
            != BRAIN_TASK_SUMMARY_SCHEMA_VERSION
        ):
            raise BrainFinalizeError(
                f"Unexpected Stage-5 summary schema for task {task_id}"
            )

        if str(summary.get("task_id")) != str(task_id):
            raise BrainFinalizeError(
                f"Stage-5 task_id mismatch for task {task_id}"
            )

        if summary.get("task_count") != expected_task_count:
            raise BrainFinalizeError(
                f"Stage-5 task_count mismatch for task {task_id}"
            )

        expected_partition = {
            "partition_scheme": (
                "contiguous_parquet_row_groups"
            ),
            "start_row_group": (
                partition.start_row_group
            ),
            "stop_row_group": (
                partition.stop_row_group
            ),
            "row_group_count": (
                partition.row_group_count
            ),
            "input_bri_chain_count": (
                partition.input_bri_chain_count
            ),
        }

        for field, expected in expected_partition.items():
            if summary.get(field) != expected:
                raise BrainFinalizeError(
                    f"Stage-5 partition mismatch in task "
                    f"{task_id}: {field}"
                )

        if summary.get("snapshot") != upstream.snapshot:
            raise BrainFinalizeError(
                f"Stage-5 snapshot mismatch in task {task_id}"
            )

        if (
            summary.get("cleaning_protocol")
            != upstream.cleaning_protocol
        ):
            raise BrainFinalizeError(
                f"Stage-5 protocol mismatch in task {task_id}"
            )

        for field, expected in provenance.items():
            if summary.get(field) != expected:
                raise BrainFinalizeError(
                    f"Stage-5 provenance mismatch in task "
                    f"{task_id}: {field}"
                )

        if summary.get("chain_accounting_valid") is not True:
            raise BrainFinalizeError(
                f"Stage-5 chain accounting failed in task {task_id}"
            )

        brain_count = _nonnegative_int(
            summary,
            "brain_chain_count",
            task_id=task_id,
        )
        undefined_count = _nonnegative_int(
            summary,
            "undefined_chain_count",
            task_id=task_id,
        )
        error_count = _nonnegative_int(
            summary,
            "processing_error_count",
            task_id=task_id,
        )

        if (
            brain_count
            + undefined_count
            + error_count
            != partition.input_bri_chain_count
        ):
            raise BrainFinalizeError(
                f"Stage-5 terminal counts do not reconcile "
                f"for task {task_id}"
            )

        _validate_schema(
            chains,
            STAGE5_BRAIN_CHAIN_SCHEMA,
        )
        _validate_schema(
            undefined,
            STAGE5_BRAIN_UNDEFINED_CHAIN_SCHEMA,
        )
        _validate_schema(
            errors,
            STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA,
        )

        observed = (
            pq.read_metadata(chains).num_rows,
            pq.read_metadata(undefined).num_rows,
            pq.read_metadata(errors).num_rows,
        )

        if observed != (
            brain_count,
            undefined_count,
            error_count,
        ):
            raise BrainFinalizeError(
                f"Stage-5 shard row-count mismatch "
                f"for task {task_id}"
            )

        artifacts.append(
            BrainTaskArtifacts(
                task_id=task_id,
                chains_path=chains,
                undefined_path=undefined,
                processing_errors_path=errors,
                summary_path=summary_path,
                summary=summary,
            )
        )

    return tuple(artifacts)


def validate_brain_global_state(
    artifacts: tuple[BrainTaskArtifacts, ...],
    *,
    upstream: UpstreamBRI,
) -> BrainGlobalValidation:
    """Validate exact Stage-3 → Stage-5 population identity."""

    identity_columns = [
        "snapshot",
        "pdb_id",
        "model_id",
        "label_chain_id",
        "retained_residue_count",
    ]

    upstream_rows = pq.read_table(
        upstream.bri_path,
        columns=identity_columns,
    ).to_pylist()

    upstream_identity_to_m = {}

    for row in upstream_rows:
        identity = (
            row["snapshot"],
            row["pdb_id"],
            row["model_id"],
            row["label_chain_id"],
        )

        if identity in upstream_identity_to_m:
            raise BrainFinalizeError(
                f"Duplicate identity in canonical Stage-3 BRI: "
                f"{identity!r}"
            )

        upstream_identity_to_m[identity] = (
            row["retained_residue_count"]
        )

    terminal_identities = set()

    brain_count = 0
    undefined_count = 0
    error_count = 0

    for artifact in artifacts:
        brain_rows = pq.read_table(
            artifact.chains_path,
            columns=identity_columns,
        ).to_pylist()

        undefined_rows = pq.read_table(
            artifact.undefined_path,
            columns=[
                *identity_columns,
                "undefined_reason",
            ],
        ).to_pylist()

        error_rows = pq.read_table(
            artifact.processing_errors_path,
            columns=identity_columns,
        ).to_pylist()

        for kind, rows in (
            ("brain", brain_rows),
            ("undefined", undefined_rows),
            ("error", error_rows),
        ):
            for row in rows:
                identity = (
                    row["snapshot"],
                    row["pdb_id"],
                    row["model_id"],
                    row["label_chain_id"],
                )

                if identity in terminal_identities:
                    raise BrainFinalizeError(
                        "Duplicate Stage-5 terminal identity: "
                        f"{identity!r}"
                    )

                expected_m = upstream_identity_to_m.get(
                    identity
                )

                if expected_m is None:
                    raise BrainFinalizeError(
                        "Stage-5 identity absent from canonical "
                        f"Stage-3 BRI: {identity!r}"
                    )

                if row["retained_residue_count"] != expected_m:
                    raise BrainFinalizeError(
                        "Stage-5 retained residue count changed "
                        f"for identity {identity!r}"
                    )

                if (
                    kind == "brain"
                    and expected_m < 2
                ):
                    raise BrainFinalizeError(
                        "Brain-defined population contains m<2"
                    )

                if kind == "undefined":
                    if expected_m != 1:
                        raise BrainFinalizeError(
                            "Brain-undefined population contains m!=1"
                        )

                    if (
                        row["undefined_reason"]
                        != BRAIN_UNDEFINED_M1_REASON
                    ):
                        raise BrainFinalizeError(
                            "Unexpected Brain undefined reason"
                        )

                terminal_identities.add(
                    identity
                )

        brain_count += len(brain_rows)
        undefined_count += len(undefined_rows)
        error_count += len(error_rows)

    if error_count != 0:
        raise BrainFinalizeError(
            "Cannot finalize Stage-5 Brain with processing errors"
        )

    upstream_identities = set(
        upstream_identity_to_m
    )

    if terminal_identities != upstream_identities:
        missing = (
            upstream_identities
            - terminal_identities
        )
        unexpected = (
            terminal_identities
            - upstream_identities
        )

        raise BrainFinalizeError(
            "Stage-5 terminal population does not exactly "
            "match canonical Stage-3 BRI: "
            f"missing={len(missing)}, "
            f"unexpected={len(unexpected)}"
        )

    if (
        brain_count
        + undefined_count
        + error_count
        != upstream.bri_chain_count
    ):
        raise BrainFinalizeError(
            "Stage-5 global population accounting failed"
        )

    return BrainGlobalValidation(
        input_bri_chain_count=(
            upstream.bri_chain_count
        ),
        brain_chain_count=brain_count,
        undefined_chain_count=undefined_count,
        processing_error_count=error_count,
        unique_input_identity_count=len(
            upstream_identities
        ),
        unique_terminal_identity_count=len(
            terminal_identities
        ),
    )


def _write_population_atomic(
    artifacts: tuple[BrainTaskArtifacts, ...],
    *,
    artifact_name: str,
    source_attribute: str,
    schema: pa.Schema,
    output_path: Path,
    expected_rows: int,
) -> Path:
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = output_path.with_suffix(
        output_path.suffix + ".tmp"
    )

    if temporary.exists():
        temporary.unlink()

    writer = None
    rows_written = 0

    try:
        writer = pq.ParquetWriter(
            temporary,
            schema,
            compression="zstd",
            version="2.6",
        )

        for artifact in artifacts:
            source = getattr(
                artifact,
                source_attribute,
            )

            parquet = pq.ParquetFile(
                source
            )

            for batch in parquet.iter_batches(
                batch_size=256,
            ):
                table = pa.Table.from_pylist(
                    batch.to_pylist(),
                    schema=schema,
                )

                writer.write_table(
                    table
                )
                rows_written += table.num_rows

        writer.close()
        writer = None

        if rows_written != expected_rows:
            raise BrainFinalizeError(
                f"Finalized {artifact_name} row-count mismatch: "
                f"written={rows_written}, "
                f"expected={expected_rows}"
            )

        _validate_schema(
            temporary,
            schema,
        )

        if (
            pq.read_metadata(
                temporary
            ).num_rows
            != expected_rows
        ):
            raise BrainFinalizeError(
                f"Finalized {artifact_name} Parquet "
                "metadata row-count mismatch"
            )

        temporary.replace(
            output_path
        )

    except Exception:
        if writer is not None:
            try:
                writer.close()
            except Exception:
                pass

        if temporary.exists():
            temporary.unlink()

        raise

    return output_path


def _write_json_atomic(
    payload: dict[str, Any],
    path: Path,
) -> Path:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = path.with_suffix(
        path.suffix + ".tmp"
    )

    temporary.write_text(
        json.dumps(
            payload,
            sort_keys=True,
            indent=2,
            ensure_ascii=True,
        )
        + "\n",
        encoding="utf-8",
    )

    temporary.replace(
        path
    )

    return path


def finalize_brain_stage(
    *,
    brain_root: str | Path,
    upstream: UpstreamBRI,
    brain_pipeline_git_commit: str,
    finalizer_pipeline_git_commit: str,
    row_groups_per_task: int = (
        DEFAULT_BRAIN_ROW_GROUPS_PER_TASK
    ),
) -> BrainFinalizePublication:
    """Validate and publish canonical Stage-5 Brain artifacts."""

    root = Path(
        brain_root
    )

    brain_pipeline_git_commit = _validated_git_commit(
        brain_pipeline_git_commit,
        field="brain_pipeline_git_commit",
    )
    finalizer_pipeline_git_commit = _validated_git_commit(
        finalizer_pipeline_git_commit,
        field="finalizer_pipeline_git_commit",
    )

    success_path = root / "_SUCCESS"

    if success_path.exists():
        success_path.unlink()

    artifacts = discover_brain_task_artifacts(
        root,
        upstream=upstream,
        brain_pipeline_git_commit=(
            brain_pipeline_git_commit
        ),
        row_groups_per_task=(
            row_groups_per_task
        ),
    )

    validation = validate_brain_global_state(
        artifacts,
        upstream=upstream,
    )

    finalized = root / "finalized"

    brain_path = _write_population_atomic(
        artifacts,
        artifact_name="Brain-defined",
        source_attribute="chains_path",
        schema=STAGE5_BRAIN_CHAIN_SCHEMA,
        output_path=finalized / "brain.parquet",
        expected_rows=(
            validation.brain_chain_count
        ),
    )

    undefined_path = _write_population_atomic(
        artifacts,
        artifact_name="Brain-undefined",
        source_attribute="undefined_path",
        schema=STAGE5_BRAIN_UNDEFINED_CHAIN_SCHEMA,
        output_path=finalized / "undefined.parquet",
        expected_rows=(
            validation.undefined_chain_count
        ),
    )

    errors_path = _write_population_atomic(
        artifacts,
        artifact_name="Brain processing-error",
        source_attribute="processing_errors_path",
        schema=STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA,
        output_path=(
            finalized
            / "processing_errors.parquet"
        ),
        expected_rows=0,
    )

    provenance = {
        "quality_pipeline_git_commit": (
            upstream.quality_pipeline_git_commit
        ),
        "geometric_validation_pipeline_git_commit": (
            upstream.geometric_validation_pipeline_git_commit
        ),
        "geometric_validation_finalizer_git_commit": (
            upstream.geometric_validation_finalizer_git_commit
        ),
        "bri_pipeline_git_commit": (
            upstream.bri_pipeline_git_commit
        ),
        "bri_finalizer_git_commit": (
            upstream.bri_finalizer_git_commit
        ),
        "brain_pipeline_git_commit": (
            brain_pipeline_git_commit
        ),
        "finalizer_pipeline_git_commit": (
            finalizer_pipeline_git_commit
        ),
    }

    summary = {
        "summary_schema_name": (
            BRAIN_GLOBAL_SUMMARY_SCHEMA_NAME
        ),
        "summary_schema_version": (
            BRAIN_GLOBAL_SUMMARY_SCHEMA_VERSION
        ),
        "snapshot": upstream.snapshot,
        "cleaning_protocol": (
            upstream.cleaning_protocol
        ),
        **provenance,
        "task_count": len(artifacts),
        "row_groups_per_task": (
            row_groups_per_task
        ),
        "input_bri_chain_count": (
            validation.input_bri_chain_count
        ),
        "brain_chain_count": (
            validation.brain_chain_count
        ),
        "undefined_chain_count": (
            validation.undefined_chain_count
        ),
        "processing_error_count": (
            validation.processing_error_count
        ),
        "unique_input_identity_count": (
            validation.unique_input_identity_count
        ),
        "unique_terminal_identity_count": (
            validation.unique_terminal_identity_count
        ),
        "chain_accounting_valid": (
            validation.input_bri_chain_count
            == validation.brain_chain_count
            + validation.undefined_chain_count
            + validation.processing_error_count
        ),
    }

    summary_path = _write_json_atomic(
        summary,
        root / "global_summary.json",
    )

    success = {
        "success_schema_name": (
            BRAIN_SUCCESS_SCHEMA_NAME
        ),
        "success_schema_version": (
            BRAIN_SUCCESS_SCHEMA_VERSION
        ),
        "snapshot": upstream.snapshot,
        "cleaning_protocol": (
            upstream.cleaning_protocol
        ),
        **provenance,
        "task_count": len(artifacts),
        "global_summary": "global_summary.json",
        "finalized_directory": "finalized",
        "brain_population": (
            "finalized/brain.parquet"
        ),
        "undefined_population": (
            "finalized/undefined.parquet"
        ),
        "processing_errors": (
            "finalized/processing_errors.parquet"
        ),
    }

    # Completion marker is deliberately written last.
    success_path = _write_json_atomic(
        success,
        success_path,
    )

    return BrainFinalizePublication(
        brain_path=brain_path,
        undefined_path=undefined_path,
        processing_errors_path=errors_path,
        global_summary_path=summary_path,
        success_path=success_path,
        global_summary=summary,
    )
