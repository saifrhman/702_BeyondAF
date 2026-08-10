"""Merge and globally validate distributed Protocol 3.2 quality outputs."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.gold import (
    QUALITY_TASK_SUMMARY_SCHEMA_NAME,
    QUALITY_TASK_SUMMARY_SCHEMA_VERSION,
)
from pdbclean.manifest import (
    ManifestError,
    manifest_partition_count,
    validate_manifest_table,
)
from pdbclean.schemas import (
    GOLD_ACCEPTED_CHAIN_SCHEMA,
    GOLD_DIRTY_RESIDUE_SCHEMA,
    GOLD_NON_CANDIDATE_CHAIN_SCHEMA,
    GOLD_PROCESSING_ERROR_SCHEMA,
    GOLD_REJECTED_CHAIN_SCHEMA,
)


class QualityMergeError(RuntimeError):
    """Raised when a distributed quality stage cannot be validated safely."""


QUALITY_SHARD_SCHEMAS: dict[str, pa.Schema] = {
    "accepted": GOLD_ACCEPTED_CHAIN_SCHEMA,
    "rejected": GOLD_REJECTED_CHAIN_SCHEMA,
    "non_candidates": GOLD_NON_CANDIDATE_CHAIN_SCHEMA,
    "dirty_residues": GOLD_DIRTY_RESIDUE_SCHEMA,
    "errors": GOLD_PROCESSING_ERROR_SCHEMA,
}

_TASK_SUMMARY_RE = re.compile(r"^task_(\d+)\.json$")
_TASK_SHARD_RE = re.compile(r"^task_(\d+)\.parquet$")


@dataclass(frozen=True)
class QualityTaskArtifacts:
    """Validated filesystem locations and summary for one quality task."""

    task_id: int
    summary_path: Path
    shard_paths: dict[str, Path]
    summary: dict[str, Any]


def expected_quality_task_ids(
    manifest_row_count: int,
    batch_size: int,
) -> tuple[int, ...]:
    """Derive expected zero-based quality task IDs dynamically."""

    partition_count = manifest_partition_count(
        manifest_row_count,
        batch_size,
    )
    return tuple(range(partition_count))


def _parse_task_id(
    path: Path,
    pattern: re.Pattern[str],
    *,
    artifact_type: str,
) -> int:
    match = pattern.fullmatch(path.name)

    if match is None:
        raise QualityMergeError(
            f"Unexpected {artifact_type} filename: {path.name!r}"
        )

    return int(match.group(1))


def _discover_task_ids(
    directory: Path,
    *,
    suffix: str,
    pattern: re.Pattern[str],
    artifact_type: str,
) -> set[int]:
    if not directory.is_dir():
        return set()

    task_ids: set[int] = set()

    for path in sorted(directory.glob(f"*{suffix}")):
        task_id = _parse_task_id(
            path,
            pattern,
            artifact_type=artifact_type,
        )

        if task_id in task_ids:
            raise QualityMergeError(
                f"Duplicate {artifact_type} task ID: {task_id}"
            )

        task_ids.add(task_id)

    return task_ids


def _validate_exact_task_ids(
    observed: set[int],
    expected: set[int],
    *,
    artifact_type: str,
) -> None:
    missing = sorted(expected - observed)
    unexpected = sorted(observed - expected)

    if missing or unexpected:
        pieces = []

        if missing:
            pieces.append(f"missing={missing}")

        if unexpected:
            pieces.append(f"unexpected={unexpected}")

        raise QualityMergeError(
            f"Invalid {artifact_type} task set: "
            + ", ".join(pieces)
        )


def _load_task_summary(
    path: Path,
    *,
    expected_task_id: int,
    expected_snapshot: str,
    expected_cleaning_protocol: str,
    expected_pipeline_git_commit: str,
) -> dict[str, Any]:
    try:
        summary = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise QualityMergeError(
            f"Cannot read quality-task summary {path}: {exc}"
        ) from exc

    if not isinstance(summary, dict):
        raise QualityMergeError(
            f"Quality-task summary must contain a JSON object: {path}"
        )

    if (
        summary.get("summary_schema_name")
        != QUALITY_TASK_SUMMARY_SCHEMA_NAME
    ):
        raise QualityMergeError(
            f"Unexpected summary schema name in {path}"
        )

    if (
        summary.get("summary_schema_version")
        != QUALITY_TASK_SUMMARY_SCHEMA_VERSION
    ):
        raise QualityMergeError(
            f"Unexpected summary schema version in {path}"
        )

    if str(summary.get("task_id", "")) != str(expected_task_id):
        raise QualityMergeError(
            f"Summary task ID does not match filename for task "
            f"{expected_task_id}"
        )

    if summary.get("snapshot") != expected_snapshot:
        raise QualityMergeError(
            f"Snapshot mismatch in task {expected_task_id} summary"
        )

    if (
        summary.get("cleaning_protocol")
        != expected_cleaning_protocol
    ):
        raise QualityMergeError(
            f"Cleaning protocol mismatch in task "
            f"{expected_task_id} summary"
        )

    if (
        summary.get("pipeline_git_commit")
        != expected_pipeline_git_commit
    ):
        raise QualityMergeError(
            f"Git commit mismatch in task {expected_task_id} summary"
        )

    if summary.get("source_object_accounting_valid") is not True:
        raise QualityMergeError(
            f"Source accounting is invalid for task {expected_task_id}"
        )

    if summary.get("selected_chain_accounting_valid") is not True:
        raise QualityMergeError(
            f"Selected-chain accounting is invalid for task "
            f"{expected_task_id}"
        )

    return summary


def _parquet_field_matches(
    observed: pa.Field,
    expected: pa.Field,
) -> bool:
    """Compare Arrow fields while tolerating Parquet list-child renaming."""

    if (
        observed.name != expected.name
        or observed.nullable != expected.nullable
        or observed.metadata != expected.metadata
    ):
        return False

    observed_type = observed.type
    expected_type = expected.type

    if pa.types.is_list(observed_type) and pa.types.is_list(expected_type):
        observed_value = observed_type.value_field
        expected_value = expected_type.value_field

        # PyArrow's Parquet round trip normalises the list child field
        # name from "item" to "element". The child name is not part of
        # the published PDBClean schema contract; its type, nullability,
        # and metadata remain strict.
        return (
            observed_value.type == expected_value.type
            and observed_value.nullable == expected_value.nullable
            and observed_value.metadata == expected_value.metadata
        )

    return observed_type == expected_type


def _validate_shard_schema(
    path: Path,
    expected_schema: pa.Schema,
) -> None:
    try:
        observed_schema = pq.read_schema(path)
    except Exception as exc:
        raise QualityMergeError(
            f"Cannot read Parquet schema from {path}: {exc}"
        ) from exc

    if observed_schema.metadata != expected_schema.metadata:
        raise QualityMergeError(
            f"Unexpected Parquet schema metadata for {path}"
        )

    if len(observed_schema) != len(expected_schema):
        raise QualityMergeError(
            f"Unexpected Parquet field count for {path}"
        )

    if any(
        not _parquet_field_matches(observed, expected)
        for observed, expected in zip(
            observed_schema,
            expected_schema,
            strict=True,
        )
    ):
        raise QualityMergeError(
            f"Unexpected Parquet schema for {path}"
        )


def discover_quality_task_artifacts(
    quality_root: str | Path,
    *,
    expected_task_ids: Iterable[int],
    expected_snapshot: str,
    expected_cleaning_protocol: str,
    expected_pipeline_git_commit: str,
) -> tuple[QualityTaskArtifacts, ...]:
    """Discover and validate the filesystem contract for all quality tasks."""

    root = Path(quality_root)
    expected_ids = tuple(sorted(set(expected_task_ids)))
    expected_set = set(expected_ids)

    if expected_ids != tuple(range(len(expected_ids))):
        raise QualityMergeError(
            "Expected task IDs must be contiguous and zero-based"
        )

    temporary_files = sorted(root.rglob("*.tmp"))

    if temporary_files:
        raise QualityMergeError(
            "Temporary quality-stage files remain: "
            + ", ".join(str(path) for path in temporary_files)
        )

    summary_dir = root / "summaries"

    observed_summary_ids = _discover_task_ids(
        summary_dir,
        suffix=".json",
        pattern=_TASK_SUMMARY_RE,
        artifact_type="summary",
    )

    _validate_exact_task_ids(
        observed_summary_ids,
        expected_set,
        artifact_type="summary",
    )

    for shard_name in QUALITY_SHARD_SCHEMAS:
        observed_shard_ids = _discover_task_ids(
            root / shard_name,
            suffix=".parquet",
            pattern=_TASK_SHARD_RE,
            artifact_type=f"{shard_name} shard",
        )

        _validate_exact_task_ids(
            observed_shard_ids,
            expected_set,
            artifact_type=f"{shard_name} shard",
        )

    artifacts: list[QualityTaskArtifacts] = []

    for task_id in expected_ids:
        summary_path = (
            summary_dir / f"task_{task_id}.json"
        )

        summary = _load_task_summary(
            summary_path,
            expected_task_id=task_id,
            expected_snapshot=expected_snapshot,
            expected_cleaning_protocol=(
                expected_cleaning_protocol
            ),
            expected_pipeline_git_commit=(
                expected_pipeline_git_commit
            ),
        )

        shard_paths: dict[str, Path] = {}

        for shard_name, expected_schema in (
            QUALITY_SHARD_SCHEMAS.items()
        ):
            shard_path = (
                root
                / shard_name
                / f"task_{task_id}.parquet"
            )

            _validate_shard_schema(
                shard_path,
                expected_schema,
            )

            shard_paths[shard_name] = shard_path

        artifacts.append(
            QualityTaskArtifacts(
                task_id=task_id,
                summary_path=summary_path,
                shard_paths=shard_paths,
                summary=summary,
            )
        )

    return tuple(artifacts)


@dataclass(frozen=True)
class QualityTaskValidation:
    """Validated accounting facts for one completed quality task."""

    task_id: int
    expected_input_source_object_count: int
    shard_row_counts: dict[str, int]


_SUMMARY_COUNT_FIELDS = (
    "input_source_object_count",
    "successful_source_object_count",
    "failed_source_object_count",
    "parsed_silver_chain_count",
    "selected_silver_chain_count",
    "candidate_entry_count",
    "candidate_chain_count",
    "non_candidate_chain_count",
    "accepted_chain_count",
    "accepted_trimmed_chain_count",
    "rejected_chain_count",
    "dirty_residue_count",
    "processing_error_count",
    "chain_level_processing_error_count",
    "source_entry_processing_error_count",
    "total_gold_record_count",
)


def _summary_nonnegative_int(
    summary: dict[str, Any],
    field: str,
    *,
    task_id: int,
) -> int:
    value = summary.get(field)

    if (
        not isinstance(value, int)
        or isinstance(value, bool)
        or value < 0
    ):
        raise QualityMergeError(
            f"Task {task_id} summary field {field!r} "
            "must be a non-negative integer"
        )

    return value


def _expected_task_source_count(
    manifest_row_count: int,
    batch_size: int,
    task_id: int,
) -> int:
    if (
        not isinstance(manifest_row_count, int)
        or isinstance(manifest_row_count, bool)
        or manifest_row_count < 0
    ):
        raise QualityMergeError(
            "manifest_row_count must be a non-negative integer"
        )

    expected_ids = expected_quality_task_ids(
        manifest_row_count,
        batch_size,
    )

    if task_id not in expected_ids:
        raise QualityMergeError(
            f"Task {task_id} is outside the expected manifest partitions"
        )

    start = task_id * batch_size
    return min(batch_size, manifest_row_count - start)


def _parquet_row_count(path: Path) -> int:
    try:
        metadata = pq.ParquetFile(path).metadata
    except Exception as exc:
        raise QualityMergeError(
            f"Cannot read Parquet metadata from {path}: {exc}"
        ) from exc

    return metadata.num_rows


def _accepted_trimmed_count(path: Path) -> int:
    try:
        values = pq.read_table(
            path,
            columns=["terminal_trimmed"],
        )["terminal_trimmed"].to_pylist()
    except Exception as exc:
        raise QualityMergeError(
            f"Cannot read terminal_trimmed values from {path}: {exc}"
        ) from exc

    return sum(bool(value) for value in values)


def _processing_error_level_counts(
    path: Path,
) -> tuple[int, int]:
    try:
        table = pq.read_table(
            path,
            columns=["model_id", "label_chain_id"],
        )
    except Exception as exc:
        raise QualityMergeError(
            f"Cannot read processing-error identities from {path}: {exc}"
        ) from exc

    model_ids = table["model_id"].to_pylist()
    chain_ids = table["label_chain_id"].to_pylist()

    chain_level = sum(
        model_id is not None and chain_id is not None
        for model_id, chain_id in zip(
            model_ids,
            chain_ids,
            strict=True,
        )
    )

    return chain_level, table.num_rows - chain_level


def validate_quality_task_accounting(
    artifacts: Iterable[QualityTaskArtifacts],
    *,
    manifest_row_count: int,
    batch_size: int,
) -> tuple[QualityTaskValidation, ...]:
    """Validate task summaries against manifest partitions and shard contents."""

    results: list[QualityTaskValidation] = []

    for artifact in artifacts:
        task_id = artifact.task_id
        summary = artifact.summary

        counts = {
            field: _summary_nonnegative_int(
                summary,
                field,
                task_id=task_id,
            )
            for field in _SUMMARY_COUNT_FIELDS
        }

        expected_source_count = _expected_task_source_count(
            manifest_row_count,
            batch_size,
            task_id,
        )

        if (
            counts["input_source_object_count"]
            != expected_source_count
        ):
            raise QualityMergeError(
                f"Task {task_id} input source count mismatch: "
                f"expected {expected_source_count}, "
                f"summary reports "
                f"{counts['input_source_object_count']}"
            )

        shard_row_counts = {
            shard_name: _parquet_row_count(path)
            for shard_name, path in artifact.shard_paths.items()
        }

        summary_shard_fields = {
            "accepted": "accepted_chain_count",
            "rejected": "rejected_chain_count",
            "non_candidates": "non_candidate_chain_count",
            "dirty_residues": "dirty_residue_count",
            "errors": "processing_error_count",
        }

        for shard_name, summary_field in summary_shard_fields.items():
            observed = shard_row_counts[shard_name]
            reported = counts[summary_field]

            if observed != reported:
                raise QualityMergeError(
                    f"Task {task_id} {shard_name} row-count mismatch: "
                    f"shard has {observed}, summary reports {reported}"
                )

        accepted_trimmed = _accepted_trimmed_count(
            artifact.shard_paths["accepted"]
        )

        if (
            accepted_trimmed
            != counts["accepted_trimmed_chain_count"]
        ):
            raise QualityMergeError(
                f"Task {task_id} accepted-trimmed count mismatch"
            )

        chain_errors, source_errors = (
            _processing_error_level_counts(
                artifact.shard_paths["errors"]
            )
        )

        if (
            chain_errors
            != counts["chain_level_processing_error_count"]
        ):
            raise QualityMergeError(
                f"Task {task_id} chain-level processing-error "
                "count mismatch"
            )

        if (
            source_errors
            != counts["source_entry_processing_error_count"]
        ):
            raise QualityMergeError(
                f"Task {task_id} source-entry processing-error "
                "count mismatch"
            )

        if (
            counts["input_source_object_count"]
            != counts["successful_source_object_count"]
            + counts["failed_source_object_count"]
        ):
            raise QualityMergeError(
                f"Task {task_id} source-object accounting mismatch"
            )

        if (
            counts["selected_silver_chain_count"]
            != counts["accepted_chain_count"]
            + counts["rejected_chain_count"]
            + counts["non_candidate_chain_count"]
            + counts["chain_level_processing_error_count"]
        ):
            raise QualityMergeError(
                f"Task {task_id} selected-chain accounting mismatch"
            )

        calculated_total_gold = (
            counts["accepted_chain_count"]
            + counts["rejected_chain_count"]
            + counts["non_candidate_chain_count"]
            + counts["dirty_residue_count"]
            + counts["processing_error_count"]
        )

        if counts["total_gold_record_count"] != calculated_total_gold:
            raise QualityMergeError(
                f"Task {task_id} total Gold record count mismatch"
            )

        results.append(
            QualityTaskValidation(
                task_id=task_id,
                expected_input_source_object_count=(
                    expected_source_count
                ),
                shard_row_counts=shard_row_counts,
            )
        )

    return tuple(results)


@dataclass(frozen=True)
class QualityGlobalValidation:
    """Global identity and provenance facts for a complete quality stage."""

    accepted_chain_count: int
    rejected_chain_count: int
    non_candidate_chain_count: int
    dirty_residue_count: int
    processing_error_count: int
    unique_outcome_chain_count: int


def _manifest_source_index(
    manifest: pa.Table,
) -> dict[str, tuple[str, str]]:
    required = {"pdb_id", "s3_key", "etag", "snapshot"}

    missing = required - set(manifest.column_names)
    if missing:
        raise QualityMergeError(
            "Manifest is missing provenance columns: "
            + ", ".join(sorted(missing))
        )

    index: dict[str, tuple[str, str]] = {}

    for row in manifest.select(
        ["pdb_id", "s3_key", "etag"]
    ).to_pylist():
        pdb_id = row["pdb_id"]
        source = (row["s3_key"], row["etag"])

        if pdb_id in index:
            raise QualityMergeError(
                f"Duplicate manifest PDB ID: {pdb_id}"
            )

        index[pdb_id] = source

    return index


def _validate_source_provenance_rows(
    table: pa.Table,
    *,
    manifest_sources: dict[str, tuple[str, str]],
    expected_snapshot: str,
    expected_cleaning_protocol: str,
    expected_pipeline_git_commit: str,
    table_name: str,
) -> None:
    columns = set(table.column_names)

    for row in table.to_pylist():
        pdb_id = row["pdb_id"]

        expected_source = manifest_sources.get(pdb_id)
        if expected_source is None:
            raise QualityMergeError(
                f"{table_name} references PDB ID absent from manifest: "
                f"{pdb_id}"
            )

        observed_source = (
            row["source_mmcif_key"],
            row["source_etag"],
        )

        if observed_source != expected_source:
            raise QualityMergeError(
                f"{table_name} source provenance mismatch for "
                f"{pdb_id}"
            )

        if row["snapshot"] != expected_snapshot:
            raise QualityMergeError(
                f"{table_name} snapshot mismatch for {pdb_id}"
            )

        if (
            "cleaning_protocol" in columns
            and row["cleaning_protocol"]
            != expected_cleaning_protocol
        ):
            raise QualityMergeError(
                f"{table_name} cleaning protocol mismatch for "
                f"{pdb_id}"
            )

        if (
            "pipeline_git_commit" in columns
            and row["pipeline_git_commit"]
            != expected_pipeline_git_commit
        ):
            raise QualityMergeError(
                f"{table_name} Git commit mismatch for {pdb_id}"
            )


def _chain_identity(
    row: dict[str, Any],
) -> tuple[str, int, str]:
    return (
        row["pdb_id"],
        row["model_id"],
        row["label_chain_id"],
    )


def _iter_parquet_tables(
    path: Path,
    *,
    columns: list[str],
    batch_size: int = 65536,
) -> Iterable[pa.Table]:
    """Stream selected Parquet columns without loading a whole shard."""

    try:
        parquet_file = pq.ParquetFile(path)

        for batch in parquet_file.iter_batches(
            batch_size=batch_size,
            columns=columns,
        ):
            yield pa.Table.from_batches([batch])

    except Exception as exc:
        raise QualityMergeError(
            f"Cannot stream Parquet data from {path}: {exc}"
        ) from exc


def _stream_outcome_shards(
    artifacts: tuple[QualityTaskArtifacts, ...],
    *,
    manifest_sources: dict[str, tuple[str, str]],
    expected_snapshot: str,
    expected_cleaning_protocol: str,
    expected_pipeline_git_commit: str,
) -> tuple[
    dict[str, int],
    dict[tuple[str, int, str], str],
]:
    """Validate accepted/rejected/non-candidate rows incrementally."""

    counts = {
        "accepted": 0,
        "rejected": 0,
        "non_candidates": 0,
    }

    seen: dict[tuple[str, int, str], str] = {}

    columns = [
        "snapshot",
        "pdb_id",
        "model_id",
        "label_chain_id",
        "source_mmcif_key",
        "source_etag",
        "cleaning_protocol",
        "pipeline_git_commit",
    ]

    # Scan one outcome family at a time so duplicate/conflicting
    # identities are detected across both tasks and categories.
    for table_name in (
        "accepted",
        "rejected",
        "non_candidates",
    ):
        for artifact in artifacts:
            path = artifact.shard_paths[table_name]

            for table in _iter_parquet_tables(
                path,
                columns=columns,
            ):
                counts[table_name] += table.num_rows

                _validate_source_provenance_rows(
                    table,
                    manifest_sources=manifest_sources,
                    expected_snapshot=expected_snapshot,
                    expected_cleaning_protocol=(
                        expected_cleaning_protocol
                    ),
                    expected_pipeline_git_commit=(
                        expected_pipeline_git_commit
                    ),
                    table_name=table_name,
                )

                identity_rows = table.select(
                    [
                        "pdb_id",
                        "model_id",
                        "label_chain_id",
                    ]
                ).to_pylist()

                for row in identity_rows:
                    identity = _chain_identity(row)

                    previous = seen.get(identity)

                    if previous is not None:
                        raise QualityMergeError(
                            "Duplicate or conflicting Gold chain "
                            f"identity {identity!r}: found in "
                            f"{previous} and {table_name}"
                        )

                    seen[identity] = table_name

    return counts, seen


def _stream_error_shards(
    artifacts: tuple[QualityTaskArtifacts, ...],
    *,
    manifest_sources: dict[str, tuple[str, str]],
    expected_snapshot: str,
    expected_cleaning_protocol: str,
    expected_pipeline_git_commit: str,
    outcome_identities: set[tuple[str, int, str]],
) -> int:
    """Validate processing-error rows incrementally."""

    count = 0
    seen_errors: set[tuple[str, int, str]] = set()

    columns = [
        "snapshot",
        "pdb_id",
        "model_id",
        "label_chain_id",
        "source_mmcif_key",
        "source_etag",
        "pipeline_git_commit",
    ]

    for artifact in artifacts:
        path = artifact.shard_paths["errors"]

        for table in _iter_parquet_tables(
            path,
            columns=columns,
        ):
            count += table.num_rows

            _validate_source_provenance_rows(
                table,
                manifest_sources=manifest_sources,
                expected_snapshot=expected_snapshot,
                expected_cleaning_protocol=(
                    expected_cleaning_protocol
                ),
                expected_pipeline_git_commit=(
                    expected_pipeline_git_commit
                ),
                table_name="errors",
            )

            for row in table.select(
                [
                    "pdb_id",
                    "model_id",
                    "label_chain_id",
                ]
            ).to_pylist():
                model_id = row["model_id"]
                chain_id = row["label_chain_id"]

                if model_id is None or chain_id is None:
                    continue

                identity = (
                    row["pdb_id"],
                    model_id,
                    chain_id,
                )

                if identity in seen_errors:
                    raise QualityMergeError(
                        "Duplicate chain-level processing error: "
                        f"{identity!r}"
                    )

                if identity in outcome_identities:
                    raise QualityMergeError(
                        "Chain has both a Gold outcome and a "
                        f"processing error: {identity!r}"
                    )

                seen_errors.add(identity)

    return count


def _stream_dirty_residue_shards(
    artifacts: tuple[QualityTaskArtifacts, ...],
    *,
    manifest_sources: dict[str, tuple[str, str]],
    expected_snapshot: str,
    expected_cleaning_protocol: str,
    expected_pipeline_git_commit: str,
) -> int:
    """Validate dirty-residue rows incrementally."""

    count = 0
    seen: set[tuple[str, int, str, int]] = set()

    columns = [
        "snapshot",
        "pdb_id",
        "model_id",
        "label_chain_id",
        "label_seq_id",
        "source_mmcif_key",
        "source_etag",
    ]

    for artifact in artifacts:
        path = artifact.shard_paths["dirty_residues"]

        for table in _iter_parquet_tables(
            path,
            columns=columns,
        ):
            count += table.num_rows

            _validate_source_provenance_rows(
                table,
                manifest_sources=manifest_sources,
                expected_snapshot=expected_snapshot,
                expected_cleaning_protocol=(
                    expected_cleaning_protocol
                ),
                expected_pipeline_git_commit=(
                    expected_pipeline_git_commit
                ),
                table_name="dirty_residues",
            )

            for row in table.select(
                [
                    "pdb_id",
                    "model_id",
                    "label_chain_id",
                    "label_seq_id",
                ]
            ).to_pylist():
                identity = (
                    row["pdb_id"],
                    row["model_id"],
                    row["label_chain_id"],
                    row["label_seq_id"],
                )

                if identity in seen:
                    raise QualityMergeError(
                        "Duplicate dirty-residue identity: "
                        f"{identity!r}"
                    )

                seen.add(identity)

    return count


def validate_quality_global_state(
    artifacts: Iterable[QualityTaskArtifacts],
    *,
    manifest: pa.Table,
    expected_snapshot: str,
    expected_cleaning_protocol: str,
    expected_pipeline_git_commit: str,
) -> QualityGlobalValidation:
    """Stream and validate global Gold identity and provenance state."""

    artifacts = tuple(
        sorted(
            artifacts,
            key=lambda artifact: artifact.task_id,
        )
    )

    try:
        validate_manifest_table(
            manifest,
            expected_snapshot=expected_snapshot,
        )
    except ManifestError as exc:
        raise QualityMergeError(
            f"Bronze manifest validation failed: {exc}"
        ) from exc

    manifest_sources = _manifest_source_index(manifest)

    outcome_counts, outcome_map = _stream_outcome_shards(
        artifacts,
        manifest_sources=manifest_sources,
        expected_snapshot=expected_snapshot,
        expected_cleaning_protocol=(
            expected_cleaning_protocol
        ),
        expected_pipeline_git_commit=(
            expected_pipeline_git_commit
        ),
    )

    outcome_identities = set(outcome_map)

    processing_error_count = _stream_error_shards(
        artifacts,
        manifest_sources=manifest_sources,
        expected_snapshot=expected_snapshot,
        expected_cleaning_protocol=(
            expected_cleaning_protocol
        ),
        expected_pipeline_git_commit=(
            expected_pipeline_git_commit
        ),
        outcome_identities=outcome_identities,
    )

    dirty_residue_count = _stream_dirty_residue_shards(
        artifacts,
        manifest_sources=manifest_sources,
        expected_snapshot=expected_snapshot,
        expected_cleaning_protocol=(
            expected_cleaning_protocol
        ),
        expected_pipeline_git_commit=(
            expected_pipeline_git_commit
        ),
    )

    return QualityGlobalValidation(
        accepted_chain_count=outcome_counts["accepted"],
        rejected_chain_count=outcome_counts["rejected"],
        non_candidate_chain_count=(
            outcome_counts["non_candidates"]
        ),
        dirty_residue_count=dirty_residue_count,
        processing_error_count=processing_error_count,
        unique_outcome_chain_count=len(outcome_identities),
    )


QUALITY_GLOBAL_SUMMARY_SCHEMA_NAME = "pdbclean_quality_global_summary"
QUALITY_GLOBAL_SUMMARY_SCHEMA_VERSION = "1.0"
QUALITY_SUCCESS_SCHEMA_NAME = "pdbclean_quality_success"
QUALITY_SUCCESS_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class QualityMergePublication:
    """Published outputs for one completely validated quality stage."""

    merged_paths: dict[str, Path]
    global_summary_path: Path
    success_path: Path
    global_summary: dict[str, Any]


def _write_json_atomic(
    payload: dict[str, Any],
    output_path: str | Path,
) -> Path:
    """Write deterministic JSON atomically."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    temporary = output.with_suffix(output.suffix + ".tmp")
    text = json.dumps(
        payload,
        sort_keys=True,
        indent=2,
        ensure_ascii=True,
    )

    temporary.write_text(
        text + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)

    return output


def _write_merged_parquet_atomic(
    shard_paths: Iterable[Path],
    *,
    schema: pa.Schema,
    output_path: str | Path,
) -> Path:
    """Stream task shards into one deterministic schema-enforced Parquet."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    temporary = output.with_suffix(output.suffix + ".tmp")

    writer: pq.ParquetWriter | None = None

    try:
        writer = pq.ParquetWriter(
            temporary,
            schema,
            compression="zstd",
            version="2.6",
        )

        for shard_path in shard_paths:
            table = pq.read_table(shard_path)
            canonical = table.select(schema.names).cast(schema)
            writer.write_table(canonical)

    except Exception:
        if writer is not None:
            writer.close()

        if temporary.exists():
            temporary.unlink()

        raise

    else:
        assert writer is not None
        writer.close()
        temporary.replace(output)

    return output


def _aggregate_summary_count(
    artifacts: Iterable[QualityTaskArtifacts],
    field: str,
) -> int:
    return sum(
        _summary_nonnegative_int(
            artifact.summary,
            field,
            task_id=artifact.task_id,
        )
        for artifact in artifacts
    )


def build_quality_global_summary(
    artifacts: Iterable[QualityTaskArtifacts],
    *,
    manifest_row_count: int,
    batch_size: int,
    snapshot: str,
    cleaning_protocol: str,
    pipeline_git_commit: str,
    global_validation: QualityGlobalValidation,
) -> dict[str, Any]:
    """Build deterministic global quality-stage accounting summary."""

    artifacts = tuple(artifacts)

    summary = {
        "summary_schema_name": QUALITY_GLOBAL_SUMMARY_SCHEMA_NAME,
        "summary_schema_version": QUALITY_GLOBAL_SUMMARY_SCHEMA_VERSION,
        "snapshot": snapshot,
        "cleaning_protocol": cleaning_protocol,
        "pipeline_git_commit": pipeline_git_commit,
        "manifest_row_count": manifest_row_count,
        "batch_size": batch_size,
        "task_count": len(artifacts),
        "input_source_object_count": _aggregate_summary_count(
            artifacts,
            "input_source_object_count",
        ),
        "successful_source_object_count": _aggregate_summary_count(
            artifacts,
            "successful_source_object_count",
        ),
        "failed_source_object_count": _aggregate_summary_count(
            artifacts,
            "failed_source_object_count",
        ),
        "parsed_silver_chain_count": _aggregate_summary_count(
            artifacts,
            "parsed_silver_chain_count",
        ),
        "selected_silver_chain_count": _aggregate_summary_count(
            artifacts,
            "selected_silver_chain_count",
        ),
        "candidate_entry_count": _aggregate_summary_count(
            artifacts,
            "candidate_entry_count",
        ),
        "candidate_chain_count": _aggregate_summary_count(
            artifacts,
            "candidate_chain_count",
        ),
        "accepted_chain_count": (
            global_validation.accepted_chain_count
        ),
        "rejected_chain_count": (
            global_validation.rejected_chain_count
        ),
        "non_candidate_chain_count": (
            global_validation.non_candidate_chain_count
        ),
        "dirty_residue_count": (
            global_validation.dirty_residue_count
        ),
        "processing_error_count": (
            global_validation.processing_error_count
        ),
        "unique_outcome_chain_count": (
            global_validation.unique_outcome_chain_count
        ),
        "source_object_accounting_valid": (
            _aggregate_summary_count(
                artifacts,
                "input_source_object_count",
            )
            == _aggregate_summary_count(
                artifacts,
                "successful_source_object_count",
            )
            + _aggregate_summary_count(
                artifacts,
                "failed_source_object_count",
            )
        ),
        "selected_chain_accounting_valid": (
            _aggregate_summary_count(
                artifacts,
                "selected_silver_chain_count",
            )
            == global_validation.accepted_chain_count
            + global_validation.rejected_chain_count
            + global_validation.non_candidate_chain_count
            + _aggregate_summary_count(
                artifacts,
                "chain_level_processing_error_count",
            )
        ),
    }

    if summary["input_source_object_count"] != manifest_row_count:
        raise QualityMergeError(
            "Global input-source accounting does not cover the "
            "complete Bronze manifest"
        )

    if summary["source_object_accounting_valid"] is not True:
        raise QualityMergeError(
            "Global source-object accounting failed"
        )

    if summary["selected_chain_accounting_valid"] is not True:
        raise QualityMergeError(
            "Global selected-chain accounting failed"
        )

    return summary


def publish_quality_merge(
    artifacts: Iterable[QualityTaskArtifacts],
    *,
    quality_root: str | Path,
    manifest_row_count: int,
    batch_size: int,
    snapshot: str,
    cleaning_protocol: str,
    pipeline_git_commit: str,
    global_validation: QualityGlobalValidation,
) -> QualityMergePublication:
    """Publish merged quality outputs and stage completion marker last."""

    artifacts = tuple(sorted(
        artifacts,
        key=lambda artifact: artifact.task_id,
    ))

    root = Path(quality_root)
    merged_root = root / "merged"

    # A previous completion marker must never survive an attempted
    # republish. Once publication begins, the stage is incomplete until
    # a fresh _SUCCESS is written strictly last.
    success_path = root / "_SUCCESS"
    if success_path.exists():
        success_path.unlink()

    merged_paths: dict[str, Path] = {}

    for shard_name, schema in QUALITY_SHARD_SCHEMAS.items():
        output = merged_root / f"{shard_name}.parquet"

        merged_paths[shard_name] = _write_merged_parquet_atomic(
            (
                artifact.shard_paths[shard_name]
                for artifact in artifacts
            ),
            schema=schema,
            output_path=output,
        )

    global_summary = build_quality_global_summary(
        artifacts,
        manifest_row_count=manifest_row_count,
        batch_size=batch_size,
        snapshot=snapshot,
        cleaning_protocol=cleaning_protocol,
        pipeline_git_commit=pipeline_git_commit,
        global_validation=global_validation,
    )

    global_summary_path = _write_json_atomic(
        global_summary,
        root / "global_summary.json",
    )

    success_payload = {
        "success_schema_name": QUALITY_SUCCESS_SCHEMA_NAME,
        "success_schema_version": QUALITY_SUCCESS_SCHEMA_VERSION,
        "snapshot": snapshot,
        "cleaning_protocol": cleaning_protocol,
        "pipeline_git_commit": pipeline_git_commit,
        "manifest_row_count": manifest_row_count,
        "batch_size": batch_size,
        "task_count": len(artifacts),
        "global_summary": "global_summary.json",
        "merged_directory": "merged",
    }

    # Stage completion is published strictly last. Its presence means that
    # discovery, task accounting, global identity/provenance validation,
    # merged Parquet publication, and global-summary publication succeeded.
    success_path = _write_json_atomic(
        success_payload,
        root / "_SUCCESS",
    )

    return QualityMergePublication(
        merged_paths=merged_paths,
        global_summary_path=global_summary_path,
        success_path=success_path,
        global_summary=global_summary,
    )


def merge_quality_stage(
    *,
    quality_root: str | Path,
    manifest: pa.Table,
    manifest_row_count: int,
    batch_size: int,
    snapshot: str,
    cleaning_protocol: str,
    pipeline_git_commit: str,
) -> QualityMergePublication:
    """Validate and publish one complete distributed quality stage.

    A pre-existing _SUCCESS marker is invalidated before any validation
    begins. A fresh marker is written only after every validation and
    publication step succeeds.
    """

    root = Path(quality_root)
    success_path = root / "_SUCCESS"

    if success_path.exists():
        success_path.unlink()

    expected_ids = expected_quality_task_ids(
        manifest_row_count,
        batch_size,
    )

    artifacts = discover_quality_task_artifacts(
        root,
        expected_task_ids=expected_ids,
        expected_snapshot=snapshot,
        expected_cleaning_protocol=cleaning_protocol,
        expected_pipeline_git_commit=pipeline_git_commit,
    )

    validate_quality_task_accounting(
        artifacts,
        manifest_row_count=manifest_row_count,
        batch_size=batch_size,
    )

    global_validation = validate_quality_global_state(
        artifacts,
        manifest=manifest,
        expected_snapshot=snapshot,
        expected_cleaning_protocol=cleaning_protocol,
        expected_pipeline_git_commit=pipeline_git_commit,
    )

    return publish_quality_merge(
        artifacts,
        quality_root=root,
        manifest_row_count=manifest_row_count,
        batch_size=batch_size,
        snapshot=snapshot,
        cleaning_protocol=cleaning_protocol,
        pipeline_git_commit=pipeline_git_commit,
        global_validation=global_validation,
    )
