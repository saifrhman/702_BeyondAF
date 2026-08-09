"""Materialize Protocol 3.2 cleaning outcomes into Gold records."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.schemas import (
    GOLD_ACCEPTED_CHAIN_SCHEMA,
    GOLD_DIRTY_RESIDUE_SCHEMA,
    GOLD_NON_CANDIDATE_CHAIN_SCHEMA,
    GOLD_PROCESSING_ERROR_SCHEMA,
    GOLD_REJECTED_CHAIN_SCHEMA,
)


@dataclass(frozen=True)
class GoldTables:
    """Schema-enforced Gold quality tables for one processing batch."""

    accepted_chains: pa.Table
    rejected_chains: pa.Table
    non_candidate_chains: pa.Table
    dirty_residues: pa.Table
    processing_errors: pa.Table


QUALITY_TASK_SUMMARY_SCHEMA_NAME = "pdbclean_quality_task_summary"
QUALITY_TASK_SUMMARY_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class QualityTaskContext:
    """Execution metadata and upstream counts for one quality task."""

    task_id: str
    snapshot: str
    cleaning_protocol: str
    pipeline_git_commit: str

    started_at_utc: str
    completed_at_utc: str
    runtime_seconds: float

    input_source_object_count: int
    successful_source_object_count: int
    failed_source_object_count: int

    parsed_silver_chain_count: int
    selected_silver_chain_count: int
    candidate_entry_count: int
    candidate_chain_count: int

    slurm_job_id: str | None = None
    slurm_array_task_id: str | None = None
    peak_memory_bytes: int | None = None


@dataclass(frozen=True)
class GoldProvenance:
    """Release and source lineage required by Gold quality outputs."""

    snapshot: str
    source_mmcif_key: str
    source_etag: str
    cleaning_protocol: str
    pipeline_git_commit: str


@dataclass(frozen=True)
class GoldChainRecords:
    """Gold records produced from one parsed Silver chain."""

    accepted_chain: dict[str, Any] | None = None
    rejected_chain: dict[str, Any] | None = None
    non_candidate_chain: dict[str, Any] | None = None
    dirty_residues: tuple[dict[str, Any], ...] = ()


def _chain_identity(chain: "ChainObservation") -> dict[str, Any]:
    """Canonical chain identity shared by Gold quality outputs."""

    return {
        "pdb_id": chain.pdb_id,
        "model_id": chain.model_id,
        "label_chain_id": chain.label_chain_id,
        "auth_chain_id": chain.auth_chain_id,
        "entity_id": chain.entity_id,
    }


def _residue_ids_in_observation_order(
    chain: "ChainObservation",
) -> list[int]:
    """Unique label_seq_id values in surviving backbone row order."""

    residue_ids: list[int] = []
    seen: set[int] = set()

    for atom in chain.atoms:
        residue_id = atom.label_seq_id

        if residue_id is None:
            raise ValueError(
                "Protocol 3.2 Gold chain contains an atom without "
                "label_seq_id"
            )

        if residue_id not in seen:
            seen.add(residue_id)
            residue_ids.append(residue_id)

    return residue_ids


def _retained_sequence(
    chain: "ChainObservation",
    residue_ids: list[int],
) -> str:
    """Build the retained sequence using pinned BRI residue mapping."""

    from bri.base.base_util import amino_acid_short

    sequence: list[str] = []

    for residue_id in residue_ids:
        deposited = next(
            atom.residue_name
            for atom in chain.atoms
            if atom.label_seq_id == residue_id
        )
        mapped = amino_acid_short.get(deposited)

        if mapped is None:
            raise ValueError(
                "Accepted Protocol 3.2 residue has no pinned BRI "
                f"mapping: {deposited!r} at label_seq_id "
                f"{residue_id}"
            )

        sequence.append(mapped)

    return "".join(sequence)


def _dirty_rule_ids(
    result: "ChainCleaningResult",
) -> list[str]:
    """Return dirty rule IDs once each, preserving cleaning order."""

    return list(
        dict.fromkeys(
            record.rule_id
            for record in result.dirty_residues
        )
    )


def _dirty_residue_records(
    chain: "ChainObservation",
    result: "ChainCleaningResult",
    provenance: GoldProvenance,
) -> tuple[dict[str, Any], ...]:
    """Materialize residue-level Protocol 3.2 decision lineage."""

    return tuple(
        {
            "snapshot": provenance.snapshot,
            **_chain_identity(chain),
            "label_seq_id": record.residue_id,
            "deposited_residue_name": record.deposited_residue_name,
            "mapped_residue_code": record.mapped_residue_code,
            "rule_id": record.rule_id,
            "dirty_type": record.dirty_type,
            "cleaning_stage": record.cleaning_stage,
            "details_json": record.details_json,
            "source_mmcif_key": provenance.source_mmcif_key,
            "source_etag": provenance.source_etag,
        }
        for record in result.dirty_residues
    )


def gold_records_to_tables(
    records: Iterable[GoldChainRecords],
    processing_errors: Iterable[dict[str, Any]] = (),
) -> GoldTables:
    """Convert Gold records into schema-enforced Arrow tables."""

    accepted_rows: list[dict[str, Any]] = []
    rejected_rows: list[dict[str, Any]] = []
    non_candidate_rows: list[dict[str, Any]] = []
    dirty_rows: list[dict[str, Any]] = []
    error_rows = list(processing_errors)

    for record in records:
        terminal_rows = (
            record.accepted_chain,
            record.rejected_chain,
            record.non_candidate_chain,
        )

        terminal_count = sum(
            row is not None
            for row in terminal_rows
        )

        if terminal_count != 1:
            raise ValueError(
                "Each GoldChainRecords object must contain exactly one "
                "terminal chain-level outcome"
            )

        if record.accepted_chain is not None:
            accepted_rows.append(record.accepted_chain)

        if record.rejected_chain is not None:
            rejected_rows.append(record.rejected_chain)

        if record.non_candidate_chain is not None:
            non_candidate_rows.append(record.non_candidate_chain)

        dirty_rows.extend(record.dirty_residues)

    return GoldTables(
        accepted_chains=pa.Table.from_pylist(
            accepted_rows,
            schema=GOLD_ACCEPTED_CHAIN_SCHEMA,
        ),
        rejected_chains=pa.Table.from_pylist(
            rejected_rows,
            schema=GOLD_REJECTED_CHAIN_SCHEMA,
        ),
        non_candidate_chains=pa.Table.from_pylist(
            non_candidate_rows,
            schema=GOLD_NON_CANDIDATE_CHAIN_SCHEMA,
        ),
        dirty_residues=pa.Table.from_pylist(
            dirty_rows,
            schema=GOLD_DIRTY_RESIDUE_SCHEMA,
        ),
        processing_errors=pa.Table.from_pylist(
            error_rows,
            schema=GOLD_PROCESSING_ERROR_SCHEMA,
        ),
    )


def _sorted_counts(values: Iterable[str]) -> dict[str, int]:
    """Return deterministic lexical-key counts."""

    counts = Counter(values)
    return {
        key: counts[key]
        for key in sorted(counts)
    }


def build_quality_task_summary(
    tables: GoldTables,
    context: QualityTaskContext,
) -> dict[str, Any]:
    """Build one deterministic quality-task observability summary."""

    accepted_count = tables.accepted_chains.num_rows
    rejected_count = tables.rejected_chains.num_rows
    non_candidate_count = tables.non_candidate_chains.num_rows
    dirty_count = tables.dirty_residues.num_rows
    error_count = tables.processing_errors.num_rows

    error_rows = tables.processing_errors.to_pylist()

    chain_level_error_count = sum(
        row["model_id"] is not None
        and row["label_chain_id"] is not None
        for row in error_rows
    )
    source_entry_error_count = (
        error_count - chain_level_error_count
    )

    accepted_trimmed_count = sum(
        bool(value)
        for value in tables.accepted_chains[
            "terminal_trimmed"
        ].to_pylist()
    )

    source_accounting_valid = (
        context.input_source_object_count
        == context.successful_source_object_count
        + context.failed_source_object_count
    )

    chain_accounting_valid = (
        context.selected_silver_chain_count
        == accepted_count
        + rejected_count
        + non_candidate_count
        + chain_level_error_count
    )

    rejected_rows = tables.rejected_chains.to_pylist()
    dirty_rows = tables.dirty_residues.to_pylist()

    total_gold_record_count = (
        accepted_count
        + rejected_count
        + non_candidate_count
        + dirty_count
        + error_count
    )

    return {
        "summary_schema_name": QUALITY_TASK_SUMMARY_SCHEMA_NAME,
        "summary_schema_version": QUALITY_TASK_SUMMARY_SCHEMA_VERSION,
        "task_id": context.task_id,
        "snapshot": context.snapshot,
        "cleaning_protocol": context.cleaning_protocol,
        "pipeline_git_commit": context.pipeline_git_commit,
        "started_at_utc": context.started_at_utc,
        "completed_at_utc": context.completed_at_utc,
        "runtime_seconds": context.runtime_seconds,
        "slurm_job_id": context.slurm_job_id,
        "slurm_array_task_id": context.slurm_array_task_id,
        "peak_memory_bytes": context.peak_memory_bytes,
        "input_source_object_count": context.input_source_object_count,
        "successful_source_object_count": (
            context.successful_source_object_count
        ),
        "failed_source_object_count": context.failed_source_object_count,
        "parsed_silver_chain_count": context.parsed_silver_chain_count,
        "selected_silver_chain_count": (
            context.selected_silver_chain_count
        ),
        "candidate_entry_count": context.candidate_entry_count,
        "candidate_chain_count": context.candidate_chain_count,
        "non_candidate_chain_count": non_candidate_count,
        "accepted_chain_count": accepted_count,
        "accepted_trimmed_chain_count": accepted_trimmed_count,
        "rejected_chain_count": rejected_count,
        "dirty_residue_count": dirty_count,
        "processing_error_count": error_count,
        "chain_level_processing_error_count": chain_level_error_count,
        "source_entry_processing_error_count": source_entry_error_count,
        "total_gold_record_count": total_gold_record_count,
        "rejected_by_terminal_reason": _sorted_counts(
            row["terminal_reason"]
            for row in rejected_rows
        ),
        "rejected_by_terminal_stage": _sorted_counts(
            row["terminal_stage"]
            for row in rejected_rows
        ),
        "dirty_residues_by_rule_id": _sorted_counts(
            row["rule_id"]
            for row in dirty_rows
        ),
        "dirty_residues_by_type": _sorted_counts(
            row["dirty_type"]
            for row in dirty_rows
        ),
        "processing_errors_by_stage": _sorted_counts(
            row["processing_stage"]
            for row in error_rows
        ),
        "processing_errors_by_type": _sorted_counts(
            row["error_type"]
            for row in error_rows
        ),
        "source_object_accounting_valid": source_accounting_valid,
        "selected_chain_accounting_valid": chain_accounting_valid,
    }


def write_quality_task_summary_atomic(
    summary: dict[str, Any],
    output_root: str | Path,
) -> Path:
    """Validate and atomically write one quality-task JSON summary."""

    if summary.get("summary_schema_name") != QUALITY_TASK_SUMMARY_SCHEMA_NAME:
        raise ValueError("Unexpected quality-task summary schema name")

    if (
        summary.get("summary_schema_version")
        != QUALITY_TASK_SUMMARY_SCHEMA_VERSION
    ):
        raise ValueError("Unexpected quality-task summary schema version")

    if summary.get("source_object_accounting_valid") is not True:
        raise ValueError(
            "Cannot publish quality-task summary: "
            "source-object accounting failed"
        )

    if summary.get("selected_chain_accounting_valid") is not True:
        raise ValueError(
            "Cannot publish quality-task summary: "
            "selected-chain accounting failed"
        )

    task_id = str(summary.get("task_id", ""))

    if not task_id:
        raise ValueError("Quality-task summary requires a task_id")

    if "/" in task_id or "\\" in task_id or task_id in {".", ".."}:
        raise ValueError(f"Unsafe quality-task task_id: {task_id!r}")

    output = (
        Path(output_root)
        / "summaries"
        / f"task_{task_id}.json"
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    temporary = output.with_suffix(output.suffix + ".tmp")

    payload = json.dumps(
        summary,
        sort_keys=True,
        indent=2,
        ensure_ascii=True,
    )

    temporary.write_text(
        payload + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)

    return output


def _write_gold_table_atomic(
    table: pa.Table,
    schema: pa.Schema,
    output_path: str | Path,
) -> Path:
    """Write one Gold table atomically using its explicit schema."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    canonical = table.select(schema.names).cast(schema)
    temporary = output.with_suffix(output.suffix + ".tmp")

    pq.write_table(
        canonical,
        temporary,
        compression="zstd",
        version="2.6",
    )

    temporary.replace(output)
    return output


def write_gold_quality_shards(
    tables: GoldTables,
    output_root: str | Path,
    task_id: str | int,
) -> dict[str, Path]:
    """Write one task's schema-enforced Gold quality shards atomically."""

    root = Path(output_root)
    filename = f"task_{task_id}.parquet"

    outputs = {
        "accepted": root / "accepted" / filename,
        "rejected": root / "rejected" / filename,
        "non_candidates": root / "non_candidates" / filename,
        "dirty_residues": root / "dirty_residues" / filename,
        "errors": root / "errors" / filename,
    }

    _write_gold_table_atomic(
        tables.accepted_chains,
        GOLD_ACCEPTED_CHAIN_SCHEMA,
        outputs["accepted"],
    )
    _write_gold_table_atomic(
        tables.rejected_chains,
        GOLD_REJECTED_CHAIN_SCHEMA,
        outputs["rejected"],
    )
    _write_gold_table_atomic(
        tables.non_candidate_chains,
        GOLD_NON_CANDIDATE_CHAIN_SCHEMA,
        outputs["non_candidates"],
    )
    _write_gold_table_atomic(
        tables.dirty_residues,
        GOLD_DIRTY_RESIDUE_SCHEMA,
        outputs["dirty_residues"],
    )
    _write_gold_table_atomic(
        tables.processing_errors,
        GOLD_PROCESSING_ERROR_SCHEMA,
        outputs["errors"],
    )

    return outputs


def materialize_gold_chain(
    chain: "ChainObservation",
    result: "ChainCleaningResult",
    provenance: GoldProvenance,
) -> GoldChainRecords:
    """Materialize one scientific cleaning outcome into Gold records."""

    if result.status == "non_candidate":
        return GoldChainRecords(
            non_candidate_chain={
                "snapshot": provenance.snapshot,
                **_chain_identity(chain),
                "terminal_status": result.status,
                "terminal_reason": result.reason,
                "terminal_stage": result.terminal_stage,
                "source_mmcif_key": provenance.source_mmcif_key,
                "source_etag": provenance.source_etag,
                "cleaning_protocol": provenance.cleaning_protocol,
                "pipeline_git_commit": provenance.pipeline_git_commit,
            }
        )

    if result.status == "rejected":
        return GoldChainRecords(
            rejected_chain={
                "snapshot": provenance.snapshot,
                **_chain_identity(chain),
                "terminal_status": result.status,
                "terminal_reason": result.reason,
                "terminal_stage": result.terminal_stage,
                "missing_label_seq_ids": list(
                    result.missing_label_seq_ids
                ),
                "dirty_residue_count": len(result.dirty_residues),
                "dirty_rule_ids": _dirty_rule_ids(result),
                "source_mmcif_key": provenance.source_mmcif_key,
                "source_etag": provenance.source_etag,
                "cleaning_protocol": provenance.cleaning_protocol,
                "pipeline_git_commit": provenance.pipeline_git_commit,
            },
            dirty_residues=_dirty_residue_records(
                chain,
                result,
                provenance,
            ),
        )

    if result.status == "accepted":
        retained = result.retained_chain

        if retained is None or not retained.atoms:
            raise ValueError(
                "Accepted Protocol 3.2 result has no retained chain"
            )

        from pdbclean.quality import protocol32_backbone_projection

        original = protocol32_backbone_projection(chain)

        original_ids = _residue_ids_in_observation_order(original)
        retained_ids = _residue_ids_in_observation_order(retained)

        if not original_ids or not retained_ids:
            raise ValueError(
                "Accepted Protocol 3.2 chain has an empty backbone "
                "residue set"
            )

        return GoldChainRecords(
            accepted_chain={
                "snapshot": provenance.snapshot,
                **_chain_identity(chain),
                "original_start_label_seq_id": min(original_ids),
                "original_end_label_seq_id": max(original_ids),
                "retained_start_label_seq_id": min(retained_ids),
                "retained_end_label_seq_id": max(retained_ids),
                "retained_residue_count": len(retained_ids),
                "retained_label_seq_ids": retained_ids,
                "retained_sequence": _retained_sequence(
                    retained,
                    retained_ids,
                ),
                "terminal_trimmed": retained_ids != original_ids,
                "dirty_residue_count": len(result.dirty_residues),
                "dirty_rule_ids": _dirty_rule_ids(result),
                "source_mmcif_key": provenance.source_mmcif_key,
                "source_etag": provenance.source_etag,
                "cleaning_protocol": provenance.cleaning_protocol,
                "pipeline_git_commit": provenance.pipeline_git_commit,
            },
            dirty_residues=_dirty_residue_records(
                chain,
                result,
                provenance,
            ),
        )

    raise NotImplementedError(
        f"Gold materialization not implemented for status {result.status!r}"
    )


from pdbclean.cleaning import ChainCleaningResult
from pdbclean.mmcif_parser import ChainObservation
