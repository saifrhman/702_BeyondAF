"""Materialize Protocol 3.2 cleaning outcomes into Gold records."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

import pyarrow as pa

from pdbclean.schemas import (
    GOLD_ACCEPTED_CHAIN_SCHEMA,
    GOLD_DIRTY_RESIDUE_SCHEMA,
    GOLD_NON_CANDIDATE_CHAIN_SCHEMA,
    GOLD_REJECTED_CHAIN_SCHEMA,
)


@dataclass(frozen=True)
class GoldTables:
    """Schema-enforced Gold quality tables for one processing batch."""

    accepted_chains: pa.Table
    rejected_chains: pa.Table
    non_candidate_chains: pa.Table
    dirty_residues: pa.Table


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
) -> GoldTables:
    """Convert per-chain Gold records into schema-enforced Arrow tables."""

    accepted_rows: list[dict[str, Any]] = []
    rejected_rows: list[dict[str, Any]] = []
    non_candidate_rows: list[dict[str, Any]] = []
    dirty_rows: list[dict[str, Any]] = []

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
    )


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
