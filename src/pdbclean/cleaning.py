"""Stateful Protocol 3.2 chain-cleaning orchestration.

Individual scientific rule semantics live in :mod:`pdbclean.quality`.
This module composes those rules in the stateful order used by pinned
BRI v1.2.2.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import json
from typing import Literal

from pdbclean.mmcif_parser import ChainObservation
from pdbclean.quality import (
    evaluate_q001_entry_protein,
    evaluate_q003_residue_continuity,
    find_q005_backbone_clashes,
    missing_internal_label_seq_ids,
    protocol32_backbone_projection,
    q002_disordered_residue_ids,
    q004_incomplete_residue_ids,
    q006_nonstandard_residue_ids,
)


ChainCleaningStatus = Literal[
    "non_candidate",
    "accepted",
    "rejected",
]


@dataclass(frozen=True)
class DirtyResidueRecord:
    """One residue-level Protocol 3.2 cleaning decision."""

    residue_id: int
    deposited_residue_name: str
    mapped_residue_code: str | None
    rule_id: str
    dirty_type: str
    cleaning_stage: str
    details_json: str = "{}"


@dataclass(frozen=True)
class ChainCleaningResult:
    """Terminal scientific outcome for one parsed Silver chain."""

    status: ChainCleaningStatus
    reason: str
    terminal_stage: str

    retained_chain: ChainObservation | None = None
    dirty_residues: tuple[DirtyResidueRecord, ...] = ()
    missing_label_seq_ids: tuple[int, ...] = ()


def _remove_residue_ids(
    chain: ChainObservation,
    residue_ids: set[int],
) -> ChainObservation:
    """Return a working-chain copy without the selected residue IDs."""

    return replace(
        chain,
        atoms=[
            atom
            for atom in chain.atoms
            if atom.label_seq_id not in residue_ids
        ],
    )


def _deposited_residue_name(
    chain: ChainObservation,
    residue_id: int,
) -> str:
    """Return the deposited residue name for one observed residue."""

    for atom in chain.atoms:
        if atom.label_seq_id == residue_id:
            return atom.residue_name

    raise ValueError(
        f"Residue {residue_id} is absent from chain {chain.canonical_key}"
    )


def _json_details(payload: dict[str, object]) -> str:
    """Serialize Gold rule evidence deterministically."""

    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    )


def _q002_details_json(
    chain: ChainObservation,
    residue_id: int,
) -> str:
    """Capture BRI Q002 evidence before the residue is removed."""

    residue_atoms = [
        atom
        for atom in chain.atoms
        if atom.label_seq_id == residue_id
    ]

    counts: dict[tuple[int | None, str, int, str], int] = {}

    for atom in chain.atoms:
        key = (
            atom.label_seq_id,
            atom.label_chain_id,
            atom.model_id,
            atom.atom_name,
        )
        counts[key] = counts.get(key, 0) + 1

    duplicate_atom_names = sorted(
        {
            atom.atom_name
            for atom in residue_atoms
            if counts[
                (
                    atom.label_seq_id,
                    atom.label_chain_id,
                    atom.model_id,
                    atom.atom_name,
                )
            ] > 1
        }
    )

    bad_occupancy_atoms = [
        {
            "atom_name": atom.atom_name,
            "occupancy_raw": atom.occupancy_raw,
            "alt_id": atom.alt_id,
        }
        for atom in residue_atoms
        if atom.occupancy_raw is None
        or (
            atom.occupancy_raw != "."
            and not atom.occupancy_raw.startswith("1")
        )
    ]

    return _json_details(
        {
            "duplicate_backbone_atom_names": duplicate_atom_names,
            "bad_occupancy_atoms": bad_occupancy_atoms,
        }
    )


def _q006_details_json(
    chain: ChainObservation,
    residue_id: int,
) -> str:
    """Capture pinned BRI Q006 residue-mapping evidence."""

    return _json_details(
        {
            "deposited_residue_name": _deposited_residue_name(
                chain,
                residue_id,
            ),
            "mapped_residue_code": _mapped_residue_code(
                chain,
                residue_id,
            ),
            "mapping_source": "bri_v1.2.2_amino_acid_short",
            "accepted_code_set": "canonical_20_one_letter",
        }
    )


def _q004_details_json(
    chain: ChainObservation,
    residue_id: int,
) -> str:
    """Capture BRI Q004 backbone-row completeness evidence."""

    residue_atoms = [
        atom
        for atom in chain.atoms
        if atom.label_seq_id == residue_id
    ]

    return _json_details(
        {
            "backbone_row_count": len(residue_atoms),
            "backbone_atom_names": sorted(
                atom.atom_name
                for atom in residue_atoms
            ),
            "required_backbone_atoms": ["N", "CA", "C"],
        }
    )


def _mapped_residue_code(
    chain: ChainObservation,
    residue_id: int,
) -> str | None:
    """Return the pinned BRI CCD-derived one-letter residue mapping."""

    from bri.base.base_util import amino_acid_short

    deposited = _deposited_residue_name(chain, residue_id)
    return amino_acid_short.get(deposited)


def _q005_details_json(
    residue_id: int,
    clashes: tuple,
) -> str:
    """Capture strict BRI Q005 clash evidence before residue removal."""

    residue_clashes = [
        {
            "atom1": clash.atom1,
            "atom2": clash.atom2,
            "distance_angstrom": clash.distance_angstrom,
        }
        for clash in clashes
        if clash.residue_id == residue_id
    ]

    return _json_details(
        {
            "threshold_angstrom": 0.01,
            "comparison": "distance < threshold",
            "clashes": residue_clashes,
        }
    )


def is_protocol32_candidate(chain: ChainObservation) -> bool:
    """Whether BRI v1.2.2 would send this chain to chainwise cleaning.

    ``on_entry(HETATM=False)`` skips chains whose ATOM N/CA/C feature
    projection is empty.
    """

    return bool(protocol32_backbone_projection(chain).atoms)


def clean_protocol32_chain(
    chain: ChainObservation,
) -> ChainCleaningResult:
    """Apply pinned BRI v1.2.2 Protocol 3.2 statefully."""

    working = protocol32_backbone_projection(chain)

    if not working.atoms:
        return ChainCleaningResult(
            status="non_candidate",
            reason="empty_protocol32_backbone_projection",
            terminal_stage="candidate_selection",
        )

    q001 = evaluate_q001_entry_protein(chain)

    if not q001.passed:
        return ChainCleaningResult(
            status="rejected",
            reason=q001.reason,
            terminal_stage="Q001",
        )

    # Q002: detect dirty disorder residues, record them, then remove the
    # complete residues from the current working chain.
    q002_ids = set(q002_disordered_residue_ids(working))

    dirty_residues = tuple(
        DirtyResidueRecord(
            residue_id=residue_id,
            deposited_residue_name=_deposited_residue_name(
                working,
                residue_id,
            ),
            mapped_residue_code=None,
            rule_id="Q002",
            dirty_type="disordered",
            cleaning_stage="Q002",
            details_json=_q002_details_json(
                working,
                residue_id,
            ),
        )
        for residue_id in sorted(q002_ids)
    )

    if q002_ids:
        working = _remove_residue_ids(working, q002_ids)

        if not working.atoms:
            return ChainCleaningResult(
                status="rejected",
                reason="all_working_residues_removed",
                terminal_stage="Q002",
                dirty_residues=dirty_residues,
            )

    # BRI checks continuity immediately after Q002 residue removal.
    q003 = evaluate_q003_residue_continuity(working)

    if not q003.passed:
        return ChainCleaningResult(
            status="rejected",
            reason=q003.reason,
            terminal_stage="Q003_after_Q002",
            dirty_residues=dirty_residues,
            missing_label_seq_ids=tuple(
                missing_internal_label_seq_ids(working)
            ),
        )

    # Q006 and Q004 are evaluated on the same current working chain.
    # BRI appends non-standard rows before incomplete rows and then
    # de-duplicates by residue identity, so Q006 wins when one residue
    # fails both checks.
    q006_ids = set(q006_nonstandard_residue_ids(working))
    q004_ids = set(q004_incomplete_residue_ids(working))

    q006_records = tuple(
        DirtyResidueRecord(
            residue_id=residue_id,
            deposited_residue_name=_deposited_residue_name(
                working,
                residue_id,
            ),
            mapped_residue_code=_mapped_residue_code(
                working,
                residue_id,
            ),
            rule_id="Q006",
            dirty_type="non-standard",
            cleaning_stage="Q006_Q004",
            details_json=_q006_details_json(
                working,
                residue_id,
            ),
        )
        for residue_id in sorted(q006_ids)
    )

    # Exact BRI dirty-record precedence: non-standard precedes incomplete.
    q004_only_ids = q004_ids - q006_ids

    q004_records = tuple(
        DirtyResidueRecord(
            residue_id=residue_id,
            deposited_residue_name=_deposited_residue_name(
                working,
                residue_id,
            ),
            mapped_residue_code=_mapped_residue_code(
                working,
                residue_id,
            ),
            rule_id="Q004",
            dirty_type="incomplete",
            cleaning_stage="Q006_Q004",
            details_json=_q004_details_json(
                working,
                residue_id,
            ),
        )
        for residue_id in sorted(q004_only_ids)
    )

    dirty_residues = (
        dirty_residues
        + q006_records
        + q004_records
    )

    q006_q004_ids = q006_ids | q004_ids

    if q006_q004_ids:
        working = _remove_residue_ids(
            working,
            q006_q004_ids,
        )

        if not working.atoms:
            return ChainCleaningResult(
                status="rejected",
                reason="all_working_residues_removed",
                terminal_stage="Q006_Q004",
                dirty_residues=dirty_residues,
            )

    # BRI immediately checks continuity after the combined Q006/Q004
    # residue-removal stage.
    q003 = evaluate_q003_residue_continuity(working)

    if not q003.passed:
        return ChainCleaningResult(
            status="rejected",
            reason=q003.reason,
            terminal_stage="Q003_after_Q006_Q004",
            dirty_residues=dirty_residues,
            missing_label_seq_ids=tuple(
                missing_internal_label_seq_ids(working)
            ),
        )

    # Q005: detect strict backbone clashes below 0.01 A.
    clashes = find_q005_backbone_clashes(
        working,
        required_atoms=("N", "CA", "C"),
        minimum_distance_angstrom=0.01,
    )

    # BRI ultimately de-duplicates dirty rows by residue identity.
    q005_ids = {
        clash.residue_id
        for clash in clashes
    }

    q005_records = tuple(
        DirtyResidueRecord(
            residue_id=residue_id,
            deposited_residue_name=_deposited_residue_name(
                working,
                residue_id,
            ),
            mapped_residue_code=_mapped_residue_code(
                working,
                residue_id,
            ),
            rule_id="Q005",
            dirty_type="clash",
            cleaning_stage="Q005",
            details_json=_q005_details_json(
                residue_id,
                clashes,
            ),
        )
        for residue_id in sorted(q005_ids)
    )

    dirty_residues = dirty_residues + q005_records

    if q005_ids:
        working = _remove_residue_ids(
            working,
            q005_ids,
        )

        if not working.atoms:
            return ChainCleaningResult(
                status="rejected",
                reason="all_working_residues_removed",
                terminal_stage="Q005",
                dirty_residues=dirty_residues,
            )

    # Final BRI continuity check after clash-residue removal.
    q003 = evaluate_q003_residue_continuity(working)

    if not q003.passed:
        return ChainCleaningResult(
            status="rejected",
            reason=q003.reason,
            terminal_stage="Q003_final",
            dirty_residues=dirty_residues,
            missing_label_seq_ids=tuple(
                missing_internal_label_seq_ids(working)
            ),
        )

    return ChainCleaningResult(
        status="accepted",
        reason="protocol32_clean",
        terminal_stage="completed",
        retained_chain=working,
        dirty_residues=dirty_residues,
    )

