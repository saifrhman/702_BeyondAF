"""Scientific quality rules for parsed PDB chain observations."""

from __future__ import annotations

from dataclasses import dataclass

from pdbclean.mmcif_parser import ChainObservation


@dataclass(frozen=True)
class RuleResult:
    """Outcome of one quality rule for one structural chain."""

    rule_id: str
    passed: bool
    reason: str


def evaluate_q001_protein_polymer(
    chain: ChainObservation,
    *,
    allowed_polymer_types: set[str],
) -> RuleResult:
    """Q001: retain only explicitly allowed protein polymer types.

    Polymer identity comes from `_entity_poly.type` via the chain's
    `_atom_site.label_entity_id`. Non-polymer entities therefore fail
    rather than being inferred from residue or atom names.
    """

    polymer_type = chain.polymer_type

    if polymer_type is None:
        return RuleResult(
            rule_id="Q001",
            passed=False,
            reason="missing_or_nonpolymer_entity_poly_type",
        )

    if polymer_type not in allowed_polymer_types:
        return RuleResult(
            rule_id="Q001",
            passed=False,
            reason=f"disallowed_polymer_type:{polymer_type}",
        )

    return RuleResult(
        rule_id="Q001",
        passed=True,
        reason="allowed_protein_polymer_type",
    )


def evaluate_q002_occupancy(
    chain: ChainObservation,
    *,
    minimum_occupancy: float,
    reject_alternate_locations: bool,
) -> RuleResult:
    """Q002: require complete occupancy for all polymer-chain atoms.

    The caller applies this rule after Q001 has established that the
    observation represents an allowed protein polymer. Missing occupancy
    is treated as unverifiable and therefore fails the rule.
    """

    missing_occupancy_count = sum(
        atom.occupancy is None
        for atom in chain.atoms
    )

    if missing_occupancy_count:
        return RuleResult(
            rule_id="Q002",
            passed=False,
            reason=(
                "missing_occupancy:"
                f"{missing_occupancy_count}"
            ),
        )

    below_minimum_count = sum(
        atom.occupancy is not None
        and atom.occupancy < minimum_occupancy
        for atom in chain.atoms
    )

    if below_minimum_count:
        return RuleResult(
            rule_id="Q002",
            passed=False,
            reason=(
                "occupancy_below_minimum:"
                f"{below_minimum_count}"
            ),
        )

    if reject_alternate_locations:
        alternate_location_count = sum(
            atom.alt_id is not None
            for atom in chain.atoms
        )

        if alternate_location_count:
            return RuleResult(
                rule_id="Q002",
                passed=False,
                reason=(
                    "alternate_locations_present:"
                    f"{alternate_location_count}"
                ),
            )

    return RuleResult(
        rule_id="Q002",
        passed=True,
        reason="occupancy_and_altloc_requirements_met",
    )


def missing_internal_label_seq_ids(
    chain: ChainObservation,
) -> list[int]:
    """Return absent label_seq_ids between observed minimum and maximum."""

    observed = sorted(
        {
            atom.label_seq_id
            for atom in chain.atoms
            if atom.label_seq_id is not None
        }
    )

    if len(observed) < 2:
        return []

    observed_set = set(observed)

    return [
        seq_id
        for seq_id in range(observed[0], observed[-1] + 1)
        if seq_id not in observed_set
    ]


def evaluate_q003_residue_continuity(
    chain: ChainObservation,
    *,
    allow_internal_gaps: bool,
) -> RuleResult:
    """Q003: require consecutive observed residues by label_seq_id.

    Missing residues before the first observed residue or after the last
    observed residue are terminal truncations and are not internal gaps.
    """

    missing_identifier_count = sum(
        atom.label_seq_id is None
        for atom in chain.atoms
    )

    if missing_identifier_count:
        return RuleResult(
            rule_id="Q003",
            passed=False,
            reason=(
                "missing_label_seq_id:"
                f"{missing_identifier_count}"
            ),
        )

    observed_ids = {
        atom.label_seq_id
        for atom in chain.atoms
        if atom.label_seq_id is not None
    }

    if not observed_ids:
        return RuleResult(
            rule_id="Q003",
            passed=False,
            reason="no_observed_label_seq_ids",
        )

    internal_gaps = missing_internal_label_seq_ids(chain)

    if internal_gaps and not allow_internal_gaps:
        return RuleResult(
            rule_id="Q003",
            passed=False,
            reason=(
                "internal_label_seq_id_gaps:"
                + ",".join(str(value) for value in internal_gaps)
            ),
        )

    if internal_gaps:
        return RuleResult(
            rule_id="Q003",
            passed=True,
            reason="internal_gaps_allowed",
        )

    return RuleResult(
        rule_id="Q003",
        passed=True,
        reason="observed_residues_consecutive",
    )


@dataclass(frozen=True)
class BackboneIssueSummary:
    """Missing and duplicate required backbone atoms for one chain."""

    missing_atom_count: int
    duplicate_atom_count: int
    missing_atoms: tuple[str, ...]
    duplicate_atoms: tuple[str, ...]


def summarize_backbone_atom_issues(
    chain: ChainObservation,
    *,
    required_atoms: tuple[str, ...],
) -> BackboneIssueSummary:
    """Summarize required backbone atom multiplicity by label_seq_id."""

    residue_ids = sorted(
        {
            atom.label_seq_id
            for atom in chain.atoms
            if atom.label_seq_id is not None
        }
    )

    counts: dict[int, dict[str, int]] = {
        residue_id: {
            atom_name: 0
            for atom_name in required_atoms
        }
        for residue_id in residue_ids
    }

    for atom in chain.atoms:
        if atom.label_seq_id is None:
            continue

        if atom.atom_name not in required_atoms:
            continue

        counts[atom.label_seq_id][atom.atom_name] += 1

    missing_atoms: list[str] = []
    duplicate_atoms: list[str] = []
    duplicate_atom_count = 0

    for residue_id in residue_ids:
        for atom_name in required_atoms:
            count = counts[residue_id][atom_name]

            if count == 0:
                missing_atoms.append(
                    f"{residue_id}:{atom_name}"
                )
            elif count > 1:
                duplicate_atoms.append(
                    f"{residue_id}:{atom_name}:{count}"
                )
                duplicate_atom_count += count - 1

    return BackboneIssueSummary(
        missing_atom_count=len(missing_atoms),
        duplicate_atom_count=duplicate_atom_count,
        missing_atoms=tuple(missing_atoms),
        duplicate_atoms=tuple(duplicate_atoms),
    )


def evaluate_q004_backbone_atoms(
    chain: ChainObservation,
    *,
    required_atoms: tuple[str, ...],
    require_exactly_one: bool,
) -> RuleResult:
    """Q004: require the configured backbone atoms in every residue."""

    if any(
        atom.label_seq_id is None
        for atom in chain.atoms
    ):
        return RuleResult(
            rule_id="Q004",
            passed=False,
            reason="missing_label_seq_id",
        )

    observed_residue_ids = {
        atom.label_seq_id
        for atom in chain.atoms
        if atom.label_seq_id is not None
    }

    if not observed_residue_ids:
        return RuleResult(
            rule_id="Q004",
            passed=False,
            reason="no_observed_residues",
        )

    issues = summarize_backbone_atom_issues(
        chain,
        required_atoms=required_atoms,
    )

    if issues.missing_atom_count:
        return RuleResult(
            rule_id="Q004",
            passed=False,
            reason=(
                "missing_backbone_atoms:"
                f"{issues.missing_atom_count}:"
                + ",".join(issues.missing_atoms)
            ),
        )

    if (
        require_exactly_one
        and issues.duplicate_atom_count
    ):
        return RuleResult(
            rule_id="Q004",
            passed=False,
            reason=(
                "duplicate_backbone_atoms:"
                f"{issues.duplicate_atom_count}:"
                + ",".join(issues.duplicate_atoms)
            ),
        )

    return RuleResult(
        rule_id="Q004",
        passed=True,
        reason="required_backbone_atoms_present_once",
    )


import math


def ordered_backbone_atoms(
    chain: ChainObservation,
    *,
    required_atoms: tuple[str, ...],
) -> list:
    """Return backbone atoms ordered by residue then configured atom order.

    Q005 assumes Q004 semantics: every observed residue must contain
    exactly one occurrence of every configured backbone atom.
    """

    residue_ids = sorted(
        {
            atom.label_seq_id
            for atom in chain.atoms
            if atom.label_seq_id is not None
        }
    )

    by_residue: dict[int, dict[str, list]] = {
        residue_id: {
            atom_name: []
            for atom_name in required_atoms
        }
        for residue_id in residue_ids
    }

    for atom in chain.atoms:
        if atom.label_seq_id is None:
            continue

        if atom.atom_name not in required_atoms:
            continue

        by_residue[atom.label_seq_id][atom.atom_name].append(atom)

    ordered = []

    for residue_id in residue_ids:
        for atom_name in required_atoms:
            atoms = by_residue[residue_id][atom_name]

            if len(atoms) != 1:
                raise ValueError(
                    "Q005 requires exactly one "
                    f"{atom_name} atom for residue {residue_id}"
                )

            ordered.append(atoms[0])

    return ordered


def minimum_consecutive_backbone_distance(
    chain: ChainObservation,
    *,
    required_atoms: tuple[str, ...],
) -> float | None:
    """Return the minimum distance between consecutive ordered backbone atoms."""

    atoms = ordered_backbone_atoms(
        chain,
        required_atoms=required_atoms,
    )

    if len(atoms) < 2:
        return None

    return min(
        math.dist(
            (left.x, left.y, left.z),
            (right.x, right.y, right.z),
        )
        for left, right in zip(atoms, atoms[1:])
    )


def evaluate_q005_backbone_distance(
    chain: ChainObservation,
    *,
    required_atoms: tuple[str, ...],
    minimum_distance_angstrom: float,
) -> RuleResult:
    """Q005: reject implausibly coincident consecutive backbone atoms."""

    try:
        minimum_distance = minimum_consecutive_backbone_distance(
            chain,
            required_atoms=required_atoms,
        )
    except ValueError:
        return RuleResult(
            rule_id="Q005",
            passed=False,
            reason="backbone_not_exactly_one_per_residue",
        )

    if minimum_distance is None:
        return RuleResult(
            rule_id="Q005",
            passed=False,
            reason="insufficient_backbone_atoms",
        )

    if minimum_distance < minimum_distance_angstrom:
        return RuleResult(
            rule_id="Q005",
            passed=False,
            reason=(
                "consecutive_backbone_distance_below_minimum:"
                f"{minimum_distance:.12g}"
            ),
        )

    return RuleResult(
        rule_id="Q005",
        passed=True,
        reason="backbone_distances_meet_minimum",
    )
