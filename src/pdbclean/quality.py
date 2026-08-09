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


def q002_disordered_residue_ids(
    chain: ChainObservation,
    *,
    backbone_atoms: tuple[str, ...] = ("N", "CA", "C"),
) -> set[int]:
    """Return residues flagged by BRI v1.2.2 `disorder_check`.

    Protocol 3.2 cleaning receives only backbone N/CA/C atoms.

    BRI flags:
      1. duplicate rows sharing residue_id, chain_id, model_id and atom name;
      2. occupancy tokens that do not start with "1", except ".".

    Alternate-location labels are not an independent rejection criterion.
    """

    backbone = [
        atom
        for atom in chain.atoms
        if atom.atom_name in backbone_atoms
    ]

    disordered: set[int] = set()

    # BRI:
    # chain.duplicated(
    #     ["residue_id", "chain_id", "model_id", "atom"],
    #     keep=False,
    # )
    counts: dict[tuple[int | None, str, int, str], int] = {}

    for atom in backbone:
        key = (
            atom.label_seq_id,
            atom.label_chain_id,
            atom.model_id,
            atom.atom_name,
        )
        counts[key] = counts.get(key, 0) + 1

    duplicate_keys = {
        key
        for key, count in counts.items()
        if count > 1
    }

    for atom in backbone:
        key = (
            atom.label_seq_id,
            atom.label_chain_id,
            atom.model_id,
            atom.atom_name,
        )

        if key in duplicate_keys and atom.label_seq_id is not None:
            disordered.add(atom.label_seq_id)

        # Preserve BRI's string-level occupancy semantics exactly.
        token = atom.occupancy_raw

        if token is None:
            # Production parser observations preserve the raw token.
            # Treat an absent raw token as unverifiable rather than
            # silently reconstructing BRI semantics from the float.
            if atom.label_seq_id is not None:
                disordered.add(atom.label_seq_id)
            continue

        if token != "." and not token.startswith("1"):
            if atom.label_seq_id is not None:
                disordered.add(atom.label_seq_id)

    return disordered


def evaluate_q002_disorder(
    chain: ChainObservation,
    *,
    backbone_atoms: tuple[str, ...] = ("N", "CA", "C"),
) -> RuleResult:
    """Q002: reproduce BRI v1.2.2 Protocol 3.2 disorder detection."""

    disordered = sorted(
        q002_disordered_residue_ids(
            chain,
            backbone_atoms=backbone_atoms,
        )
    )

    if disordered:
        return RuleResult(
            rule_id="Q002",
            passed=False,
            reason=(
                "disordered_backbone_residues:"
                + ",".join(str(residue_id) for residue_id in disordered)
            ),
        )

    return RuleResult(
        rule_id="Q002",
        passed=True,
        reason="no_disordered_backbone_residues",
    )

def missing_internal_label_seq_ids(
    chain: ChainObservation,
) -> list[int]:
    """Return residue-index breaks using BRI v1.2.2 continuity semantics.

    BRI constructs the integer range from the minimum to maximum residue_id
    and removes the observed residue IDs plus the special residue ID 0.
    """

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
        and seq_id != 0
    ]


def evaluate_q003_residue_continuity(
    chain: ChainObservation,
) -> RuleResult:
    """Q003: reproduce BRI v1.2.2 residue continuity checking."""

    # BRI's residue_id column is integer-typed. A missing label_seq_id in
    # our parsed representation therefore violates the executable input
    # precondition rather than representing an allowed gap.
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

    if internal_gaps:
        return RuleResult(
            rule_id="Q003",
            passed=False,
            reason=(
                "internal_label_seq_id_gaps:"
                + ",".join(str(value) for value in internal_gaps)
            ),
        )

    return RuleResult(
        rule_id="Q003",
        passed=True,
        reason="observed_residues_consecutive",
    )


def q004_incomplete_residue_ids(
    chain: ChainObservation,
    *,
    required_atoms: tuple[str, ...] = ("N", "CA", "C"),
) -> set[int]:
    """Return residues failing BRI v1.2.2 `residue_completeness_check`.

    `integrated_chainwise_filter` receives only backbone feature rows
    (N, CA and C). BRI counts rows per residue and marks a residue
    incomplete whenever that count is not exactly three.

    In the integrated Protocol 3.2 pipeline, duplicate backbone rows have
    already been handled by Q002 before this check is reached.
    """

    backbone = [
        atom
        for atom in chain.atoms
        if atom.atom_name in required_atoms
    ]

    counts: dict[int, int] = {}

    for atom in backbone:
        if atom.label_seq_id is None:
            continue

        counts[atom.label_seq_id] = (
            counts.get(atom.label_seq_id, 0) + 1
        )

    return {
        residue_id
        for residue_id, count in counts.items()
        if count != 3
    }


def evaluate_q004_backbone_completeness(
    chain: ChainObservation,
    *,
    required_atoms: tuple[str, ...] = ("N", "CA", "C"),
) -> RuleResult:
    """Q004: reproduce BRI v1.2.2 backbone completeness checking."""

    if any(
        atom.label_seq_id is None
        for atom in chain.atoms
        if atom.atom_name in required_atoms
    ):
        return RuleResult(
            rule_id="Q004",
            passed=False,
            reason="missing_label_seq_id",
        )

    observed_residue_ids = {
        atom.label_seq_id
        for atom in chain.atoms
        if atom.atom_name in required_atoms
        and atom.label_seq_id is not None
    }

    if not observed_residue_ids:
        return RuleResult(
            rule_id="Q004",
            passed=False,
            reason="no_observed_backbone_residues",
        )

    incomplete = sorted(
        q004_incomplete_residue_ids(
            chain,
            required_atoms=required_atoms,
        )
    )

    if incomplete:
        return RuleResult(
            rule_id="Q004",
            passed=False,
            reason=(
                "incomplete_backbone_residues:"
                + ",".join(str(residue_id) for residue_id in incomplete)
            ),
        )

    return RuleResult(
        rule_id="Q004",
        passed=True,
        reason="backbone_residue_row_counts_complete",
    )


import math


@dataclass(frozen=True)
class BackboneClash:
    """One Protocol 3.2 backbone clash.

    residue_id follows the BRI v1.2.2 convention: for an inter-residue
    comparison with N(i+1), the clash is assigned to residue i.
    """

    residue_id: int
    atom1: str
    atom2: str
    distance_angstrom: float


def _index_q005_backbone_atoms(
    chain: ChainObservation,
    *,
    required_atoms: tuple[str, ...],
) -> dict[int, dict[str, object]]:
    """Index exactly one configured backbone atom per observed residue."""

    if "N" not in required_atoms:
        raise ValueError("Q005 requires N among the configured backbone atoms")

    if any(atom.label_seq_id is None for atom in chain.atoms):
        raise ValueError("Q005 requires label_seq_id for all atoms")

    residue_ids = sorted(
        {
            atom.label_seq_id
            for atom in chain.atoms
            if atom.label_seq_id is not None
        }
    )

    if not residue_ids:
        raise ValueError("Q005 requires at least one observed residue")

    indexed: dict[int, dict[str, list]] = {
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

        indexed[atom.label_seq_id][atom.atom_name].append(atom)

    result: dict[int, dict[str, object]] = {}

    for residue_id in residue_ids:
        result[residue_id] = {}

        for atom_name in required_atoms:
            atoms = indexed[residue_id][atom_name]

            if len(atoms) != 1:
                raise ValueError(
                    "Q005 requires exactly one "
                    f"{atom_name} atom for residue {residue_id}"
                )

            result[residue_id][atom_name] = atoms[0]

    return result


def _atom_distance(left: object, right: object) -> float:
    """Euclidean distance between two parsed atom observations."""

    return math.dist(
        (left.x, left.y, left.z),
        (right.x, right.y, right.z),
    )


def protocol32_backbone_distances(
    chain: ChainObservation,
    *,
    required_atoms: tuple[str, ...],
) -> tuple[BackboneClash, ...]:
    """Return all distances tested by BRI v1.2.2 `clash_check`.

    For N, CA, C this reproduces the six `min_test_set` comparisons:

      N(i)-CA(i)
      N(i)-C(i)
      N(i)-N(i+1)
      CA(i)-C(i)
      CA(i)-N(i+1)
      C(i)-N(i+1)

    The returned objects represent tested distances, not necessarily clashes.
    """

    indexed = _index_q005_backbone_atoms(
        chain,
        required_atoms=required_atoms,
    )

    residue_ids = sorted(indexed)
    distances: list[BackboneClash] = []

    # Same-residue pairs. For N, CA, C these are:
    # N-CA, N-C, CA-C.
    for residue_id in residue_ids:
        residue_atoms = indexed[residue_id]

        for left_index, atom1 in enumerate(required_atoms):
            for atom2 in required_atoms[left_index + 1:]:
                distances.append(
                    BackboneClash(
                        residue_id=residue_id,
                        atom1=atom1,
                        atom2=atom2,
                        distance_angstrom=_atom_distance(
                            residue_atoms[atom1],
                            residue_atoms[atom2],
                        ),
                    )
                )

    # Inter-residue pairs against N(i+1). BRI assigns these to residue i.
    for residue_id, next_residue_id in zip(
        residue_ids,
        residue_ids[1:],
    ):
        if next_residue_id != residue_id + 1:
            raise ValueError(
                "Q005 requires consecutive residue ids; "
                f"found {residue_id} followed by {next_residue_id}"
            )

        next_n = indexed[next_residue_id]["N"]

        for atom1 in required_atoms:
            distances.append(
                BackboneClash(
                    residue_id=residue_id,
                    atom1=atom1,
                    atom2="N+1",
                    distance_angstrom=_atom_distance(
                        indexed[residue_id][atom1],
                        next_n,
                    ),
                )
            )

    return tuple(distances)


def minimum_protocol32_backbone_distance(
    chain: ChainObservation,
    *,
    required_atoms: tuple[str, ...],
) -> float | None:
    """Minimum distance over every Protocol 3.2 clash comparison."""

    distances = protocol32_backbone_distances(
        chain,
        required_atoms=required_atoms,
    )

    if not distances:
        return None

    return min(item.distance_angstrom for item in distances)


def find_q005_backbone_clashes(
    chain: ChainObservation,
    *,
    required_atoms: tuple[str, ...],
    minimum_distance_angstrom: float,
) -> tuple[BackboneClash, ...]:
    """Find Protocol 3.2 backbone clashes strictly below the threshold."""

    return tuple(
        item
        for item in protocol32_backbone_distances(
            chain,
            required_atoms=required_atoms,
        )
        if item.distance_angstrom < minimum_distance_angstrom
    )


def evaluate_q005_backbone_distance(
    chain: ChainObservation,
    *,
    required_atoms: tuple[str, ...],
    minimum_distance_angstrom: float,
) -> RuleResult:
    """Q005: reproduce the BRI v1.2.2 Protocol 3.2 clash criterion."""

    try:
        distances = protocol32_backbone_distances(
            chain,
            required_atoms=required_atoms,
        )
    except ValueError as exc:
        return RuleResult(
            rule_id="Q005",
            passed=False,
            reason=f"q005_precondition_failed:{exc}",
        )

    if not distances:
        return RuleResult(
            rule_id="Q005",
            passed=False,
            reason="insufficient_backbone_atoms",
        )

    clashes = tuple(
        item
        for item in distances
        if item.distance_angstrom < minimum_distance_angstrom
    )

    if clashes:
        minimum_distance = min(
            item.distance_angstrom
            for item in clashes
        )

        return RuleResult(
            rule_id="Q005",
            passed=False,
            reason=(
                "backbone_clashes:"
                f"{len(clashes)}:"
                f"{minimum_distance:.12g}"
            ),
        )

    return RuleResult(
        rule_id="Q005",
        passed=True,
        reason="backbone_clash_threshold_met",
    )

