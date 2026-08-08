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
