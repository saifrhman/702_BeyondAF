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
