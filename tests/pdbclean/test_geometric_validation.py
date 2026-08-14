"""Tests for post-cleaning geometric validation."""

from __future__ import annotations

from dataclasses import replace
import math

import pytest

from pdbclean.geometric_validation import (
    GeometricValidationConfig,
    reconstruct_retained_backbone_chain,
    validate_post_cleaning_geometry,
)
from pdbclean.mmcif_parser import AtomObservation, ChainObservation


def _atom(
    residue_id: int,
    atom_name: str,
    xyz: tuple[float, float, float],
) -> AtomObservation:
    return AtomObservation(
        model_id=1,
        label_chain_id="A",
        auth_chain_id="A",
        entity_id="1",
        label_seq_id=residue_id,
        auth_seq_id=str(residue_id),
        residue_name="ALA",
        atom_name=atom_name,
        alt_id=None,
        occupancy=1.0,
        x=xyz[0],
        y=xyz[1],
        z=xyz[2],
        group_pdb="ATOM",
        occupancy_raw="1.00",
    )


def _chain(
    residues: list[
        tuple[
            int,
            tuple[float, float, float],
            tuple[float, float, float],
            tuple[float, float, float],
        ]
    ],
) -> ChainObservation:
    atoms = []

    for residue_id, n, ca, c in residues:
        atoms.extend(
            [
                _atom(residue_id, "N", n),
                _atom(residue_id, "CA", ca),
                _atom(residue_id, "C", c),
            ]
        )

    return ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        auth_chain_id="A",
        entity_id="1",
        polymer_type="polypeptide(L)",
        entry_has_polypeptide=True,
        atoms=atoms,
    )


def test_valid_cleaned_chain_passes_default_validation() -> None:
    chain = _chain(
        [
            (
                1,
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
            ),
            (
                2,
                (2.0, 1.0, 0.0),
                (2.0, 2.0, 0.0),
                (3.0, 2.0, 0.0),
            ),
        ]
    )

    result = validate_post_cleaning_geometry(chain)

    assert result.passed is True
    assert result.violations == ()
    assert result.minimum_observed_backbone_distance_angstrom is not None
    assert result.minimum_observed_backbone_distance_angstrom >= 0.010
    assert result.minimum_observed_triangle_angle_degrees is not None
    assert result.minimum_observed_triangle_angle_degrees >= 3.0
    assert result.minimum_observed_basis_h_norm_angstrom is not None
    assert result.minimum_observed_basis_h_norm_angstrom > 0.0


def test_exact_backbone_distance_threshold_is_allowed() -> None:
    chain = _chain(
        [
            (
                1,
                (0.0, 0.0, 0.0),
                (0.010, 0.0, 0.0),
                (0.010, 1.0, 0.0),
            ),
        ]
    )

    result = validate_post_cleaning_geometry(
        chain,
        config=GeometricValidationConfig(
            minimum_backbone_distance_angstrom=0.010,
            minimum_triangle_angle_degrees=0.0,
        ),
    )

    assert not [
        violation
        for violation in result.violations
        if violation.violation_type
        == "backbone_distance_below_minimum"
    ]


def test_backbone_distance_threshold_is_configurable() -> None:
    chain = _chain(
        [
            (
                1,
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
            ),
        ]
    )

    default_result = validate_post_cleaning_geometry(chain)

    stricter_result = validate_post_cleaning_geometry(
        chain,
        config=GeometricValidationConfig(
            minimum_backbone_distance_angstrom=1.1,
            minimum_triangle_angle_degrees=3.0,
        ),
    )

    assert default_result.passed is True

    assert any(
        violation.violation_type
        == "backbone_distance_below_minimum"
        for violation in stricter_result.violations
    )


def test_two_degree_triangle_fails_default_three_degree_threshold() -> None:
    angle = math.radians(2.0)

    chain = _chain(
        [
            (
                1,
                (1.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (math.cos(angle), math.sin(angle), 0.0),
            ),
        ]
    )

    result = validate_post_cleaning_geometry(chain)

    assert result.passed is False

    assert any(
        violation.violation_type
        == "triangle_angle_below_minimum"
        for violation in result.violations
    )


def test_triangle_angle_threshold_is_configurable() -> None:
    angle = math.radians(2.0)

    chain = _chain(
        [
            (
                1,
                (1.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (math.cos(angle), math.sin(angle), 0.0),
            ),
        ]
    )

    result = validate_post_cleaning_geometry(
        chain,
        config=GeometricValidationConfig(
            minimum_backbone_distance_angstrom=0.010,
            minimum_triangle_angle_degrees=1.0,
        ),
    )

    assert result.passed is True


def test_collinear_residue_is_explicitly_rejected_before_bri() -> None:
    chain = _chain(
        [
            (
                1,
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (2.0, 0.0, 0.0),
            ),
        ]
    )

    result = validate_post_cleaning_geometry(chain)

    assert result.passed is False

    types = {
        violation.violation_type
        for violation in result.violations
    }

    assert "triangle_angle_below_minimum" in types
    assert "definition_3_4_undefined_h" in types


def test_nonconsecutive_retained_residue_ids_are_recorded() -> None:
    chain = _chain(
        [
            (
                1,
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
            ),
            (
                3,
                (2.0, 1.0, 0.0),
                (2.0, 2.0, 0.0),
                (3.0, 2.0, 0.0),
            ),
        ]
    )

    result = validate_post_cleaning_geometry(chain)

    assert result.passed is False

    assert any(
        violation.violation_type
        == "nonconsecutive_residue_ids"
        for violation in result.violations
    )


def test_reconstruct_retained_backbone_chain_uses_exact_gold_ids() -> None:
    chain = _chain(
        [
            (
                1,
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
            ),
            (
                2,
                (2.0, 1.0, 0.0),
                (2.0, 2.0, 0.0),
                (3.0, 2.0, 0.0),
            ),
            (
                3,
                (4.0, 2.0, 0.0),
                (4.0, 3.0, 0.0),
                (5.0, 3.0, 0.0),
            ),
        ]
    )

    # These must be removed by the same Protocol 3.2 projection
    # semantics used during cleaning.
    chain.atoms.append(
        replace(
            _atom(2, "CA", (99.0, 99.0, 99.0)),
            atom_name="CB",
        )
    )
    chain.atoms.append(
        replace(
            _atom(2, "CA", (88.0, 88.0, 88.0)),
            group_pdb="HETATM",
        )
    )

    retained = reconstruct_retained_backbone_chain(
        chain,
        [2, 3],
    )

    observed_ids = []
    seen = set()

    for atom in retained.atoms:
        residue_id = atom.label_seq_id

        if residue_id not in seen:
            seen.add(residue_id)
            observed_ids.append(residue_id)

    assert observed_ids == [2, 3]
    assert len(retained.atoms) == 6

    assert {
        atom.atom_name
        for atom in retained.atoms
    } == {"N", "CA", "C"}

    assert all(
        atom.group_pdb == "ATOM"
        for atom in retained.atoms
    )


def test_reconstruct_retained_backbone_chain_rejects_missing_gold_id() -> None:
    chain = _chain(
        [
            (
                1,
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
            ),
            (
                2,
                (2.0, 1.0, 0.0),
                (2.0, 2.0, 0.0),
                (3.0, 2.0, 0.0),
            ),
        ]
    )

    with pytest.raises(
        ValueError,
        match="do not match Gold lineage",
    ):
        reconstruct_retained_backbone_chain(
            chain,
            [1, 2, 3],
        )


def test_reconstruct_retained_backbone_chain_requires_gold_order() -> None:
    chain = _chain(
        [
            (
                1,
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
            ),
            (
                2,
                (2.0, 1.0, 0.0),
                (2.0, 2.0, 0.0),
                (3.0, 2.0, 0.0),
            ),
        ]
    )

    with pytest.raises(
        ValueError,
        match="do not match Gold lineage",
    ):
        reconstruct_retained_backbone_chain(
            chain,
            [2, 1],
        )


def test_reconstruct_retained_backbone_chain_rejects_duplicate_gold_ids() -> None:
    chain = _chain(
        [
            (
                1,
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
            ),
        ]
    )

    with pytest.raises(
        ValueError,
        match="must not contain duplicates",
    ):
        reconstruct_retained_backbone_chain(
            chain,
            [1, 1],
        )
