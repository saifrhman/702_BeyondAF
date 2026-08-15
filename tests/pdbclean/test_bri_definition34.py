"""Tests for MATCH Definition 3.4 Backbone Rigid Invariant."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from bri.invariant import get_invariant as reference_get_invariant

from pdbclean.bri import BRI_COLUMNS, compute_bri
from pdbclean.mmcif_parser import AtomObservation, ChainObservation


def _atom(
    residue_id: int,
    atom_name: str,
    xyz: tuple[float, float, float],
    *,
    residue_name: str = "ALA",
) -> AtomObservation:
    return AtomObservation(
        model_id=1,
        label_chain_id="A",
        auth_chain_id="A",
        entity_id="1",
        label_seq_id=residue_id,
        auth_seq_id=str(residue_id),
        residue_name=residue_name,
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


def _reference_dataframe(
    chain: ChainObservation,
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "model_id": atom.model_id,
                "chain_id": atom.label_chain_id,
                "residue_id": atom.label_seq_id,
                "residue_label": atom.residue_name,
                "atom": atom.atom_name,
                "x": atom.x,
                "y": atom.y,
                "z": atom.z,
            }
            for atom in chain.atoms
            if atom.group_pdb == "ATOM"
            and atom.atom_name in {"N", "CA", "C"}
        ]
    )


def _reference_bri(
    chain: ChainObservation,
) -> np.ndarray:
    reference = reference_get_invariant(
        _reference_dataframe(chain).copy(deep=True)
    )

    return reference.loc[
        :,
        list(BRI_COLUMNS),
    ].to_numpy(dtype=np.float64)


def _example_chain() -> ChainObservation:
    return _chain(
        [
            (
                1,
                (0.1, 0.2, 0.3),
                (1.2, -0.1, 0.5),
                (1.6, 1.1, 0.8),
            ),
            (
                2,
                (2.4, 1.5, 1.2),
                (2.9, 2.4, 1.6),
                (4.1, 2.2, 2.0),
            ),
            (
                3,
                (4.6, 3.0, 2.5),
                (5.2, 3.7, 3.1),
                (6.0, 3.3, 4.0),
            ),
        ]
    )


def test_bri_matches_pinned_reference() -> None:
    chain = _example_chain()

    observed = compute_bri(chain)
    expected = _reference_bri(chain)

    np.testing.assert_array_equal(
        observed,
        expected,
    )


def test_bri_shape_and_dtype() -> None:
    observed = compute_bri(_example_chain())

    assert observed.shape == (3, 9)
    assert observed.dtype == np.float64


def test_first_row_has_definition34_zero_structure() -> None:
    observed = compute_bri(_example_chain())

    np.testing.assert_array_equal(
        observed[0, [1, 2, 3, 4, 5, 8]],
        np.zeros(6, dtype=np.float64),
    )


def test_single_residue_chain_is_valid() -> None:
    chain = _chain(
        [
            (
                1,
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
            )
        ]
    )

    observed = compute_bri(chain)
    expected = _reference_bri(chain)

    assert observed.shape == (1, 9)

    np.testing.assert_array_equal(
        observed,
        expected,
    )


def test_bri_is_translation_invariant() -> None:
    original = _example_chain()
    translation = np.asarray(
        [17.25, -9.5, 31.125],
        dtype=np.float64,
    )

    translated_atoms = [
        AtomObservation(
            **{
                **atom.__dict__,
                "x": atom.x + translation[0],
                "y": atom.y + translation[1],
                "z": atom.z + translation[2],
            }
        )
        for atom in original.atoms
    ]

    translated = ChainObservation(
        **{
            **original.__dict__,
            "atoms": translated_atoms,
        }
    )

    np.testing.assert_array_equal(
        compute_bri(original),
        compute_bri(translated),
    )


def test_bri_is_rotation_invariant() -> None:
    original = _example_chain()

    theta = np.deg2rad(37.0)

    rotation = np.asarray(
        [
            [np.cos(theta), -np.sin(theta), 0.0],
            [np.sin(theta), np.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )

    rotated_atoms = []

    for atom in original.atoms:
        coordinate = np.asarray(
            [atom.x, atom.y, atom.z],
            dtype=np.float64,
        )

        rotated = rotation @ coordinate

        rotated_atoms.append(
            AtomObservation(
                **{
                    **atom.__dict__,
                    "x": rotated[0],
                    "y": rotated[1],
                    "z": rotated[2],
                }
            )
        )

    rotated_chain = ChainObservation(
        **{
            **original.__dict__,
            "atoms": rotated_atoms,
        }
    )

    np.testing.assert_array_equal(
        compute_bri(original),
        compute_bri(rotated_chain),
    )


def test_bri_rejects_nonconsecutive_residues() -> None:
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

    with pytest.raises(
        ValueError,
        match="consecutive retained residue IDs",
    ):
        compute_bri(chain)


def test_bri_rejects_undefined_definition34_basis() -> None:
    chain = _chain(
        [
            (
                1,
                (1.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (2.0, 0.0, 0.0),
            )
        ]
    )

    with pytest.raises(
        ValueError,
        match=r"\|h\|",
    ):
        compute_bri(chain)
