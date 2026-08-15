"""Backbone Rigid Invariant (BRI) computation.

This module implements the complete nine-coordinate Backbone Rigid
Invariant of MATCH Definition 3.4 for an already-cleaned,
post-geometric-validation protein backbone.

It does not perform Protocol 3.2 cleaning and does not decide which
residues belong to a retained chain.
"""

from __future__ import annotations

import math

import numpy as np

from pdbclean.mmcif_parser import AtomObservation, ChainObservation
from pdbclean.quality import protocol32_backbone_projection


BRI_COLUMNS = (
    "x(N)",
    "y(N)",
    "z(N)",
    "x(A)",
    "y(A)",
    "z(A)",
    "x(C)",
    "y(C)",
    "z(C)",
)

BACKBONE_ATOMS = ("N", "CA", "C")


def _xyz(atom: AtomObservation) -> np.ndarray:
    """Return one atomic coordinate as float64."""

    coordinate = np.asarray(
        [atom.x, atom.y, atom.z],
        dtype=np.float64,
    )

    if not np.isfinite(coordinate).all():
        raise ValueError(
            "BRI computation requires finite N/CA/C coordinates"
        )

    return coordinate


def _index_backbone(
    chain: ChainObservation,
) -> list[
    tuple[
        int,
        np.ndarray,
        np.ndarray,
        np.ndarray,
    ]
]:
    """Return ordered residues containing exactly one N, CA and C."""

    backbone = protocol32_backbone_projection(
        chain,
        backbone_atoms=BACKBONE_ATOMS,
    )

    indexed: dict[
        int,
        dict[str, list[AtomObservation]],
    ] = {}

    for atom in backbone.atoms:
        if atom.label_seq_id is None:
            raise ValueError(
                "BRI computation requires label_seq_id"
            )

        residue_id = atom.label_seq_id

        indexed.setdefault(
            residue_id,
            {name: [] for name in BACKBONE_ATOMS},
        )

        indexed[residue_id][atom.atom_name].append(atom)

    if not indexed:
        raise ValueError(
            "BRI computation requires at least one residue"
        )

    residue_ids = sorted(indexed)

    for left, right in zip(
        residue_ids,
        residue_ids[1:],
    ):
        if right != left + 1:
            raise ValueError(
                "BRI computation requires consecutive retained "
                f"residue IDs; found {left} followed by {right}"
            )

    result = []

    for residue_id in residue_ids:
        atoms = indexed[residue_id]

        selected: dict[str, AtomObservation] = {}

        for atom_name in BACKBONE_ATOMS:
            observations = atoms[atom_name]

            if len(observations) != 1:
                raise ValueError(
                    "BRI computation requires exactly one "
                    f"{atom_name} atom for residue {residue_id}; "
                    f"found {len(observations)}"
                )

            selected[atom_name] = observations[0]

        result.append(
            (
                residue_id,
                _xyz(selected["N"]),
                _xyz(selected["CA"]),
                _xyz(selected["C"]),
            )
        )

    return result


def _residue_basis(
    n: np.ndarray,
    ca: np.ndarray,
    c: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Construct the Definition 3.4 orthonormal residue basis."""

    ca_to_n = n - ca
    ca_to_c = c - ca

    ca_to_n_squared = float(
        np.dot(ca_to_n, ca_to_n)
    )

    if (
        not math.isfinite(ca_to_n_squared)
        or ca_to_n_squared <= 0.0
    ):
        raise ValueError(
            "Definition 3.4 is undefined because |CA->N| is zero "
            "or non-finite"
        )

    h = (
        ca_to_c
        - (
            float(np.dot(ca_to_n, ca_to_c))
            / ca_to_n_squared
        )
        * ca_to_n
    )

    ca_to_n_norm = float(np.linalg.norm(ca_to_n))
    h_norm = float(np.linalg.norm(h))

    if (
        not math.isfinite(ca_to_n_norm)
        or ca_to_n_norm <= 0.0
    ):
        raise ValueError(
            "Definition 3.4 is undefined because |CA->N| is zero "
            "or non-finite"
        )

    if (
        not math.isfinite(h_norm)
        or h_norm <= 0.0
    ):
        raise ValueError(
            "Definition 3.4 is undefined because |h| is zero "
            "or non-finite"
        )

    u = ca_to_n / ca_to_n_norm
    v = h / h_norm
    w = np.cross(u, v)

    if not (
        np.isfinite(u).all()
        and np.isfinite(v).all()
        and np.isfinite(w).all()
    ):
        raise ValueError(
            "Definition 3.4 produced a non-finite residue basis"
        )

    return u, v, w


def _coordinates_in_basis(
    vector: np.ndarray,
    basis: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> np.ndarray:
    """Represent one vector in an orthonormal residue basis."""

    return np.asarray(
        [
            np.dot(vector, basis[0]),
            np.dot(vector, basis[1]),
            np.dot(vector, basis[2]),
        ],
        dtype=np.float64,
    )


def compute_bri(
    chain: ChainObservation,
) -> np.ndarray:
    """Compute the canonical complete BRI of one retained chain.

    The returned array has shape ``(m, 9)``.  Its columns are
    ``BRI_COLUMNS``.  Values are canonicalized with
    ``numpy.around(..., 3)`` to reproduce the pinned BRI v1.2.2
    representation.

    The input chain must already represent the exact retained
    Stage-2-eligible backbone.  No cleaning or residue removal occurs
    here.
    """

    residues = _index_backbone(chain)
    length = len(residues)

    result = np.empty(
        (length, len(BRI_COLUMNS)),
        dtype=np.float64,
    )

    _, first_n, first_ca, first_c = residues[0]

    first_basis = _residue_basis(
        first_n,
        first_ca,
        first_c,
    )

    result[0, 0:3] = _coordinates_in_basis(
        first_n - first_ca,
        first_basis,
    )
    result[0, 3:6] = _coordinates_in_basis(
        first_ca - first_ca,
        first_basis,
    )
    result[0, 6:9] = _coordinates_in_basis(
        first_c - first_ca,
        first_basis,
    )

    for index in range(1, length):
        (
            _residue_id,
            current_n,
            current_ca,
            current_c,
        ) = residues[index]

        (
            _previous_residue_id,
            previous_n,
            previous_ca,
            previous_c,
        ) = residues[index - 1]

        previous_basis = _residue_basis(
            previous_n,
            previous_ca,
            previous_c,
        )

        result[index, 0:3] = _coordinates_in_basis(
            current_n - previous_c,
            previous_basis,
        )
        result[index, 3:6] = _coordinates_in_basis(
            current_ca - current_n,
            previous_basis,
        )
        result[index, 6:9] = _coordinates_in_basis(
            current_c - current_ca,
            previous_basis,
        )

    if not np.isfinite(result).all():
        raise ValueError(
            "BRI computation produced non-finite coordinates"
        )

    return np.around(result, 3)
