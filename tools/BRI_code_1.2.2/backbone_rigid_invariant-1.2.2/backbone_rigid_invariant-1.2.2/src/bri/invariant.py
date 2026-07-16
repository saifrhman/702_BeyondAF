# -*- coding = utf-8 -*-
# @File: invariant.PY
"""Backbone Rigid Invariant (BRI) computation module.

This module provides functions for computing and analyzing backbone rigid invariants
of protein structures. It includes functionality for:

- Extending invariants with bond and torsion angles
- Generating invariant summaries
- Computing coordinate transformations and basis vectors

The module implements the mathematical foundations of the BRI method as described
in the associated research papers.
"""


from typing import Dict, List, Tuple, Union, Optional

import numpy as np
import pandas as pd

from .base.math_base import (
    get_distance,
    dot_product,
    cross_product,
    vector_norm,
    vector_round,
    get_angle,
    get_dihedral_angle,
    get_third_side_length,
)
from .base.base_util import amino_acid_short

dic: Dict[str, List[Optional[int]]] = {
    "N": [4, 7],
    "CA": [7, 10],
    "C": [10, 13],
    "O": [13, None],
}

BTP_ATTR = [
    "|C-N-A|",
    "|N-A-C|",
    "|A-C-N|",
    "TP(NA)_x",
    "TP(NA)_y",
    "TP(AC)_x",
    "TP(AC)_y",
    "TP(CN)_x",
    "TP(CN)_y",
]


def get_invariant_BTP(invariant: pd.DataFrame) -> pd.DataFrame:
    """Computes bond-torsion-point (BTP) invariant (from BRI).

    Bond-torsion-point (BTP) invariant is a complete invariant simpler than both BRI and
    AL. It consists of three 2-bond lengths (diagonal length skipping one intermediate
    atom) :math:`C_{i-1}A_i, N_iC_i, A_iN_{i+1}`, and three torsion points TP(bond)
    = (bond_length*cos(torsion_angle), bond_length*sin(torsion_angle)) for 3 bonds
    :math:`C_{i-1}N_i, N_iA_i, A_iC_i` and corresponding torsion angles at these bonds.

    This BTP invariant saves 9 mentioned values per residue as: ``|C-N-A|``, ``|N-A-C|``,
    ``|A-C-N|``, ``TP(CN)_x``, ``TP(CN)_y``, ``TP(NA)_x``, ``TP(NA)_y``, ``TP(AC)_x``,
    ``TP(AC)_y``.

    :param invariant: BRI of a chain with angles
    :return: bond-torsion-point (BTP) invariant
    """
    columns = [
        "length(N)",
        "length(A)",
        "length(C)",
        "angle(N)",
        "angle(A)",
        "angle(C)",
        "tau(NA)",
        "tau(AC)",
        "tau(CN)",
    ]
    tmp = invariant.copy()[columns].astype("float")
    property_this = ["length(N)", "angle(C)", "tau(CN)"]
    property_next = ["length(N_1)", "angle(C_1)", "tau(CN_1)"]
    tmp[property_next] = tmp[property_this].shift(-1)

    # |C-N-A|, |N-A-C|, |A-C-N|
    invariant["|C-N-A|"] = tmp.apply(
        lambda row: get_third_side_length(
            row["length(N)"], row["length(A)"], row["angle(N)"]
        ),
        axis=1,
    )
    invariant["|N-A-C|"] = tmp.apply(
        lambda row: get_third_side_length(
            row["length(A)"], row["length(C)"], row["angle(A)"]
        ),
        axis=1,
    )
    invariant["|A-C-N|"] = tmp.apply(
        lambda row: get_third_side_length(
            row["length(C)"], row["length(N_1)"], row["angle(C)"]
        ),
        axis=1,
    )

    # TP(CN)_x,TP(CN)_y, TP(NA)_x,TP(NA)_y,TP(AC)_x,TP(AC)_y
    # TP(bond) = (|bond_length|*cos(torsion_angle), |bond_length|*sin(torsion_angle))
    invariant["TP(NA)_x"] = tmp.apply(
        lambda row: row["length(A)"] * np.cos(np.deg2rad(row["tau(NA)"])), axis=1
    )
    invariant["TP(NA)_y"] = tmp.apply(
        lambda row: row["length(A)"] * np.sin(np.deg2rad(row["tau(NA)"])), axis=1
    )
    invariant["TP(AC)_x"] = tmp.apply(
        lambda row: row["length(C)"] * np.cos(np.deg2rad(row["tau(AC)"])), axis=1
    )
    invariant["TP(AC)_y"] = tmp.apply(
        lambda row: row["length(C)"] * np.sin(np.deg2rad(row["tau(AC)"])), axis=1
    )
    invariant["TP(CN)_x"] = tmp.apply(
        lambda row: row["length(N_1)"] * np.cos(np.deg2rad(row["tau(CN)"])), axis=1
    )
    invariant["TP(CN)_y"] = tmp.apply(
        lambda row: row["length(N_1)"] * np.sin(np.deg2rad(row["tau(CN)"])), axis=1
    )

    inv_cols = list(invariant.columns)
    invariant.iloc[0, inv_cols.index("TP(NA)_x")] = tmp.iloc[0]["length(A)"]
    invariant.iloc[-1, inv_cols.index("TP(AC)_x")] = tmp.iloc[-1]["length(C)"]

    invariant[BTP_ATTR] = invariant[BTP_ATTR].apply(vector_round, decimal=2)
    return invariant


def invariant_ext(invariant: pd.DataFrame) -> pd.DataFrame:
    """Extend an input invariant by computing bond and torsion angles from it.

    :param invariant: A backbone-rigid-invariant with complete 3 weaker invariants and 9 stronger invariants at least.
    :return: Extended backbone-rigid-invariant with bond and torsion angles.
    """

    columns = [
        "x(AN)",
        "x(AC)",
        "y(AC)",
        "x(N)",
        "y(N)",
        "z(N)",
        "x(A)",
        "y(A)",
        "z(A)",
        "x(C)",
        "y(C)",
        "z(C)",
    ]
    tmp = invariant[columns].astype("float")
    tmp[["y(AN)", "z(AN)", "z(AC)"]] = [0.0, 0.0, 0.0]

    AN_this = ["x(AN)", "y(AN)", "z(AN)"]
    AC_this = ["x(AC)", "y(AC)", "z(AC)"]
    C_0N = ["x(N)", "y(N)", "z(N)"]
    NA = ["x(A)", "y(A)", "z(A)"]
    AC = ["x(C)", "y(C)", "z(C)"]
    CN_2 = ["x(N_1)", "y(N_1)", "z(N_1)"]
    N_2A_2 = ["x(A_1)", "y(A_1)", "z(A_1)"]
    tmp[CN_2 + N_2A_2] = tmp[C_0N + NA].shift(-1)

    # angle(A), angle(N), angle(C)
    invariant["angle(N)"] = tmp.apply(
        lambda row: get_angle(-1 * row[C_0N], row[NA]), axis=1
    )
    invariant["angle(A)"] = tmp.apply(
        lambda row: get_angle(row[AN_this], row[AC_this]), axis=1
    )
    invariant["angle(C)"] = tmp.apply(
        lambda row: get_angle(-1 * row[AC_this], row[CN_2]), axis=1
    )
    invariant[["angle(A)", "angle(N)", "angle(C)"]] = invariant[
        ["angle(A)", "angle(N)", "angle(C)"]
    ].apply(vector_round, decimal=2)

    # tau(NA), tau(AC), tau(CN)
    invariant["tau(NA)"] = tmp.apply(
        lambda row: get_dihedral_angle(row[C_0N], row[NA], row[AC]), axis=1
    )
    invariant["tau(AC)"] = tmp.apply(
        lambda row: get_dihedral_angle(-1 * row[AN_this], row[AC_this], row[CN_2]),
        axis=1,
    )
    invariant["tau(CN)"] = tmp.apply(
        lambda row: get_dihedral_angle(row[AC_this], row[CN_2], row[N_2A_2]), axis=1
    )
    invariant[["tau(NA)", "tau(AC)", "tau(CN)"]] = invariant[
        ["tau(NA)", "tau(AC)", "tau(CN)"]
    ].apply(vector_round, decimal=2)
    return invariant


def get_invariant_summary(
    invariant: pd.DataFrame, bond_length: bool = True, min_and_max: bool = False
) -> pd.DataFrame:
    """Compute backbone-rigid-invariant summary consisting of averages and standard deviations.

    :param invariant: Backbone-Rigid-Invariant of a chain
    :param bond_length: If input BRI includes bond length. Default True.
    :param min_and_max: If to compute minimum and maximum of input BRI. Default False.
    :return: Summary of input invariant.
    """
    columns = [
        "x(AN)",
        "x(AC)",
        "y(AC)",
        "x(N)",
        "y(N)",
        "z(N)",
        "x(A)",
        "y(A)",
        "z(A)",
        "x(C)",
        "y(C)",
        "z(C)",
    ]
    if bond_length:
        columns = [
            "x(AN)",
            "x(AC)",
            "y(AC)",
            "x(N)",
            "y(N)",
            "z(N)",
            "x(A)",
            "y(A)",
            "z(A)",
            "x(C)",
            "y(C)",
            "z(C)",
            "length(N)",
            "length(A)",
            "length(C)",
        ]
    invariant_data = invariant[columns].values.astype("float64")

    mean = vector_round(np.average(invariant_data[1:, 3:12], axis=0))
    dev = vector_round(np.std(invariant_data[1:, 3:12], axis=0))
    weak_mean = vector_round(np.average(invariant_data[:, 2:5], axis=0))
    weak_dev = vector_round(np.std(invariant_data[:, 2:5], axis=0))
    res = np.hstack((weak_mean, mean, weak_dev, dev)).reshape(1, -1)

    if bond_length:
        length_mean = np.average(invariant_data[:, 13:15], axis=0)
        length_mean = vector_round(
            np.hstack((np.average(invariant_data[1:, 12], axis=0), length_mean))
        )
        length_dev = np.std(invariant_data[:, 13:15], axis=0)
        length_dev = vector_round(
            np.hstack((np.std(invariant_data[1:, 12], axis=0), length_dev))
        )
        res = np.hstack(
            (weak_mean, mean, length_mean, weak_dev, dev, length_dev)
        ).reshape(1, -1)

    res_col = [col + "_mean" for col in columns] + [col + "_dev" for col in columns]

    if min_and_max:
        invariant_min = np.min(invariant_data[1:, 3:12], axis=0)
        invariant_max = np.max(invariant_data[1:, 3:12], axis=0)
        weak_invariant_min = np.min(invariant_data[:, 2:5], axis=0)
        weak_invariant_max = np.max(invariant_data[:, 2:5], axis=0)
        res = np.hstack(
            (res, weak_invariant_min, invariant_min, weak_invariant_max, invariant_max)
        ).reshape(1, -1)
        res_col = (
            res_col
            + [col + "_min" for col in columns]
            + [col + "_max" for col in columns]
        )

    summary = pd.DataFrame(data=res, columns=res_col)
    summary["chain_length"] = len(invariant)
    return summary


def get_invariant(chain: pd.DataFrame, ext: bool = False) -> pd.DataFrame:
    """Computes backbone-rigid-invariant

    :param chain: A :class:`pandas.DataFrame` that represents a chain with atomic coordinates.
    :param ext: Include bond angles and torsion angles. Default False.
    :return: 3-decimal alpha-invariant of a chain
    """

    # reconstruct DataFrame format
    chain["residue_label"] = chain["residue_label"].map(amino_acid_short)
    chain["atom"] = chain["atom"].map({"CA": "A", "N": "N", "C": "C"})
    residue_table = chain.pivot(
        index=["model_id", "chain_id", "residue_id", "residue_label"],
        columns=["atom"],
        values=["x", "y", "z"],
    ).reset_index()
    residue_table.columns = [
        f"{col[0]}({col[1]})" if col[1] else col[0] for col in residue_table.columns
    ]
    residue_table = (
        residue_table.sort_values("residue_id")
        .iloc[:, [0, 1, 2, 3, 6, 9, 12, 4, 7, 10, 5, 8, 11]]
        .reset_index(drop=True)
    )
    columns = list(residue_table.columns)
    columns = columns[:4] + ["x(AN)", "x(AC)", "y(AC)"] + columns[4:]
    # data value extraction for computation
    length = len(residue_table)
    residue_table = residue_table.to_numpy()
    invariant = np.copy(residue_table)
    extra_info = np.copy(residue_table)

    # first row
    first = residue_table[0]
    invariant[0] = row_trans(first, first, orthogonal=True)
    extra_info[0] = invariant[0]
    # rest row
    if length > 1:
        invariant[1:] = [
            row_trans(r, lr, orthogonal=True)
            for r, lr in zip(residue_table[1:, :], residue_table[0:-1, :])
        ]
        extra_info[1:] = [
            row_trans(r, r, orthogonal=True) for r in residue_table[1:, :]
        ]

    # construct DataFrame
    invariant = np.hstack(
        (invariant[:, :4], extra_info[:, [4, 10, 11]], invariant[:, 4:])
    )
    invariant[:, 4:] = vector_round(invariant[:, 4:].astype(float), 3)
    invariant = pd.DataFrame(invariant, columns=columns)

    # bond length
    for a in ("N", "A", "C"):
        invariant[f"length({a})"] = invariant[[f"x({a})", f"y({a})", f"z({a})"]].apply(
            get_distance, axis=1
        )
        invariant[f"length({a})"] = vector_round(invariant[f"length({a})"])

    invariant.iloc[0, invariant.columns.get_loc("length(A)")] = invariant.iloc[
        0, invariant.columns.get_loc("length(N)")
    ]
    invariant.iloc[0, invariant.columns.get_loc("length(N)")] = np.NAN
    invariant["chain_length"] = length

    if ext:
        return invariant_ext(invariant)
    return invariant


def row_trans(
    row: pd.Series,
    last_row: Optional[pd.Series] = None,
    update_atoms: List[str] = ["N", "CA", "C"],
    o: str = "CA",
    orthogonal: bool = True,
) -> pd.Series:
    """Transform coordinates to a new basis.

    :param row: Current row of coordinates
    :param last_row: Previous row of coordinates (optional)
    :param update_atoms: List of atoms to update, defaults to ["N", "CA", "C"]
    :param o: Origin atom, defaults to "CA"
    :param orthogonal: Whether to use orthogonal basis, defaults to True
    :return: Transformed coordinates
    """

    backbone_order = ("N", "CA", "C")
    # get the coordinate basis from the previous residue
    orthonormal_basis = get_basis(last_row, orthogonal=orthogonal)
    # atoms to represent
    res_row = np.copy(row)
    if np.array_equal(row, last_row):
        # fix origin as CA
        origin = last_row[dic[o][0] : dic[o][1]]
        # iterate atoms N, CA, C
        for a in update_atoms:
            # compute vectors originate from CA, i.e., CA-N, CA-CA, CA-N
            p = (row[dic[a][0] : dic[a][1]] - origin).astype("float64")
            # compute new coordinates and construct results
            c1, c2, c3 = [dot_product(p, v) for v in orthonormal_basis]
            res_row[dic[a][0] : dic[a][1]] = [c1, c2, c3]
        return res_row

    for a in update_atoms:
        # select vector origin as the atom before the current one
        origin_atom = backbone_order[backbone_order.index(a) - 1]
        if a == "N":
            # compute vector origin for C-N+1
            origin = last_row[dic[origin_atom][0] : dic[origin_atom][1]]
        else:
            # select origins for vectors N+1-CA+1, CA+1-C+1
            origin = row[dic[origin_atom][0] : dic[origin_atom][1]]
        # compute vector
        p = (row[dic[a][0] : dic[a][1]] - origin).astype("float64")
        # compute new coordinates and construct results
        c1, c2, c3 = [dot_product(p, v) for v in orthonormal_basis]
        res_row[dic[a][0] : dic[a][1]] = [c1, c2, c3]

    return res_row


def get_basis(
    row: pd.Series,
    o: str = "CA",
    v1_p: str = "N",
    v2_p: str = "C",
    orthogonal: bool = False,
) -> List[np.ndarray]:
    """Get basis vectors for coordinate transformation.

    :param row: Row of coordinates
    :param o: Origin atom, defaults to "CA"
    :param v1_p: First vector point, defaults to "N"
    :param v2_p: Second vector point, defaults to "C"
    :param orthogonal: Whether to use orthogonal basis, defaults to True
    :return: List of 3 basis vectors
    """

    # coordinates of CA
    origin = row[dic[o][0] : dic[o][1]]
    # vector CA-N and CA-C
    v1 = (row[dic[v1_p][0] : dic[v1_p][1]] - origin).astype("float64")
    v2 = (row[dic[v2_p][0] : dic[v2_p][1]] - origin).astype("float64")
    if orthogonal:
        v2 = v2 - (dot_product(v1, v2) / dot_product(v1, v1)) * v1
    # normalise v1, v2
    v1 = v1 / vector_norm(v1)
    v2 = v2 / vector_norm(v2)
    v3 = cross_product(v1, v2)
    return v1, v2, v3
