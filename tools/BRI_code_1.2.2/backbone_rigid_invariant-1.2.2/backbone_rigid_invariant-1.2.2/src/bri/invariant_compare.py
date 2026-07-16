# -*- coding = utf-8 -*-
# @File: pdb_util.PY
"""Backbone Rigid Invariants comparison module.

This module provides functionality for comparing protein structures using their
backbone rigid invariants (BRI). It implements efficient nearest neighbor search
algorithms to:

- Compute pairwise distances between BRI descriptors
- Compare protein sequences when requested
- Find structurally similar chains based on BRI similarity

The module enables the L-infinity (Chebyshev) and RMSD distance for BRI comparison and
Hamming distance for sequence comparison.
"""
from typing import Optional

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.neighbors import BallTree

from .pdbx2df import Chain
from .invariant import invariant_ext
from .base.math_base import get_RMSD


INVARIANT_VALUE_COLS = [
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


def extract_chain_info(df: pd.DataFrame) -> pd.DataFrame:
    """Extract chain identification and residue sequence

    :param df: Set of chains
    :return: Unique chain identificaiton with residue sequence
    """

    df["residue_label_ascii"] = df["residue_label"].apply(ord)
    chains = (
        df.groupby(["pdb_id", "model_id", "chain_id", "chain_length"], sort=False)
        .agg({"residue_label": lambda x: "".join(x)})
        .reset_index()
    )
    chains = chains.rename(columns={"residue_label": "seq"})
    return chains


def organise_dup_output(neighbours: tuple, chains: pd.DataFrame) -> pd.DataFrame:
    """Formulate duplicate search outcomes.

    :param neighbours: Tuple of object index and distance
    :param chains: Table of chain identifications
    :return: Organized distance outcomes
    """
    res = []
    for distance, idx in zip(neighbours[0][0], neighbours[1][0]):
        row = {}
        chain = chains.iloc[idx, :]

        row["distance"] = distance
        row["pdb_id"] = chain["pdb_id"]
        row["model_id"] = chain["model_id"]
        row["chain_id"] = chain["chain_id"]

        res.append(row)

    df_output = pd.DataFrame(res)
    return df_output


def generate_index_table(neighbours: tuple, reduce: bool = True) -> pd.DataFrame:
    """Generate index table from nearest neighbors query results.

    Converts the (distances, indices) tuple from kd-tree/ball-tree query into
    a DataFrame with columns: distance, idx_1, idx_2.

    :param neighbours: Tuple of (distances, indices) arrays from tree.query()
                      distances[i] contains distances for query point i
                      indices[i] contains neighbor indices for query point i
    :param reduce: If to reduce the distance matrix when generating the index table
    :return: DataFrame with columns: distance, idx_1, idx_2, filtered to idx_1 < idx_2
    """
    distances, indices = neighbours
    n_chains = len(distances)

    # Build list of rows, filtering to idx_1 < idx_2 during construction for efficiency
    rows = []
    for i in range(n_chains):
        for dist, idx in zip(distances[i], indices[i]):
            if dist == np.inf:
                break
            if reduce:  # Only include pairs where idx_1 < idx_2
                if i < idx:
                    rows.append({"distance": dist, "idx_1": i, "idx_2": idx})
            else:
                rows.append({"distance": dist, "idx_1": i, "idx_2": idx})

    return pd.DataFrame(rows)


def convert_index_table(
    idx_table: pd.DataFrame, chains: pd.DataFrame, other: pd.DataFrame, keys: list
) -> pd.DataFrame:
    """Convert index table to information table with chain identifiers.

    Converts index-based table to a table with actual chain identifiers
    using vectorized pandas operations.

    :param idx_table: DataFrame with columns: distance, idx_1, idx_2
    :param chains: DataFrame with chain information (rows indexed by idx_1)
    :param other: DataFrame with chain information (rows indexed by idx_2)
    :param keys: List of column names to extract from chains and other
    :return: DataFrame with distance and chain information from both tables.
             Columns: distance, {key}1 for each key in keys, {key}2 for each key in keys
    """
    # Use merge operations for efficiency instead of iterrows()
    result = idx_table[["distance", "idx_1", "idx_2"]].copy()

    # Merge with chains to get information for idx_1
    chains_subset = chains[keys].reset_index(drop=True)
    chains_subset.columns = [f"{key}1" for key in keys]
    chains_subset["idx_1"] = chains_subset.index
    result = result.merge(chains_subset, on="idx_1", how="left")

    # Merge with other to get information for idx_2
    other_subset = other[keys].reset_index(drop=True)
    other_subset.columns = [f"{key}2" for key in keys]
    other_subset["idx_2"] = other_subset.index
    result = result.merge(other_subset, on="idx_2", how="left")

    # Reorder columns: distance first, then keys with suffix 1, then keys with suffix 2
    col_order = ["distance"] + [f"{key}1" for key in keys] + [f"{key}2" for key in keys]
    return result[col_order]


def coordinate_value_reshape(
    df: pd.DataFrame, chain_len: int, cols: list = INVARIANT_VALUE_COLS
) -> np.ndarray:
    """Reshape atom coordinates of chains into long vectors

    :param df: Set of chain coordinates
    :param chain_len: Chian length
    :return: Long vectors of atom coordinates
    """
    coordinate_value = df.loc[:, cols].values
    n_chains = len(coordinate_value) // chain_len
    vector_len = len(cols) * chain_len
    shaped_value = np.reshape(coordinate_value, (n_chains, vector_len))
    return shaped_value


def group_invariant_compare(
    x: pd.DataFrame, y: Optional[pd.DataFrame] = None, seq_compare=False
):
    """Compare backbone rigid invariants.

    This function performs pairwise comparison of BRI descriptors. It uses k-nearest neighbors with the Chebyshev distance metric
    for structural comparison, and optionally the Hamming distance for sequence
    comparison.

    :param x: DataFrame of BRIs as base with columns:
              - pdb_id: PDB identifier
              - model_id: Model number
              - chain_id: Chain identifier
              - residue_id: Residue number
              - residue_label: One-letter amino acid code
              - x_N to z_C: 9 coordinate columns
              - chain_length: Number of residues
    :param y: DataFrame of BRIs as comparison target. If `None`, pairwise comparison will be conducted inside `x`.
    :param seq_compare: Whether to also compare sequences (default: False)
    :return: DataFrame containing pairwise comparison results with columns:
             - pdb_id1, model_id1, chain_id1: First chain identifiers
             - pdb_id2, model_id2, chain_id2: Second chain identifiers
             - L_inf_invariant: BRI distance
             If seq_compare is True, additional columns:
             - seq_diff: Number of different residues
             - seq1, seq2: Sequences being compared
             Returns None if no comparisons can be made
    :raises: None, returns None for invalid inputs
    """
    same_col = ["pdb_id", "model_id", "chain_id"]
    seq_compare_col = ["pdb_id", "model_id", "chain_id", "seq"]
    reduce_matrix = False
    if y is None:
        y = x
        reduce_matrix = True

    chains_x = extract_chain_info(x)
    chains_y = extract_chain_info(y)
    n_base = len(chains_x)
    n_target = len(chains_y)
    if (n_base < 1) or (n_target < 1):  # no enough chains
        return None
    chain_lengths = list(chains_x["chain_length"].unique())
    if len(chain_lengths) > 1:  # chain length not equal
        return None
    chain_length = chain_lengths[0]

    # reshape atom coordinates in a chain into long vectors
    base_value = coordinate_value_reshape(x, chain_length)
    search_value = coordinate_value_reshape(y, chain_length)
    # get nn with kd_tree using Chebyshev distance
    tree = cKDTree(base_value)
    nn_res = tree.query(search_value, k=len(base_value), p=np.inf, workers=-1)

    # Compute structural distance (BRI) comparisons
    idx_table = generate_index_table(nn_res, reduce_matrix)
    distance_res = convert_index_table(idx_table, chains_x, chains_x, same_col)

    if not seq_compare:
        return distance_res

    # Also compute sequence distance comparisons
    base_seq_value = coordinate_value_reshape(x, chain_length, ["residue_label_ascii"])
    targ_seq_value = coordinate_value_reshape(y, chain_length, ["residue_label_ascii"])
    seq_tree = BallTree(base_seq_value, metric="hamming")
    nn_res_seq = seq_tree.query(targ_seq_value, k=n_base, return_distance=True)

    idx_table_seq = generate_index_table(nn_res_seq, reduce_matrix)
    result_seq = convert_index_table(idx_table_seq, chains_x, chains_x, seq_compare_col)

    # Convert normalized Hamming distance to absolute sequence difference
    result_seq = result_seq.rename(columns={"distance": "seq_diff"})
    result_seq["seq_diff"] = (
        (result_seq["seq_diff"] * chain_length).round().astype("int")
    )

    # Merge structural and sequence distance results
    merge_keys = [f"{col}1" for col in same_col] + [f"{col}2" for col in same_col]
    output_full = distance_res.merge(result_seq, on=merge_keys, how="inner")
    return output_full


def neighbour_distance_search_ckdtree(
    target: pd.DataFrame, df: pd.DataFrame, seq_compare: bool = False
) -> pd.DataFrame:
    """Query Chebyshev(L_infinite) distance of one BRI against a set of BRIs

    :param target: BRI of a chain
    :param df: A set of BRIs
    :param seq_compare: Whether to also compare sequences, defaults to False
    :return: L_inf distance of the target to chains in the set
    """
    chains = extract_chain_info(df)
    n_neighbor = len(chains)
    if n_neighbor < 1:  # no neighbors
        return None
    chain_lengths = list(chains["chain_length"].unique())
    if len(chain_lengths) > 1:  # chain length not equal
        return None
    chain_length = chain_lengths[0]

    # reshape atom coordinates in a chain into long vectors
    coordinate_value = coordinate_value_reshape(df, chain_length)
    target_value = coordinate_value_reshape(target, chain_length)

    # Build kd_tree using L_infinite distance
    kd_tree = cKDTree(coordinate_value)
    neighbours = kd_tree.query(target_value, k=n_neighbor, p=np.inf)

    # organise output
    df_output = organise_dup_output(neighbours, chains)

    if not seq_compare:
        return df_output

    else:  # KNN with kd_tree and hamming distance
        df_seq_value = coordinate_value_reshape(
            df, chain_length, ["residue_label_ascii"]
        )
        target["residue_label_ascii"] = target["residue_label"].apply(ord)
        target_seq_value = coordinate_value_reshape(
            target, chain_length, ["residue_label_ascii"]
        )
        tree = BallTree(df_seq_value, metric="hamming")
        # make distance matrix with NN
        neighbours = tree.query(target_seq_value, k=n_neighbor, return_distance=True)

        df_seq_diff = organise_dup_output(neighbours, chains)
        df_seq_diff.columns = ["seq_diff", "pdb_id", "model_id", "chain_id"]
        df_seq_diff["seq_diff"] = np.around(df_seq_diff["seq_diff"] * chain_length)
        df_seq_diff["seq_diff"] = df_seq_diff["seq_diff"].astype("int")
        output = df_output.merge(df_seq_diff, on=["pdb_id", "model_id", "chain_id"])
        return output


def neighbour_distance_search_RMSD(
    target: pd.DataFrame, df: pd.DataFrame, seq_compare: bool = False
) -> pd.DataFrame:
    """Query RMSD distance of one BRI against a set of BRIs

    :param target: BRI of a chain
    :param df: A set of BRIs
    :param seq_compare: Whether to also compare sequences, defaults to False
    :return: RMSD distance of the target to chains in the set
    """
    chains = extract_chain_info(df)
    n_neighbor = len(chains)
    if n_neighbor < 1:  # no neighbors
        return None
    chain_lengths = list(chains["chain_length"].unique())
    if len(chain_lengths) > 1:  # chain length not equal
        return None
    chain_length = chain_lengths[0]

    # reshape atom coordinates in a chain into long vectors
    coordinate_value = df[INVARIANT_VALUE_COLS].values
    coordinate_value = coordinate_value.reshape((n_neighbor, 9 * chain_length))
    target_value = target[INVARIANT_VALUE_COLS].values.reshape((-1, 9 * chain_length))

    # make distance matrix with BallTree using RMSD
    tree = BallTree(coordinate_value, metric=get_RMSD)
    neighbours = tree.query(target_value, k=n_neighbor, return_distance=True)

    # organise output
    df_output = organise_dup_output(neighbours, chains)
    if not seq_compare:
        return df_output

    else:  # KNN with kd_tree and hamming distance
        df_seq_value = coordinate_value_reshape(
            df, chain_length, ["residue_label_ascii"]
        )
        target["residue_label_ascii"] = target["residue_label"].apply(ord)
        target_seq_value = coordinate_value_reshape(
            target, chain_length, ["residue_label_ascii"]
        )
        tree = BallTree(df_seq_value, metric="hamming")
        # make distance matrix with NN
        neighbours = tree.query(target_seq_value, k=n_neighbor, return_distance=True)

        df_seq_diff = organise_dup_output(neighbours, chains)
        df_seq_diff.columns = ["seq_diff", "pdb_id", "model_id", "chain_id"]
        df_seq_diff["seq_diff"] = np.around(df_seq_diff["seq_diff"] * chain_length)
        df_seq_diff["seq_diff"] = df_seq_diff["seq_diff"].astype("int")
        output = df_output.merge(df_seq_diff, on=["pdb_id", "model_id", "chain_id"])
        return output


def get_theo_lipschitz_constant(chain: Chain, chain_perturbed: Chain) -> float:
    """Compute theoretical Lipschitz constant.

    :param chain: A :class:`bri.Chain <.pdbx2df.Chain>` object as a reference.
    :param chain_perturbed: A (perturbed) :class:`bri.Chain <.pdbx2df.Chain>` object.
    :return: The theoretical Lipschitz constant.
    """

    bri_ext = invariant_ext(chain.invariant)
    bri_ext_perturb = invariant_ext(chain_perturbed.invariant)
    bri_ext["h_C"] = bri_ext["length(C)"] * np.sin(np.deg2rad(bri_ext["angle(A)"]))
    bri_ext_perturb["h_C"] = bri_ext_perturb["length(C)"] * np.sin(
        np.deg2rad(bri_ext_perturb["angle(A)"])
    )

    length_cols = ["length(N)", "length(A)", "length(C)"]
    L = max(bri_ext[length_cols].max().max(), bri_ext_perturb[length_cols].max().max())

    L_AC = max(bri_ext["length(C)"].max(), bri_ext_perturb["length(C)"].max())
    l_NA = min(bri_ext["length(A)"].min(), bri_ext_perturb["length(A)"].min())
    l_h_C = min(bri_ext["h_C"].min(), bri_ext_perturb["h_C"].min())
    K = 1 / l_NA + 2 / l_h_C * (1 + 2 * L_AC / l_NA)

    return 2 * (1 + 2 * L * K)


def get_e_lipschitz_constant(chain: Chain, chain_perturbed: Chain) -> float:
    """Compute experimental Lipschitz constant.

    :param chain: A :class:`bri.Chain <.pdbx2df.Chain>` object as a reference.
    :param chain_perturbed: A (perturbed) :class:`bri.Chain <.pdbx2df.Chain>` object.
    :return: The experimental Lipschitz constant.
    """

    inv_cols = ["x(N)", "y(N)", "z(N)", "x(A)", "y(A)", "z(A)", "x(C)", "y(C)", "z(C)"]
    coordinate_cols = ["x", "y", "z"]

    linf_bri = (
        (chain.invariant[inv_cols] - chain_perturbed.invariant[inv_cols])
        .abs()
        .max()
        .max()
    )
    eps_max = (
        (
            chain.coordinates[coordinate_cols]
            - chain_perturbed.coordinates[coordinate_cols]
        )
        .abs()
        .max()
        .max()
    )
    return linf_bri / eps_max
