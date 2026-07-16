"""Protein structure filtering and validation module.

This module provides comprehensive validation and filtering functionality for protein structures.
It performs the following checks:

1. Eligibility check - verifies if all important information is included
2. Disorder detection - checks for disordered atoms
3. Clash detection - identifies atomic clashes
4. Gap detection - identifies atom gaps
5. Residue continuity - checks for missing or non-contiguous residues
6. Atom completeness - verifies if all important atoms are present in residues
7. Standard residue validation - checks for non-standard amino acids

The module supports both entry-level and chain-level operations through decorators
and provides integrated cleaning pipelines for both mini entries and full entries.
"""

import itertools
from typing import List, Tuple, Dict, Optional, Union, Any, Set

import pandas as pd
import numpy as np

from .base.base_util import (
    basic_amino_acid_20,
    basic_amino_acid_20_s,
    amino_acid_short,
)
from .base.math_base import get_distance, get_angle
from .pdbx2df import MiniChain, MiniEntry, Entry, on_entry


atom: List[str] = ["N", "CA", "C", "N+1"]
min_test_set: List[Tuple[str, str]] = list(itertools.combinations(atom, 2))
atom_position: Dict[str, int] = {"N": 0, "CA": 1, "C": 2, "N+1": 3, "CA+1": 4, "C+1": 5}
atom_combination = itertools.combinations(atom_position.keys(), 2)
max_test_set: List[Tuple[str, str]] = [
    (a1, a2)
    for a1, a2 in atom_combination
    if (abs(atom_position[a2] - atom_position[a1]) < 4)
]
max_test_set = [
    (a1, a2) for a1, a2 in max_test_set if not (a1.endswith("+1") and a2.endswith("+1"))
]
max_test_dict: Dict[Tuple[str, str], int] = {
    (a1, a2): abs(atom_position[a2] - atom_position[a1]) for a1, a2 in max_test_set
}
ENTRY_CLEAN_COL = [
    "pdb_id",
    "entity_id",
    "model_id",
    "chain_id",
    "start_residue",
    "chain_length",
    "auth_chain_id",
    "auth_seq_id_start",
    "auth_seq_id_end",
    "seq",
]
MINI_ENTRY_CLEAN_COL = [
    "pdb_id",
    "model_id",
    "chain_id",
    "start_residue",
    "chain_length",
]


@on_entry()
def entry_disorder_check(pdb: MiniEntry) -> Optional[pd.DataFrame]:
    """Check for disordered atoms in the entry.

    :param pdb: The entry to check
    :return: DataFrame of disordered atoms if found, None otherwise
    """
    return disorder_check(pdb)


@on_entry()
def entry_clash_check(pdb: MiniEntry, lower_bound: float = 1) -> Optional[pd.DataFrame]:
    """Check for atomic clashes in the entry.

    :param pdb: The entry to check
    :return: DataFrame of clashing atoms if found, None otherwise
    """
    return clash_check(pdb, lower_bound)


@on_entry()
def entry_gap_check(
    pdb: MiniEntry, upper_bound_coefficient: float = 2
) -> Optional[pd.DataFrame]:
    """Check for gaps in the entry.

    :param pdb: The entry to check
    :param upper_bound_coefficient: Coefficient for gap detection
    :return: DataFrame of gaps if found, None otherwise
    """
    return gap_check(pdb, upper_bound_coefficient)


@on_entry()
def entry_residue_continuity_check(pdb: MiniEntry) -> Optional[pd.DataFrame]:
    """Check for residue continuity in the entry.

    :param pdb: The entry to check
    :return: DataFrame of discontinuous residues if found, None otherwise
    """
    return residue_continuity_check(pdb)


@on_entry()
def entry_residue_completeness_check(pdb: MiniEntry) -> Optional[pd.DataFrame]:
    """Check for residue completeness in the entry.

    :param pdb: The entry to check
    :return: DataFrame of incomplete residues if found, None otherwise
    """
    return residue_completeness_check(pdb)


@on_entry()
def entry_standard_residue_check(
    pdb: MiniEntry, label_len: int = 3
) -> Optional[pd.DataFrame]:
    """Check for standard residues in the entry.

    :param pdb: The entry to check
    :param label_len: Length of residue labels to check
    :return: DataFrame of non-standard residues if found, None otherwise
    """
    return standard_residue_check(pdb, label_len)


@on_entry()
def entry_residue_angel_check(
    pdb: MiniEntry, origin: str = "CA", points: Tuple[str, str] = ("C", "N")
) -> Optional[pd.DataFrame]:
    """Check residue angles in the entry.

    :param pdb: The entry to check
    :param origin: Origin atom for angle calculation
    :param points: Points to calculate angle between
    :return: DataFrame of residues with invalid angles if found, None otherwise
    """
    return residue_angel_check(pdb, origin, points)


def minientry_integrated_cleaning(
    pdb_id: str,
) -> Tuple[Optional[pd.DataFrame], pd.DataFrame]:
    """Perform integrated cleaning on a mini entry.

    This function performs a complete cleaning pipeline on a mini entry, including:

    - Disorder detection
    - Clash detection
    - Gap detection
    - Residue continuity check
    - Atom completeness check
    - Standard residue validation

    :param pdb_id: PDB ID of the structure to clean
    :return: Tuple of (clean_set, dirty_set) DataFrames, where clean_set contains
            the cleaned chains and dirty_set contains the filtered residues
    """
    entry = MiniEntry(pdb_id)
    if not isinstance(entry, MiniEntry):
        return None, pd.DataFrame(entry, index=[0])

    entry_integrated_filter = on_entry()(integrated_chainwise_filter)
    cleaning_res = entry_integrated_filter(entry)
    if cleaning_res is None:
        return None

    dirty_set = cleaning_res[cleaning_res["type"] != "clean"]
    if "chain_length" in dirty_set.columns:
        dirty_set = dirty_set[
            [
                "pdb_id",
                "model_id",
                "chain_id",
                "residue_id",
                "type",
                "residue_label",
                "chain_length",
            ]
        ]
    else:
        dirty_set = dirty_set[
            ["pdb_id", "model_id", "chain_id", "residue_id", "residue_label", "type"]
        ]
        dirty_set["chain_length"] = np.nan

    clean_set = cleaning_res[cleaning_res["type"] == "clean"]
    clean_set = (
        clean_set.groupby(["pdb_id", "model_id", "chain_id"])
        .agg(
            start_residue=("residue_id", "min"),
            seq_max=("residue_id", "max"),
            residue_label=("residue_label", lambda x: "".join(x[::3])),
        )
        .reset_index()
    )
    clean_set["chain_length"] = clean_set.apply(
        lambda row: row.seq_max - row.start_residue + 1, axis=1
    )
    clean_set.start_residue = clean_set.start_residue.astype("int")
    clean_set.chain_length = clean_set.chain_length.astype("int")

    clean_set = clean_set.drop(columns=["seq_max"])
    clean_set = clean_set.rename(columns={"residue_label": "seq"})
    return clean_set, dirty_set


def entry_integrated_cleaning(
    pdb_id: str,
) -> Tuple[Optional[pd.DataFrame], pd.DataFrame]:
    """Perform integrated cleaning on a full entry.

    Similar to minientry_integrated_cleaning but processes a full entry with
    additional entity information.

    :param pdb_id: PDB ID of the structure to clean
    :return: Tuple of (clean_set, dirty_set) DataFrames, where clean_set contains
            the cleaned chains and dirty_set contains the filtered residues
    """
    entry = field_check(pdb_id)
    if not isinstance(entry, Entry):
        return None, pd.DataFrame(entry, index=[0])

    entry_integrated_filter = on_entry()(integrated_chainwise_filter)
    cleaning_res = entry_integrated_filter(entry)
    if cleaning_res is None:
        return None

    dirty_set = cleaning_res[cleaning_res["type"] != "clean"]
    if "chain_length" in dirty_set.columns:
        dirty_set = dirty_set[
            [
                "pdb_id",
                "model_id",
                "chain_id",
                "residue_id",
                "type",
                "residue_label",
                "chain_length",
            ]
        ]
    else:
        dirty_set = dirty_set[
            ["pdb_id", "model_id", "chain_id", "residue_id", "residue_label", "type"]
        ]
        dirty_set["chain_length"] = np.nan

    clean_set = cleaning_res[cleaning_res["type"] == "clean"]
    clean_set = (
        clean_set.groupby(["pdb_id", "model_id", "chain_id", "auth_chain_id"])
        .agg(
            start_residue=("residue_id", "min"),
            seq_max=("residue_id", "max"),
            auth_seq_id_start=("auth_residue_id", "min"),
            auth_seq_id_end=("auth_residue_id", "max"),
            residue_label=("residue_label", lambda x: "".join(x[::3])),
        )
        .reset_index()
    )
    clean_set["chain_length"] = clean_set.apply(
        lambda row: row.seq_max - row.start_residue + 1, axis=1
    )
    clean_set.start_residue = clean_set.start_residue.astype("int")
    clean_set.chain_length = clean_set.chain_length.astype("int")
    clean_set.auth_seq_id_start = clean_set.auth_seq_id_start.astype("int")
    clean_set.auth_seq_id_end = clean_set.auth_seq_id_end.astype("int")

    chain_entity_dict = {i["id"]: int(i["entity_id"]) for i in entry.label_entity_info}
    clean_set["entity_id"] = clean_set["chain_id"].map(chain_entity_dict)
    clean_set = clean_set.drop(columns=["seq_max"])
    clean_set = clean_set.rename(columns={"residue_label": "seq"})
    clean_set = clean_set[ENTRY_CLEAN_COL]
    return clean_set, dirty_set


def integrated_chainwise_filter(chain: MiniChain) -> Optional[pd.DataFrame]:
    """Integrated chain cleaning

    1. Drop disordered atoms - labeled 'disordered'
    2. Drop chains with chain breaks - labeled 'chain-break'
    3. Drop non-typical residues - labeled 'non-standard'
    4. Drop residues with missing alpha group atoms (check chain breaks by the end) - labeled 'incomplete'
    5. Drop residues with atom clashes - labeled 'clash'
    6. Drop residues with atom gaps - labeled 'gap'

    Finalized by check chain breaks - labeled 'chain-break'

    :param chain: Receives a DataFrame representing atoms from one single chain
    :return: DataFrame including atom geometric data that passes all checks, otherwise None
    """

    dirty_result = None
    chain["residue_label"] = chain["residue_label"].map(
        amino_acid_short, na_action="ignore"
    )
    full_seq = ",".join(chain.drop_duplicates("residue_id")["residue_label"])

    disordered_atoms = disorder_check(chain)
    if disordered_atoms is not None:
        # remove defect residues
        disordered_atoms["type"] = "disordered"
        dirty_result = sort_dirty_residues([dirty_result, disordered_atoms])
        chain = chain[~chain["residue_id"].isin(dirty_result["residue_id"])]
        if chain.empty:
            return dirty_result

    # check chain break
    chain_break = residue_continuity_check(chain)
    if chain_break is not None:
        chain_break["type"] = "chain-break"
        chain_break["residue_label"] = full_seq
        return sort_dirty_residues([dirty_result, chain_break])

    non_typical_atoms = standard_residue_check(chain, 1)
    defect_residue_atoms = residue_completeness_check(chain)
    if non_typical_atoms is not None:
        non_typical_atoms["type"] = "non-standard"
        dirty_result = sort_dirty_residues([dirty_result, non_typical_atoms])
    if defect_residue_atoms is not None:
        defect_residue_atoms["type"] = "incomplete"
        dirty_result = sort_dirty_residues([dirty_result, defect_residue_atoms])
    # remove defect residues
    if non_typical_atoms is not None or defect_residue_atoms is not None:
        chain = chain[~chain["residue_id"].isin(dirty_result["residue_id"])]
        # check chain break
        if chain.empty:
            return dirty_result
        chain_break = residue_continuity_check(chain)
        if chain_break is not None:
            chain_break["type"] = "chain-break"
            chain_break["residue_label"] = full_seq
            return sort_dirty_residues([dirty_result, chain_break])

    clashed_atoms = clash_check(chain, 0.01)
    # gap_atoms = gap_check(chain)
    # remove defect residues
    if clashed_atoms is not None:
        clashed_atoms.drop(columns=["x", "y", "z", "occupancy"], inplace=True)
        clashed_atoms["type"] = "clash"
        dirty_result = sort_dirty_residues([dirty_result, clashed_atoms])
    # if gap_atoms is not None:
    #     gap_atoms.drop(columns=['x', 'y', 'z', 'occupancy'], inplace=True)
    #     gap_atoms['type'] = 'gap'
    #     dirty_result = sort_dirty_residues([dirty_result, gap_atoms])
    if clashed_atoms is not None:
        # remove defect residues
        chain = chain[~chain["residue_id"].isin(dirty_result["residue_id"])]
        if chain.empty:
            return dirty_result

    # final chain break check
    chain_break = residue_continuity_check(chain)
    if chain_break is not None:
        chain_break["type"] = "chain-break"
        chain_break["residue_label"] = full_seq
        return sort_dirty_residues([dirty_result, chain_break])

    chain["type"] = "clean"
    return pd.concat([chain, dirty_result], ignore_index=True)


def disorder_check(chain: MiniChain) -> Optional[pd.DataFrame]:
    """Check for disordered atoms in a chain.

    :param chain: Chain to check for disorders
    :return: DataFrame of disordered atoms if found, None otherwise
    """
    dup = chain.loc[
        chain.duplicated(["residue_id", "chain_id", "model_id", "atom"], keep=False), :
    ]
    partial = chain[chain["occupancy"] != 1.0]
    if partial.empty and dup.empty:
        return None

    res = pd.concat([dup, partial], ignore_index=True).drop_duplicates()
    return res


def clash_check(chain: MiniChain, lower_bound: float = 1) -> Optional[pd.DataFrame]:
    """Check for atomic clashes in a chain.

    :param chain: Chain to check for clashes
    :param lower_bound: Minimum distance between atoms, defaults to 1
    :return: DataFrame of clashing atoms if found, None otherwise
    """
    output = []
    for symbols in min_test_set:
        symbol_next = []
        symbol_names = []
        for s in symbols:
            if s.endswith("+1"):
                s = s[:-2]
                symbol_next.append(s)
            symbol_names.append(s)
        # select chain with target atoms
        selected_atoms = chain[chain["atom"].isin(symbol_names)]
        if len(symbol_next) and symbol_names[0] != symbol_names[1]:
            selected_atoms.loc[
                (selected_atoms["atom"].isin(symbol_next)), "residue_id"
            ] -= 1
        atom1, atom2 = symbols
        # get length
        chain_len = len(selected_atoms["residue_id"].unique())
        if chain_len == 0:
            continue

        # distance
        if symbol_names[0] != symbol_names[1]:
            chain_group = selected_atoms.groupby("residue_id")[["x", "y", "z"]]
            selected_atoms.loc[:, "distance"] = chain_group.diff(-1).apply(
                get_distance, axis=1
            )
            del chain_group
        else:
            selected_atoms.loc[:, "distance"] = (
                selected_atoms[["x", "y", "z"]].diff(-1).apply(get_distance, axis=1)
            )
        selected_atoms.loc[:, "atom1"] = atom1
        selected_atoms.loc[:, "atom2"] = atom2

        selected_atoms = selected_atoms.dropna(how="any")
        selected_atoms = selected_atoms[selected_atoms["distance"] < lower_bound]

        selected_atoms = selected_atoms.drop_duplicates(keep="first")
        output.append(selected_atoms)
    output = pd.concat(output, ignore_index=True)
    if output.empty:
        return None
    return output


def gap_check(
    chain: MiniChain, upper_bound_coefficient: float = 2
) -> Optional[pd.DataFrame]:
    """Check for gaps between atoms in a chain.

    :param chain: Chain to check for gaps
    :param upper_bound_coefficient: Coefficient for gap detection, defaults to 2
    :return: DataFrame of gaps if found, None otherwise
    """
    output = []
    for symbols, distance in max_test_dict.items():
        symbol_next = []
        symbol_names = []
        for s in symbols:
            if s.endswith("+1"):
                s = s[:-2]
                symbol_next.append(s)
            symbol_names.append(s)
        # select chain with target atoms
        selected_atoms = chain[chain["atom"].isin(symbol_names)]
        if len(symbol_next) and symbol_names[0] != symbol_names[1]:
            selected_atoms.loc[
                (selected_atoms["atom"].isin(symbol_next)), "residue_id"
            ] -= 1
        atom1, atom2 = symbols
        # get length
        chain_len = len(selected_atoms["residue_id"].unique())
        if chain_len == 0:
            continue

        # distance
        if symbol_names[0] != symbol_names[1]:
            chain_group = selected_atoms.groupby("residue_id")[["x", "y", "z"]]
            selected_atoms.loc[:, "distance"] = chain_group.diff(-1).apply(
                get_distance, axis=1
            )
            del chain_group
        else:
            selected_atoms.loc[:, "distance"] = (
                selected_atoms[["x", "y", "z"]].diff(-1).apply(get_distance, axis=1)
            )
        selected_atoms.loc[:, "atom1"] = atom1
        selected_atoms.loc[:, "atom2"] = atom2

        selected_atoms = selected_atoms.dropna(how="any")
        selected_atoms = selected_atoms[
            selected_atoms["distance"] > distance * upper_bound_coefficient
        ]

        selected_atoms = selected_atoms.drop_duplicates(keep="first")
        output.append(selected_atoms)
    output = pd.concat(output, ignore_index=True)
    if output.empty:
        return None
    return output


def residue_continuity_check(chain: MiniChain) -> Optional[pd.DataFrame]:
    """Check for residue continuity in a chain.

    :param chain: Chain to check for continuity
    :return: DataFrame of discontinuous residues if found, None otherwise
    """
    start, end = chain["residue_id"].min(), chain["residue_id"].max() + 1
    reference = set(range(int(start), int(end)))
    missing = reference - set(chain["residue_id"].unique()) - {0}
    missing = np.array(list(missing), dtype="int")
    missing.sort()
    missing = missing.astype("str")
    if len(missing) == 0:
        return None
    sample = chain.iloc[[0]][["model_id", "chain_id"]]
    sample["missed_residues"] = ",".join(missing)
    sample["residue_id"] = 0
    sample["chain_length"] = len(reference) - len(missing)
    return sample


def residue_completeness_check(chain: MiniChain) -> Optional[pd.DataFrame]:
    """Check for residue completeness in a chain.

    :param chain: Chain to check for completeness
    :return: DataFrame of incomplete residues if found, None otherwise
    """
    count_dict = chain.value_counts("residue_id", sort=True)
    residue_id = [k for k, v in count_dict.items() if v != 3]
    if len(residue_id) > 0:
        res = chain[chain["residue_id"].isin(residue_id)]
        return res
    return None


def standard_residue_check(
    chain: MiniChain, label_length: int = 3
) -> Optional[pd.DataFrame]:
    """Check for standard residues in a chain.

    :param chain: Chain to check for standard residues
    :param label_length: Length of residue labels to check
    :return: DataFrame of non-standard residues if found, None otherwise
    """
    if label_length == 3:
        standard_residues = basic_amino_acid_20
    else:
        standard_residues = basic_amino_acid_20_s

    non_standard = chain[~chain["residue_label"].isin(standard_residues)]
    if non_standard.empty:
        return None
    return non_standard.iloc[:, :-4]


def residue_angel_check(
    chain: MiniChain, origin: str = "CA", points: Tuple[str, str] = ("C", "N")
) -> Optional[pd.DataFrame]:
    """Check residue angles in a chain.

    :param chain: Chain to check for angles
    :param origin: Origin atom for angle calculation, defaults to "CA"
    :param points: Points to calculate angle between, defaults to ("C", "N")
    :return: DataFrame of residues with invalid angles if found, None otherwise
    """
    dic = {"N": [4, 7], "CA": [7, 10], "C": [10, 13], "O": [13, None]}
    residue_table = chain.pivot(
        index=["model_id", "residue_label", "chain_id", "residue_id"],
        columns=["atom"],
        values=["x", "y", "z"],
    ).reset_index()
    residue_table.columns = [
        "_".join(col) if col[1] else col[0] for col in residue_table.columns
    ]
    residue_table = (
        residue_table.sort_values("residue_id")
        .iloc[:, [0, 1, 2, 3, 6, 10, 14, 5, 9, 13, 4, 8, 12, 7, 11, 15]]
        .reset_index(drop=True)
    )
    residue_data = residue_table.to_numpy()
    origin_indices = dic[origin]
    point1_indices = dic[points[0]]
    point2_indices = dic[points[1]]
    v_1 = (
        residue_data[:, point1_indices[0] : point1_indices[1]]
        - residue_data[:, origin_indices[0] : origin_indices[1]]
    )
    v_2 = (
        residue_data[:, point2_indices[0] : point2_indices[1]]
        - residue_data[:, origin_indices[0] : origin_indices[1]]
    )
    angles = [get_angle(v_1[i], v_2[i]) for i in range(len(v_1))]
    residue_table = residue_table[["model_id", "chain_id"]].drop_duplicates()
    residue_table["min_angle_NAC"] = min(angles)
    residue_table["max_angle_NAC"] = max(angles)

    return residue_table


@on_entry()
def residue_count(chain: MiniChain) -> pd.DataFrame:
    """Count residues in a chain.

    :param chain: Chain to count residues in
    :return: DataFrame with residue counts
    """

    chain = chain.drop_duplicates(subset=["residue_id"])
    residue_count_dic = chain.value_counts("residue_label")
    res = chain.iloc[0][["model_id", "chain_id"]].to_dict()
    res.update(dict(residue_count_dic))
    res = {k: [v] for k, v in res.items()}
    del chain
    output = pd.DataFrame(res)
    output["chain_length"] = output.iloc[:, 2:].sum(axis=1)
    return output


def field_check(entry_id: str) -> Union[Entry, str]:
    """Check if an entry has all required fields.

    :param entry_id: PDB ID to check
    :return: Entry object if valid, error message if not
    """
    try:
        entry = Entry(entry_id)
        if not entry.peptide:
            return {"pdb_id": entry_id, "type": "Non-protein"}
        return entry
    except KeyError:
        return {"pdb_id": entry_id, "type": "Non-protein"}


def main_chain_check(file_dict: Dict[str, Any]) -> bool:
    """Check if a file contains main chain information.

    :param file_dict: Dictionary containing file information
    :return: True if file contains main chain info, False otherwise
    """
    if "entity_poly" not in file_dict:
        return False
    if "pdbx_struct_mod_residue" not in file_dict:
        return True
    if not file_dict["pdbx_struct_mod_residue"]:
        return True
    return False


@on_entry()
def residue_contiguity_wills_check(chain: MiniChain) -> Optional[pd.DataFrame]:
    """Check residue contiguity using Will's method.

    :param chain: Chain to check for contiguity
    :return: DataFrame of discontinuous residues if found, None otherwise
    """
    residue_numbers = chain["residue_id"].unique()
    all_consecutive = np.all(residue_numbers[1:] - residue_numbers[:-1] == 1)
    rear = [False, *(~(residue_numbers[1:] - residue_numbers[:-1] == 1))]
    front = [*(~(residue_numbers[1:] - residue_numbers[:-1] == 1)), False]
    res = list(residue_numbers[front]) + list(residue_numbers[rear])
    res = sorted(res)
    if not all_consecutive:
        return chain[chain["residue_id"].isin(res)]
    return None


def sort_dirty_residues(
    dirty_set: List[Optional[pd.DataFrame]],
) -> Optional[pd.DataFrame]:
    """Sort and combine dirty residue sets.

    :param dirty_set: List of DataFrames containing dirty residues
    :return: Combined and sorted DataFrame of dirty residues
    """
    dirty_set = [df for df in dirty_set if df is not None]
    if not dirty_set:
        return None
    dirty = pd.concat(dirty_set, ignore_index=True)
    if "chain_length" in dirty.columns:
        dirty = dirty[
            [
                "model_id",
                "chain_id",
                "residue_id",
                "residue_label",
                "type",
                "chain_length",
            ]
        ]
    else:
        dirty = dirty[["model_id", "chain_id", "residue_id", "residue_label", "type"]]
        dirty["chain_length"] = np.nan
    dirty = dirty.drop_duplicates(subset=["model_id", "chain_id", "residue_id"])
    dirty.sort_values("residue_id", inplace=True)
    return dirty
