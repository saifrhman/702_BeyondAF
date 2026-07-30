# -*- coding = utf-8 -*-
# @Time: 2023/9/22 11:23
# @File: pdbx2df.PY
"""PDB/mmCIF file parser and data structure module.

This module provides classes and functions for parsing and handling protein structure
data from PDB/mmCIF files. It includes:

- Classes for representing protein chains and entries with different levels of detail
- Functions for computing and analyzing backbone rigid invariants
- Visualization tools for generating Backbone Invariant Diagrams (BID) and Barcodes (BIB)
- Utilities for coordinate transformation and perturbation analysis

The module serves as the foundation for the BRI (Backbone Rigid Invariant) analysis
pipeline by providing data structures and basic operations.
"""

import os
from typing import List, Dict, Optional, Any, Callable
import warnings

from biotite.structure.io.pdb import PDBFile
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from .base.base_util import StructureBase, COLUMN_MAP
from .invariant import get_invariant, get_invariant_BTP, BTP_ATTR
from .base.math_base import random_matrix

pd.options.mode.chained_assignment = None
column_map: Dict[str, str] = {
    "group_PDB": "cate",
    "label_comp_id": "residue_label",
    "auth_atom_id": "atom",
    "label_entity_id": "entity_id",
    "label_asym_id": "chain_id",
    "label_seq_id": "residue_id",
    "Cartn_x": "x",
    "Cartn_y": "y",
    "Cartn_z": "z",
    "pdbx_PDB_model_num": "model_id",
    "auth_asym_id": "auth_chain_id",
    "auth_seq_id": "auth_residue_id",
}


class MiniChain(StructureBase):
    """A minimal representation of a protein chain for BRI analysis.

    This class provides essential functionality for:

    - Computing backbone rigid invariants
    - Applying coordinate perturbations
    - Generating visualization plots (BIB, BID)

    The minimal representation includes only the necessary atomic coordinates
    and metadata for BRI computation.
    """

    extract_cols: Dict[str, str] = {
        "cate": "str",
        "residue_label": "str",
        "atom": "str",
        "chain_id": "str",
        "residue_id": "int",
        "x": "float64",
        "y": "float64",
        "z": "float64",
        "occupancy": "str",
        "model_id": "int",
    }

    def __init__(
        self,
        pdb_id: str,
        model_id: int,
        chain_id: str,
        start_residue: int,
        chain_length: int,
        extra_keys: Optional[Dict[str, str]] = None,
        data: Optional[pd.DataFrame] = None,
    ) -> None:
        """Initialize a MiniChain instance.

        :param pdb_id: PDB entry identifier
        :param model_id: Model number in the structure
        :param chain_id: Chain identifier
        :param start_residue: Starting residue number
        :param chain_length: Number of residues to include
        :param extra_keys: Additional metadata keys
        :param data: Pre-loaded coordinate data
        """

        super().__init__(pdb_id, extra_keys, data)

        # self.pdb_id: str = pdb_id  #: PDB entry id
        self._model_id: int = model_id  #: Model id
        self._chain_id: str = chain_id  #: Chain id
        self._start_residue: int = (
            start_residue  #: First PDB labeled residue id to start
        )
        self._chain_length: int = chain_length  #: Number of residues
        self._perturb_radius: float = 0.0  #: Perturbation radius

        if data is None:
            # load mmCIF file as dict
            self._coordinates: pd.DataFrame = self._extract_coordinates(
                model_id, chain_id, start_residue, chain_length
            )  #: Atomic coordinates
        self.__base_coordinates: pd.DataFrame = self._coordinates.loc[:]

    def _extract_coordinates(
        self, model_id: int, chain_id: str, start_index: int, chain_length: int
    ) -> pd.DataFrame:
        """Extract atomic coordinates for the specified chain segment.

        :param model_id: Model number
        :param chain_id: Chain identifier
        :param start_index: Starting residue number
        :param chain_length: Number of residues
        :return: DataFrame containing atomic coordinates
        """
        res = self._coordinates.query(
            f'model_id == {model_id} and chain_id == "{chain_id}" and residue_id >= {start_index} and '
            f"residue_id < {start_index + chain_length}"
        )
        return res.loc[:, list(self.extract_cols.keys())]

    def get_chain_invariant_BTP(self) -> pd.DataFrame:
        """Compute bond-torsion point invariants for the chain.

        :return: DataFrame containing bond-torsion point invariants.
        """

        invariant_bri = get_invariant(self.get_feature(), ext=True)
        invariant_btp = get_invariant_BTP(invariant_bri)
        return invariant_btp.loc[
            :, ["model_id", "chain_id", "residue_id", "residue_label"] + BTP_ATTR
        ]

    def get_chain_invariant(self, angles: bool = False) -> pd.DataFrame:
        """Compute backbone rigid invariants for the chain.

        :param angles: Whether to include bond and torsion angles
        :return: DataFrame containing backbone rigid invariants
        """

        invariant_res = get_invariant(self.get_feature(), ext=angles)
        return invariant_res

    @property
    def invariant(self) -> Optional[pd.DataFrame]:
        """Get the backbone rigid invariant.

        Compute the backbone rigid invariant of a chain.

        :return: DataFrame containing backbone rigid invariant.
        """

        try:
            return self.get_chain_invariant()
            """Backbone-Rigid-Invariant of the chain computed from `self.coordinates`"""
        except Exception:
            warnings.warn(
                f"Chain {self.__repr__()} includes unexpected atoms or misses essential atoms, invariants cannot be computed"
            )
            return None

    @property
    def perturb_radius(self) -> float:
        """Get the current perturbation radius.

        :return: Current perturbation radius
        """
        return self._perturb_radius

    @perturb_radius.setter
    def perturb_radius(self, epsilon: float) -> None:
        """Apply random perturbation to atomic coordinates.

        :param epsilon: Perturbation radius
        """
        columns = ["x", "y", "z"]
        perturbation = random_matrix((len(self._coordinates), 3), epsilon)

        # update coordinates and invariant with set perturbation
        self._coordinates[columns] = self.__base_coordinates[columns] + perturbation
        self._perturb_radius = epsilon

    def generate_BID(self) -> Figure:
        """Generate a Backbone Invariant Diagram.

        Creates a plot showing the values of backbone rigid invariants
        along the protein chain.

        :return: Matplotlib :class:`Figure` containing the BID plot
        """
        tick_step = 10
        colors = ("red", "green", "blue")
        fig, axes = plt.subplots(9, 1, sharex="all", figsize=(16, 8))
        if self.invariant is None:
            return fig
        residue_id = self.invariant["residue_id"]
        bri = self.invariant.loc[
            :, ["x(N)", "y(N)", "z(N)", "x(A)", "y(A)", "z(A)", "x(C)", "y(C)", "z(C)"]
        ]

        left, right = int(residue_id.min()), int(residue_id.max())
        tick_l = ((left + tick_step - 1) // tick_step) * tick_step
        tick_r = (right // tick_step + 1) * tick_step

        for i, col in enumerate(bri.columns):
            axes[i].plot(residue_id, bri[col], marker=".", color=colors[i % 3])

            axes[i].set_xlim(left - 0.5, right + 0.5)
            axes[i].set_xticks(range(tick_l, tick_r, tick_step))
            axes[i].set_ylim(-2, 2)
            axes[i].grid(axis="y", alpha=0.6)

        fig.tight_layout()
        fig.show()
        return fig

    def generate_BIB(self) -> Figure:
        """Generate a Backbone Invariant Barcode.

        Creates a barcode visualization of backbone rigid invariants
        using color coding.

        :return: Matplotlib :class:`Figure` containing the BIB plot
        """
        tick_step = 10
        fig, axes = plt.subplots(3, 1, sharex="all", figsize=(16, 3))
        if self.invariant is None:
            return fig
        residue_id = self.invariant["residue_id"]
        bri = self.invariant[
            ["x(N)", "y(N)", "z(N)", "x(A)", "y(A)", "z(A)", "x(C)", "y(C)", "z(C)"]
        ]
        bri = bri.astype("float")

        bri = (bri + 2) / 4
        bri.iloc[0] = [1] * 9
        left, right = int(residue_id.min()), int(residue_id.max())
        tick_l = ((left + tick_step - 1) // tick_step) * tick_step
        tick_r = (right // tick_step + 1) * tick_step

        for i in range(3):
            data = bri.iloc[:, i * 3 : i * 3 + 3]

            axes[i].get_yaxis().set_visible(False)
            axes[i].set_xticks(range(tick_l, tick_r, tick_step))
            axes[i].set_xlim(left - 0.5, right + 0.5)
            axes[i].imshow(
                data.values.reshape(1, -1, 3),
                aspect="auto",
                extent=(left - 0.5, right + 0.5, 0.5, -0.5),
                norm="linear",
            )

        fig.tight_layout()
        fig.show()
        return fig

    def get_feature(
        self,
        patten: str | list[str] = "features",
        by: Optional[str] = None,
        HETATM: bool = False,
    ) -> Optional[pd.DataFrame]:
        """Extract specific features from the chain.

        :param patten: Feature pattern to extract
        :param by: Grouping criterion, 'atom' for atom-wise features or 'comp' for residue-wise features
        :param HETATM: Whether to include HETATM records, defaults to False
        :return: :class:`pandas.DataFrame` containing requested features
        """
        model: pd.DataFrame = self._coordinates
        if not HETATM:
            model: pd.DataFrame = model.loc[model["cate"] == "ATOM"]
        # ATOM record is used to identify proteins or nucleic acid atoms, and the HETATM record is used to identify
        # atoms in small molecules (or ligand).
        if patten == "features":
            res = model[model["atom"].isin(["C", "CA", "N"])]
        elif by == "atom" and isinstance(patten, str):
            res = model[model["atom"] == patten]
        elif by == "atom" and isinstance(patten, list):
            res = model[model["atom"].isin(patten)]
        elif by == "comp":
            res = model[model["residue_label"] == patten]
        else:
            return None

        return pd.DataFrame(res.reset_index(drop=True))

    @classmethod
    def from_pdb(cls, path: str):
        """Basic method that converts PDB files into standard data

        This method cannot properly handle file that includes multiple models. One should be defined in :class:`Entries` to achieve that.

        :param path: PDB file path
        :return: A structured instance of protein data
        """
        pdb_file_annot_cate_map = {
            "chain_id": "label_asym_id",
            "res_id": "label_seq_id",
            "res_name": "label_comp_id",
            "ins_code": "pdbx_PDB_ins_code",
            "hetero": "hetero",
            "atom_name": "label_atom_id",
            "element": "type_symbol",
            "occupancy": "str",
            "atom_id": "id",
            "altloc_id": "label_alt_id",
            "charge": "pdbx_formal_charge",
        }

        file = PDBFile.read(path)
        pdb_data = file.get_structure(
            extra_fields=["occupancy", "atom_id", "charge"], altloc="all"
        )
        data_dict = {
            cate: pdb_data.__dict__["_annot"][cate]
            for cate in pdb_data.get_annotation_categories()
        }
        df = pd.DataFrame(data_dict)

        # Dataframe alignment
        df["group_PDB"] = df["hetero"].apply(lambda x: "HETATM" if x else "ATOM")
        df = df.rename(columns=pdb_file_annot_cate_map).rename(columns=COLUMN_MAP)
        df[["x", "y", "z"]] = pdb_data.coord[0]

        # fill missing values
        for k, type in cls.extract_cols.items():
            if k not in df.columns:
                if type == "str":
                    df[k] = " "
                elif type.startswith("float"):
                    df[k] = -1.0
                elif type.startswith("int"):
                    df[k] = -1

        start: int = df["residue_id"].min()
        length = len(df["residue_id"].unique())
        return cls(path, -1, "A", start, length, data=df)

    def __len__(self) -> int:
        """Count the number of residues in this chain

        :return: Number of unique residue id
        """

        return len(self.coordinates.residue_id.unique())

    def __repr__(self) -> str:
        """Get string representation of the chain

        :return: String representation
        """
        return f"{self.pdb_id}-{self._model_id}-{self._chain_id}-{self._start_residue}-{self._chain_length}"

    def __eq__(self, value):
        """Compare chain with another value

        :param value: Value to compare with
        :return: True if equal, False otherwise
        """
        return self.__repr__() == value


class Chain(MiniChain):
    """An instance to represent a specified chain. This class enables invariant computation, perturbation, and
    generation of defined BIB, BID plots.

    This class extends MiniChain with additional entity information and functionality:

    - Entity ID tracking
    - Enhanced coordinate extraction
    - Full chain representation

    """

    extract_cols: Dict[str, str] = {
        "cate": "str",
        "residue_label": "str",
        "atom": "str",
        "entity_id": "int",
        "chain_id": "str",
        "auth_chain_id": "str",
        "residue_id": "int",
        "auth_residue_id": "int",
        "x": "float64",
        "y": "float64",
        "z": "float64",
        "occupancy": "str",
        "model_id": "int",
    }

    def __init__(
        self,
        pdb_id: str,
        entity_id: int,
        model_id: int,
        chain_id: str,
        start_residue: int,
        chain_length: int,
        extra_keys: Optional[Dict[str, str]] = None,
        data: Optional[pd.DataFrame] = None,
    ):
        """Initialize a Chain instance.

        :param pdb_id: PDB entry identifier
        :param entity_id: Entity identifier
        :param model_id: Model number in the structure
        :param chain_id: Chain identifier
        :param start_residue: Starting residue number
        :param chain_length: Number of residues to include
        :param extra_keys: Additional metadata keys
        :param data: Pre-loaded coordinate data
        """

        super().__init__(
            pdb_id, model_id, chain_id, start_residue, chain_length, extra_keys, data
        )

        # self.pdb_id: str = pdb_id  #: PDB entry id
        self._entity_id: int = entity_id  #: Entity id
        self._model_id: int = model_id  #: Model id
        self._chain_id: str = chain_id  #: Chain id
        self._start_residue: int = (
            start_residue  #: First PDB labeled residue id to start
        )
        self._chain_length: int = chain_length  #: Number of residues
        # self._perturb_radius: int = 0  #: Perturbation radius

        # load mmCIF file as dict
        self._coordinates: pd.DataFrame = self._extract_coordinates(
            model_id, chain_id, start_residue, chain_length
        )  #: Atomic coordinates
        self.__base_coordinates: pd.DataFrame = self._coordinates.loc[:]

        # TODO consider better override
        self._entity_info: dict[str, Any] = (
            self._entity_info.get(self._entity_id, {})
            if self._entity_info is not None
            else {}
        )  #: Keep this entity only

        if self._entity_info.get("chain_id") is not None and (
            not self._chain_id.isalpha()
        ):
            # self._chain_id = self.entity_info.get("chain_id")
            self._coordinates["chain_id"] = self._entity_info.get("chain_id")

        # self._invariant: pd.DataFrame = self.get_chain_invariant()
        """Backbone-Rigid-Invariant of the chain computed from `self.coordinates`"""

        if data is None:
            del self._file

    def _extract_coordinates(
        self, model_id: int, chain_id: str, start_index: int, chain_length: int
    ) -> pd.DataFrame:
        """Extract atomic coordinates for the specified chain segment.

        :param model_id: Model number
        :param chain_id: Chain identifier
        :param start_index: Starting residue number
        :param chain_length: Number of residues
        :return: :class:`pandas.DataFrame` containing atomic coordinates
        """
        res = self._coordinates.query(
            f'model_id == {model_id} and chain_id == "{chain_id}" and residue_id >= {start_index} and '
            f"residue_id < {start_index + chain_length}"
        )
        return res.loc[:, list(self.extract_cols.keys())]

    def __repr__(self) -> str:
        """Get string representation of the chain.

        :return: String representation
        """
        return f"{self.pdb_id}-{self._entity_id}-{self._model_id}-{self._chain_id}-{self._start_residue}-{self._chain_length}"


class MiniEntry(StructureBase):
    """An instance to represent a PDB entry with minimal info, given a path to a PDB/mmCIF file.

    This class provides functionality for:

    - Loading and parsing PDB/mmCIF files
    - Extracting and managing protein chains
    - Computing invariants for the entire entry
    - Saving invariant data

    """

    extract_cols: Dict[str, str] = {
        "cate": "str",
        "residue_label": "str",
        "atom": "str",
        "chain_id": "str",
        "residue_id": "int",
        "x": "float64",
        "y": "float64",
        "z": "float64",
        "occupancy": "str",
        "model_id": "int",
    }  #: Map of column names and their types to extract from the PDB/mmCIF file

    def __init__(
        self,
        path: str,
        extra_keys: Optional[Dict[str, str]] = None,
        data: Optional[pd.DataFrame] = None,
    ):
        """Initialize a MiniEntry instance.

        :param path: PDB entry ID or path to mmCIF file
        :param extra_keys: Additional metadata keys
        :param data: Pre-loaded coordinate data
        """

        super().__init__(path, extra_keys, data)

        self.model_num: int = 0
        self.add_extra_fields(
            {"label_entity_info": "struct_asym"}
        )  #: PDB labeled entity-chain mapping
        self.peptide: bool = False  #: Entry includes peptide
        if not self.entity_info:
            raise KeyError(
                "Cannot find entity data from the entry input. Please check if this entry is a protein or DNA/RNA molecule."
            )
        for v in self.entity_info.values():
            if v["entity_poly_type"].startswith("polypeptide"):
                self.peptide = True

        self.chains: list[MiniChain] = (
            self._extract_chains()
        )  #: An array of protein chains

        self.__base_coordinates: pd.DataFrame = self._coordinates.loc[:]

    def _extract_chains(self) -> list[MiniChain]:
        """Extract protein chains from the entry.

        :return: Tuple of :class:`MiniChain` instances
        """
        coordinates = self.coordinates
        if coordinates.empty:
            return []

        self.model_num = coordinates.iloc[-1]["model_id"]
        chains: list[MiniChain] = []

        for (model_id, chain_id), data in coordinates.groupby(["model_id", "chain_id"]):
            start_residue = data.residue_id.min()
            chain_length = data.residue_id.max() - start_residue + 1
            chains.append(
                MiniChain(
                    self.pdb_id,
                    model_id,
                    chain_id,
                    start_residue,
                    chain_length,
                    data=data,
                )
            )

        return chains

    def get_entry_invariant(
        self, chain_ids: Optional[List[str]] = None, angles: bool = False
    ) -> pd.DataFrame:
        """Compute backbone-rigid-invariant for the whole entry.

        :param chain_ids: List of chain identifiers to process, defaults to `None` for all chains
        :param angles: Whether to include bond and torsion angles
        :return: :class:`pandas.DataFrame` containing invariants for all specified chains
        """

        invariant_func = on_entry(chain_ids=chain_ids)(get_invariant)

        invariant_res = invariant_func(self, angles)
        return invariant_res

    def save_invariants(
        self,
        target_path: Optional[str] = None,
        chain_ids: Optional[List[str]] = None,
        angles: bool = False,
    ):
        """Compute invariants and save as CSV files

        :param target_path: Directory to save invariant files, defaults to `None`
        :param chain_ids: List of chain identifiers to process, defaults to `None` to save all chains in the entry as a single file
        :param angles: Whether to include bond and torsion angles
        """
        if target_path is not None and not os.path.isdir(target_path):
            os.makedirs(target_path)

        if not chain_ids:  # save all chains as a single file
            if target_path is not None:
                file = os.path.join(target_path, f"{self.pdb_id}.csv")
            else:
                file = f"{self.pdb_id}.csv"

            invariants = self.get_entry_invariant(chain_ids, angles)
            invariants.to_csv(file, index=False)
            return

        for chain_name in chain_ids:
            if target_path is not None:
                file = os.path.join(target_path, f"{chain_name}.csv")
            else:
                file = f"{chain_name}.csv"

            chain_invariant = self[chain_name].get_chain_invariant(angles)
            chain_invariant.to_csv(file, index=False)  # save results

    def get_feature(
        self,
        model_index: int,
        patten: str | list[str] = "features",
        by: Optional[str] = None,
        HETATM: bool = False,
    ) -> Optional[pd.DataFrame]:
        """Extract feature atoms from coordinates.

        :param model_index: Model number
        :param patten: Feature pattern to extract, defaults to `features`
        :param by: Grouping criterion, 'atom' for atom-wise features or 'comp' for residue-wise features, defaults to `None`
        :param HETATM: Whether to include HETATM records, defaults to `False`
        :return: :class:`pandas.DataFrame` containing extracted features
        """

        # check model id
        model_index = model_index - 1 if model_index > 0 else model_index
        model = self.chains[model_index]
        return model.get_feature(patten, by, HETATM)

    @property
    def invariant(self) -> Optional[pd.DataFrame]:
        """Get the backbone rigid invariant.

        Compute the backbone rigid invariant of all loaded chains inside the entry

        :return: DataFrame containing backbone rigid invariant
        """

        try:
            return self.get_entry_invariant()
        except Exception:
            warnings.warn(
                f"BRI computation not available due to some exceptions for atoms in {self.__repr__()}."
            )
            return None

    @property
    def perturb_radius(self) -> float:
        """Get the current perturbation radius.

        :return: Current perturbation radius
        """
        return self._perturb_radius

    @perturb_radius.setter
    def perturb_radius(self, epsilon: float) -> None:
        """Apply random perturbation to atomic coordinates and update the invariant

        :param epsilon: Perturbation radius
        """
        # update coordinates and invariant with set perturbation
        for chain in self.chains:
            chain.perturb_radius = epsilon

        self._coordinates: pd.DataFrame = pd.concat(
            [chain.coordinates for chain in self.chains]
        )
        self._perturb_radius: float = epsilon

    @classmethod
    def from_pdb(cls, path: str):
        """Basic method that converts PDB files into standard data

        This method cannot properly handle file that includes multiple models. One should be defined in :class:`Entries` to achieve that.

        :param path: PDB file path
        :return: A structured instance of protein data
        """
        ...
        pdb_file_annot_cate_map = {
            "chain_id": "label_asym_id",
            "res_id": "label_seq_id",
            "res_name": "label_comp_id",
            "ins_code": "pdbx_PDB_ins_code",
            "hetero": "hetero",
            "atom_name": "label_atom_id",
            "element": "type_symbol",
            "occupancy": "str",
            "atom_id": "id",
            "altloc_id": "label_alt_id",
            "charge": "pdbx_formal_charge",
        }

        file = PDBFile.read(path)
        pdb_data = file.get_structure(
            extra_fields=["occupancy", "atom_id", "charge"], altloc="all"
        )
        data_dict = {
            cate: pdb_data.__dict__["_annot"][cate]
            for cate in pdb_data.get_annotation_categories()
        }
        df = pd.DataFrame(data_dict)

        # Dataframe alignment
        df["group_PDB"] = df["hetero"].apply(lambda x: "HETATM" if x else "ATOM")
        df = df.rename(columns=pdb_file_annot_cate_map).rename(columns=COLUMN_MAP)
        df[["x", "y", "z"]] = pdb_data.coord[0]

        # fill missing values
        for k, type in cls.extract_cols.items():
            if k not in df.columns:
                if type == "str":
                    df[k] = " "
                elif type.startswith("float"):
                    df[k] = -1.0
                elif type.startswith("int"):
                    df[k] = -1

        return cls(path, data=df)

    def __repr__(self):
        """Get string representation of the entry.

        :return: String representation
        """
        return f"{self.pdb_id}"

    def __getitem__(self, key):
        """Get a chain by its identifier.

        :param key: Chain identifier
        :return: Chain instance
        :raises ValueError: If chain not found
        """
        try:
            index = self.chains.index(key)
            return self.chains[index]
        except Exception:
            raise ValueError(f"Chain {key} not found in entry")


class Entry(MiniEntry):
    """An instance to represent a PDB entry with more information, given an entry ID or a path to a PDB/mmCIF file.

    This class extends MiniEntry with additional functionality:

    - Entity information tracking
    - Enhanced chain management
    - Full structure representation

    """

    extract_cols: Dict[str, str] = {
        "cate": "str",
        "residue_label": "str",
        "atom": "str",
        "entity_id": "int",
        "chain_id": "str",
        "auth_chain_id": "str",
        "residue_id": "int",
        "auth_residue_id": "int",
        "x": "float64",
        "y": "float64",
        "z": "float64",
        "occupancy": "str",
        "model_id": "int",
    }  #: Map of column names and their types to extract from the PDB/mmCIF file

    def __init__(
        self,
        path: str,
        extra_keys: Optional[Dict[str, Any]] = None,
        data: Optional[pd.DataFrame] = None,
    ):
        """Initialize an Entry instance.

        :param path: PDB entry ID or path to mmCIF file
        :param extra_keys: Additional metadata keys
        :param data: Pre-loaded coordinate data
        """

        super().__init__(path, extra_keys, data)

        self.add_extra_fields(
            {"label_entity_info": "struct_asym"}
        )  #: PDB labeled entity-chain mapping
        self.peptide: bool = False  #: Entry includes peptide
        if not self.entity_info:
            raise KeyError(
                "Cannot find entity data from the entry input. Please check if this entry is a protein or DNA/RNA molecule."
            )
        for v in self.entity_info.values():
            if v["entity_poly_type"].startswith("polypeptide"):
                self.peptide = True

        self.chains: list[Chain] = self._extract_chains()  #: An array of protein chains

        self.__base_coordinates: pd.DataFrame = self._coordinates[:]

        self._invariant: Optional[pd.DataFrame] = None
        """Backbone-Rigid-Invariant of the chain computed from `self.coordinates`"""

        del self._file

    def _extract_chains(self) -> list[Chain]:
        """Extract protein chains from the entry.

        :return: Tuple of :class:`Chain` instances
        """
        coordinates = self.coordinates
        if coordinates.empty:
            return []

        self.model_num: int = coordinates.iloc[-1]["model_id"]
        chains: list[Chain] = []

        for (model_id, chain_id), data in coordinates.groupby(["model_id", "chain_id"]):
            entity_id = data.iloc[0]["entity_id"]
            start_residue = data.residue_id.min()
            chain_length = data.residue_id.max() - start_residue + 1
            chains.append(
                Chain(
                    self.pdb_id,
                    entity_id,
                    model_id,
                    chain_id,
                    start_residue,
                    chain_length,
                    data=data,
                )
            )

        return chains


def on_entry(HETATM: bool = False, chain_ids: Optional[list[str]] = None) -> Callable:
    """Decorator for applying functions to all chains in an entry.

    This decorator enables functions to operate on all chains in a PDB entry,
    with options for filtering by chain IDs and including HETATM records.

    :param HETATM: Whether to include HETATM records, defaults to `False`
    :param chain_ids: List of chain identifiers to process, defaults to `None` to process all chains
    :return: Decorated function
    """

    def inner(func_on_chain):
        def wrapper(*args, **kwargs):
            pdb_entry: MiniEntry = args[0]
            args = args[1:]
            output_list_on_entry = []
            # start iterating chains
            chains: list[MiniChain] | tuple[MiniChain] = (
                [pdb_entry[i] for i in chain_ids] if chain_ids else pdb_entry.chains
            )
            for i in chains:
                chain = i.get_feature("features", HETATM=HETATM)
                # find no target atoms, go to next chain
                if chain is None or chain.empty:
                    continue
                # extract one single chain
                output_list_on_entry.append(func_on_chain(chain, *args, **kwargs))

            # if no objects in output list, return None
            if len(output_list_on_entry) == 0:
                return None
            # organize output DataFrame
            output_on_entry = pd.concat(
                output_list_on_entry, ignore_index=True
            ).drop_duplicates()
            output_on_entry["pdb_id"] = pdb_entry.pdb_id
            # oder columns
            order = [-1] + list(range(len(output_on_entry.columns) - 1))
            output_on_entry = output_on_entry.iloc[:, order]

            return output_on_entry

        return wrapper

    return inner
