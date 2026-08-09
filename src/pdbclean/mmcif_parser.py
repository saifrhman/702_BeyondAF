"""Parse coordinate mmCIF files into chain-level structural observations."""

from __future__ import annotations

from dataclasses import dataclass, field


class MMCIFParseError(RuntimeError):
    """Raised when a coordinate mmCIF cannot be parsed safely."""


@dataclass(frozen=True)
class AtomObservation:
    """One observed atom from the mmCIF `_atom_site` category."""

    model_id: int
    label_chain_id: str
    auth_chain_id: str | None
    entity_id: str | None

    label_seq_id: int | None
    auth_seq_id: str | None

    residue_name: str
    atom_name: str
    alt_id: str | None

    occupancy: float | None

    x: float
    y: float
    z: float

    # Preserve the exact mmCIF token because BRI v1.2.2 distinguishes
    # "." from "?" during Protocol 3.2 disorder checking.
    occupancy_raw: str | None = None


@dataclass
class ChainObservation:
    """All parsed atom observations belonging to one structural chain."""

    pdb_id: str
    model_id: int
    label_chain_id: str

    auth_chain_id: str | None = None
    entity_id: str | None = None
    polymer_type: str | None = None

    atoms: list[AtomObservation] = field(default_factory=list)

    @property
    def canonical_key(self) -> tuple[str, int, str]:
        return (
            self.pdb_id,
            self.model_id,
            self.label_chain_id,
        )

    @property
    def atom_count(self) -> int:
        return len(self.atoms)


import gzip

import gemmi


_MISSING_CIF_VALUES = {"", ".", "?"}


def _optional_text(value: str) -> str | None:
    """Normalize an optional textual mmCIF value."""

    value = value.strip()

    if value in _MISSING_CIF_VALUES:
        return None

    if (
        len(value) >= 2
        and value[0] == value[-1]
        and value[0] in {"'", '"'}
    ):
        value = value[1:-1]

    return value


def _optional_int(value: str) -> int | None:
    value = value.strip()

    if value in _MISSING_CIF_VALUES:
        return None

    try:
        return int(value)
    except ValueError as exc:
        raise MMCIFParseError(
            f"Expected integer mmCIF value, found {value!r}"
        ) from exc


def _required_float(value: str, field: str) -> float:
    value = value.strip()

    if value in _MISSING_CIF_VALUES:
        raise MMCIFParseError(
            f"Required coordinate field {field} is missing"
        )

    try:
        return float(value)
    except ValueError as exc:
        raise MMCIFParseError(
            f"Invalid floating-point value for {field}: {value!r}"
        ) from exc


def _optional_float(value: str) -> float | None:
    value = value.strip()

    if value in _MISSING_CIF_VALUES:
        return None

    try:
        return float(value)
    except ValueError as exc:
        raise MMCIFParseError(
            f"Invalid floating-point mmCIF value: {value!r}"
        ) from exc


def _polymer_type_map(
    block: gemmi.cif.Block,
) -> dict[str, str]:
    """Map entity IDs to `_entity_poly.type`."""

    entity_ids = list(
        block.find_values("_entity_poly.entity_id")
    )
    polymer_types = list(
        block.find_values("_entity_poly.type")
    )

    if len(entity_ids) != len(polymer_types):
        raise MMCIFParseError(
            "_entity_poly entity/type columns have different lengths"
        )

    result: dict[str, str] = {}

    for entity_id, polymer_type in zip(
        entity_ids,
        polymer_types,
    ):
        clean_entity_id = _optional_text(entity_id)
        clean_polymer_type = _optional_text(polymer_type)

        if clean_entity_id is None or clean_polymer_type is None:
            continue

        result[clean_entity_id] = clean_polymer_type

    return result


def parse_coordinate_mmcif_bytes(
    compressed_bytes: bytes,
    *,
    pdb_id: str,
) -> list[ChainObservation]:
    """Parse one gzipped coordinate mmCIF into chain observations.

    All `_atom_site` rows are parsed. Scientific acceptance/rejection is
    deliberately left to the quality-filter stage.
    """

    try:
        raw_bytes = gzip.decompress(compressed_bytes)
    except (OSError, EOFError) as exc:
        raise MMCIFParseError(
            f"{pdb_id}: invalid or truncated gzip stream"
        ) from exc

    try:
        document = gemmi.cif.read_string(
            raw_bytes.decode("utf-8")
        )
    except (UnicodeDecodeError, RuntimeError) as exc:
        raise MMCIFParseError(
            f"{pdb_id}: invalid mmCIF content"
        ) from exc

    if len(document) == 0:
        raise MMCIFParseError(
            f"{pdb_id}: mmCIF contains no data block"
        )

    block = document.sole_block()

    required_tags = [
        "_atom_site.label_asym_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_atom_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
    ]

    for tag in required_tags:
        if len(block.find_values(tag)) == 0:
            raise MMCIFParseError(
                f"{pdb_id}: required atom-site column missing: {tag}"
            )

    tags = [
        "_atom_site.pdbx_PDB_model_num",
        "_atom_site.label_asym_id",
        "_atom_site.auth_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.auth_seq_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.occupancy",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
    ]

    columns: dict[str, list[str]] = {}

    atom_count: int | None = None

    for tag in tags:
        values = list(block.find_values(tag))

        if not values:
            # Optional columns are represented as missing values for every row.
            columns[tag] = []
            continue

        if atom_count is None:
            atom_count = len(values)
        elif len(values) != atom_count:
            raise MMCIFParseError(
                f"{pdb_id}: atom-site column length mismatch for {tag}"
            )

        columns[tag] = values

    if atom_count is None or atom_count == 0:
        raise MMCIFParseError(
            f"{pdb_id}: atom-site table contains no atom rows"
        )

    def value(tag: str, row: int, default: str = ".") -> str:
        values = columns[tag]
        return values[row] if values else default

    polymer_types = _polymer_type_map(block)

    chains: dict[
        tuple[int, str],
        ChainObservation,
    ] = {}

    for row in range(atom_count):
        model_text = value(
            "_atom_site.pdbx_PDB_model_num",
            row,
            "1",
        )

        model_id = _optional_int(model_text)

        if model_id is None:
            model_id = 1

        label_chain_id = value(
            "_atom_site.label_asym_id",
            row,
        ).strip()

        if label_chain_id in _MISSING_CIF_VALUES:
            raise MMCIFParseError(
                f"{pdb_id}: atom row {row + 1} has no label_asym_id"
            )

        entity_id = _optional_text(
            value("_atom_site.label_entity_id", row)
        )

        chain_key = (model_id, label_chain_id)

        if chain_key not in chains:
            chains[chain_key] = ChainObservation(
                pdb_id=pdb_id.lower(),
                model_id=model_id,
                label_chain_id=label_chain_id,
                auth_chain_id=_optional_text(
                    value("_atom_site.auth_asym_id", row)
                ),
                entity_id=entity_id,
                polymer_type=(
                    polymer_types.get(entity_id)
                    if entity_id is not None
                    else None
                ),
            )

        chain = chains[chain_key]

        atom = AtomObservation(
            model_id=model_id,
            label_chain_id=label_chain_id,
            auth_chain_id=_optional_text(
                value("_atom_site.auth_asym_id", row)
            ),
            entity_id=entity_id,
            label_seq_id=_optional_int(
                value("_atom_site.label_seq_id", row)
            ),
            auth_seq_id=_optional_text(
                value("_atom_site.auth_seq_id", row)
            ),
            residue_name=value(
                "_atom_site.label_comp_id",
                row,
            ).strip(),
            atom_name=value(
                "_atom_site.label_atom_id",
                row,
            ).strip(),
            alt_id=_optional_text(
                value("_atom_site.label_alt_id", row)
            ),
            occupancy=_optional_float(
                value("_atom_site.occupancy", row)
            ),
            x=_required_float(
                value("_atom_site.Cartn_x", row),
                "Cartn_x",
            ),
            y=_required_float(
                value("_atom_site.Cartn_y", row),
                "Cartn_y",
            ),
            z=_required_float(
                value("_atom_site.Cartn_z", row),
                "Cartn_z",
            ),
            occupancy_raw=value(
                "_atom_site.occupancy",
                row,
            ).strip(),
        )

        chain.atoms.append(atom)

    return sorted(
        chains.values(),
        key=lambda chain: (
            chain.model_id,
            chain.label_chain_id,
        ),
    )
