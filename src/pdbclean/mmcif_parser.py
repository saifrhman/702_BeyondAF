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
