"""Project a source mmCIF chain onto exactly the Gold-retained residues.

OpenFold builds every structural feature from ``mmcif_object.chain_to_seqres``:
``make_mmcif_features`` takes ``num_res = len(chain_to_seqres[chain_id])`` and
``get_atom_coords`` then walks ``range(num_res)`` looking up
``seqres_to_structure[chain_id][i]``.  That sequence is the *deposited*
polymer, read from ``_entity_poly_seq``, and it includes residues Protocol 3.2
removed -- unobserved ones among them.

Our MSAs were searched on the retained sequence.  Feeding OpenFold the
deposited sequence therefore produces structure features of one length and MSA
features of another, and the pipeline does not raise: it simply returns a
feature dict whose rows no longer correspond.  Training on that would silently
learn against misaligned evolutionary signal, which is worse than a crash.

This module fixes the input view and nothing else.  It rebuilds one chain's
``chain_to_seqres`` and ``seqres_to_structure`` so they describe exactly the
residues Gold retained, re-indexed to ``0..L-1``.  Coordinates are untouched --
``ResidueAtPosition`` still points at the same author residue numbers, so
``get_atom_coords`` reads the same atoms it always would.  The model, the loss,
the MSA semantics and the retained-chain identities are all unchanged.

The residue set is not recomputed here.  It is taken from Gold's
``retained_label_seq_ids``, the same lineage
:func:`pdbclean.geometric_validation.reconstruct_retained_backbone_chain` uses,
so there is one definition of "retained" in the project rather than two.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence


@dataclass(frozen=True)
class ProjectionError(Exception):
    """Raised when a chain cannot be projected onto its Gold lineage."""

    message: str

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.message


def seqres_start_number(mmcif_object: Any, entity_id: str) -> int:
    """The mmCIF ``label_seq_id`` that OpenFold's seqres index 0 refers to.

    OpenFold computes ``seq_idx = label_seq_id - min(_entity_poly_seq.num)``
    for the chain's entity.  That minimum is conventionally 1, but convention
    is not a contract and an off-by-one here would shift every residue against
    its coordinates without changing any length, so nothing downstream would
    notice.  It is therefore read from the deposited file.

    ``MmcifObject.raw_string`` is the parsed mmCIF dict, so the loop is still
    available after parsing; ``entity_id`` comes from Gold's own manifest,
    which avoids re-deriving OpenFold's chain-to-entity mapping here.
    """

    parsed_info = mmcif_object.raw_string

    nums = parsed_info.get("_entity_poly_seq.num")
    entities = parsed_info.get("_entity_poly_seq.entity_id")

    if not nums or not entities:
        raise ProjectionError(
            f"{mmcif_object.file_id} has no _entity_poly_seq loop; cannot "
            "establish the seqres numbering origin"
        )

    wanted = str(entity_id)

    matching = [
        int(num) for num, ent in zip(nums, entities) if str(ent) == wanted
    ]

    if not matching:
        raise ProjectionError(
            f"{mmcif_object.file_id} has no _entity_poly_seq rows for entity "
            f"{wanted!r}"
        )

    return min(matching)


def projected_indices(
    retained_label_seq_ids: Sequence[int],
    seq_start_num: int,
) -> list[int]:
    """Gold ``label_seq_id`` values as OpenFold seqres indices."""

    requested = tuple(retained_label_seq_ids)

    if not requested:
        raise ProjectionError("retained_label_seq_ids must not be empty")

    if len(set(requested)) != len(requested):
        raise ProjectionError("retained_label_seq_ids must not contain duplicates")

    return [int(residue_id) - seq_start_num for residue_id in requested]


def project_mmcif_object(
    mmcif_object: Any,
    auth_chain_id: str,
    retained_label_seq_ids: Sequence[int],
    entity_id: str | None = None,
    seq_start_num: int | None = None,
):
    """Return a copy of ``mmcif_object`` restricted to the retained residues.

    Only ``chain_to_seqres`` and ``seqres_to_structure`` change, and only for
    ``auth_chain_id``.  Every other chain is dropped, so a caller cannot
    accidentally build features for a chain Gold did not retain.
    """

    # Imported here so this module stays importable in the pdbclean
    # environment, which has no OpenFold.
    from openfold.data import mmcif_parsing

    if seq_start_num is None:
        if entity_id is None:
            raise ProjectionError(
                "one of entity_id or seq_start_num is required to place Gold "
                "label_seq_ids on OpenFold's seqres index"
            )

        seq_start_num = seqres_start_number(mmcif_object, entity_id)

    seqres = mmcif_object.chain_to_seqres.get(auth_chain_id)

    if seqres is None:
        raise ProjectionError(
            f"auth chain {auth_chain_id!r} absent from {mmcif_object.file_id}"
        )

    structure_map = mmcif_object.seqres_to_structure.get(auth_chain_id)

    if structure_map is None:
        raise ProjectionError(
            f"no seqres_to_structure for {mmcif_object.file_id}:{auth_chain_id}"
        )

    indices = projected_indices(retained_label_seq_ids, seq_start_num)

    out_of_range = [i for i in indices if i < 0 or i >= len(seqres)]

    if out_of_range:
        raise ProjectionError(
            f"{mmcif_object.file_id}:{auth_chain_id} Gold lineage falls "
            f"outside the deposited sequence of length {len(seqres)}: "
            f"indices {out_of_range[:5]}"
        )

    missing = [i for i in indices if i not in structure_map]

    if missing:
        raise ProjectionError(
            f"{mmcif_object.file_id}:{auth_chain_id} has no structure mapping "
            f"for seqres indices {missing[:5]}"
        )

    projected_sequence = "".join(seqres[i] for i in indices)
    projected_map = {new: structure_map[old] for new, old in enumerate(indices)}

    return mmcif_parsing.MmcifObject(
        file_id=mmcif_object.file_id,
        header=mmcif_object.header,
        structure=mmcif_object.structure,
        chain_to_seqres={auth_chain_id: projected_sequence},
        seqres_to_structure={auth_chain_id: projected_map},
        raw_string=mmcif_object.raw_string,
    )


def contiguous_ranges(values: Iterable[int]) -> list[tuple[int, int]]:
    """Collapse sorted residue ids into inclusive ``(start, end)`` runs.

    Terminal trimming produces one run; internal excisions produce more.  The
    representation matters because storing 499,770 explicit residue lists costs
    hundreds of megabytes per dataloader worker, while runs cost almost
    nothing when the lineage is contiguous.
    """

    runs: list[tuple[int, int]] = []

    start = previous = None

    for value in values:
        value = int(value)

        if start is None:
            start = previous = value
            continue

        if value == previous + 1:
            previous = value
            continue

        runs.append((start, previous))
        start = previous = value

    if start is not None:
        runs.append((start, previous))

    return runs


def expand_ranges(runs: Sequence[Sequence[int]]) -> list[int]:
    """Inverse of :func:`contiguous_ranges`."""

    values: list[int] = []

    for start, end in runs:
        values.extend(range(int(start), int(end) + 1))

    return values


def projection_entry(retained_label_seq_ids: Sequence[int]) -> list[list[int]]:
    """The compact on-disk form of one chain's Gold lineage."""

    return [list(run) for run in contiguous_ranges(retained_label_seq_ids)]


def load_projection_index(path) -> Mapping[str, list[list[int]]]:
    import json

    with open(path) as handle:
        return json.load(handle)
