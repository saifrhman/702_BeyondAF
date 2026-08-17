"""Snapshot-consistent metadata extraction for downstream duplicate analysis.

This module extracts metadata only.  It deliberately performs no
resolution, PanDDA, same-deposition, virus, ribosome, or other
downstream filtering.
"""

from __future__ import annotations

from dataclasses import dataclass
import gzip

import gemmi


class DownstreamMetadataError(RuntimeError):
    """Raised when snapshot mmCIF metadata cannot be parsed safely."""


_MISSING = {"", ".", "?"}


def _clean_text(value: str) -> str | None:
    value = str(value).strip()

    if value in _MISSING:
        return None

    if (
        len(value) >= 2
        and value[0] == value[-1]
        and value[0] in {"'", '"'}
    ):
        value = value[1:-1]

    return value


def _texts(
    block: gemmi.cif.Block,
    tag: str,
) -> tuple[str, ...]:
    result: list[str] = []

    for value in block.find_values(tag):
        cleaned = _clean_text(value)

        if cleaned is not None:
            result.append(cleaned)

    return tuple(result)


def _floats(
    block: gemmi.cif.Block,
    tag: str,
) -> tuple[float, ...]:
    result: list[float] = []

    for value in block.find_values(tag):
        cleaned = _clean_text(value)

        if cleaned is None:
            continue

        try:
            parsed = float(cleaned)
        except ValueError as exc:
            raise DownstreamMetadataError(
                f"Invalid floating-point metadata "
                f"{tag}={cleaned!r}"
            ) from exc

        result.append(parsed)

    return tuple(result)


def _first_text(
    block: gemmi.cif.Block,
    tag: str,
) -> str | None:
    values = _texts(block, tag)

    if not values:
        return None

    if len(set(values)) > 1:
        raise DownstreamMetadataError(
            f"Expected at most one distinct value for {tag}; "
            f"found {values!r}"
        )

    return values[0]


@dataclass(frozen=True)
class EntryMetadata:
    """Unfiltered metadata extracted from one frozen snapshot mmCIF."""

    pdb_id: str

    experimental_methods: tuple[str, ...]
    refine_ls_d_res_high: tuple[float, ...]
    em_3d_reconstruction_resolution: tuple[float, ...]

    initial_deposition_date: str | None
    struct_title: str | None
    struct_keywords_text: tuple[str, ...]

    deposit_group_ids: tuple[str, ...]
    deposit_group_titles: tuple[str, ...]
    deposit_group_descriptions: tuple[str, ...]
    deposit_group_types: tuple[str, ...]

    has_deposit_group: bool
    deposit_group_mentions_pandda: bool
    entry_mentions_pandda: bool


def parse_entry_metadata_bytes(
    compressed_bytes: bytes,
    *,
    pdb_id: str,
) -> EntryMetadata:
    """Parse downstream metadata from one gzipped coordinate mmCIF."""

    try:
        raw = gzip.decompress(compressed_bytes)
    except (OSError, EOFError) as exc:
        raise DownstreamMetadataError(
            f"{pdb_id}: invalid gzip stream"
        ) from exc

    try:
        document = gemmi.cif.read_string(
            raw.decode("utf-8")
        )
    except (UnicodeDecodeError, RuntimeError) as exc:
        raise DownstreamMetadataError(
            f"{pdb_id}: invalid mmCIF"
        ) from exc

    if len(document) == 0:
        raise DownstreamMetadataError(
            f"{pdb_id}: mmCIF has no data block"
        )

    block = document.sole_block()

    methods = _texts(
        block,
        "_exptl.method",
    )

    refine = _floats(
        block,
        "_refine.ls_d_res_high",
    )

    em_resolution = _floats(
        block,
        "_em_3d_reconstruction.resolution",
    )

    deposition_date = _first_text(
        block,
        "_pdbx_database_status.recvd_initial_deposition_date",
    )

    title = _first_text(
        block,
        "_struct.title",
    )

    keywords = _texts(
        block,
        "_struct_keywords.text",
    )

    group_ids = _texts(
        block,
        "_pdbx_deposit_group.group_id",
    )
    group_titles = _texts(
        block,
        "_pdbx_deposit_group.group_title",
    )
    group_descriptions = _texts(
        block,
        "_pdbx_deposit_group.group_description",
    )
    group_types = _texts(
        block,
        "_pdbx_deposit_group.group_type",
    )

    group_text = "\n".join(
        (
            *group_ids,
            *group_titles,
            *group_descriptions,
            *group_types,
        )
    ).lower()

    entry_text = "\n".join(
        (
            title or "",
            *keywords,
            group_text,
        )
    ).lower()

    has_deposit_group = bool(
        group_ids
        or group_titles
        or group_descriptions
        or group_types
    )

    return EntryMetadata(
        pdb_id=pdb_id.lower(),
        experimental_methods=methods,
        refine_ls_d_res_high=refine,
        em_3d_reconstruction_resolution=em_resolution,
        initial_deposition_date=deposition_date,
        struct_title=title,
        struct_keywords_text=keywords,
        deposit_group_ids=group_ids,
        deposit_group_titles=group_titles,
        deposit_group_descriptions=group_descriptions,
        deposit_group_types=group_types,
        has_deposit_group=has_deposit_group,
        deposit_group_mentions_pandda=(
            "pandda" in group_text
        ),
        entry_mentions_pandda=(
            "pandda" in entry_text
        ),
    )
