"""Snapshot discovery and selection for a pipeline run.

The generic pipeline default is the latest complete snapshot.  Any concrete
snapshot may be selected instead -- interactively, from the CLI, from a
configuration file, or from the UI -- and the frozen COMP702 result is
reproduced by explicitly selecting ``2026-01-01``.

Whatever the selection route, the snapshot is resolved to a concrete immutable
identity (``YYYYMMDD`` plus the verified S3 layout and a sample coordinate
object) *before* any processing begins, and that identity is what is written to
provenance.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from typing import Callable

from pdbclean.runconfig import (
    ResolvedRunConfig,
    RunConfigError,
    with_resolved_snapshot,
)
from pdbclean.snapshot import (
    ResolvedSnapshot,
    SnapshotError,
    discover_snapshot_ids,
    resolve_snapshot,
)


LATEST_SENTINELS = frozenset({"", "latest", "latest_complete", "default"})

_DATE_WITH_DASHES = re.compile(r"^(\d{4})-(\d{2})-(\d{2})$")
_COMPACT_DATE = re.compile(r"^\d{8}$")


class SnapshotSelectionError(ValueError):
    """Raised when a requested snapshot cannot be interpreted."""


def normalise_snapshot_id(value: str) -> str:
    """Accept ``YYYY-MM-DD`` or ``YYYYMMDD`` and return ``YYYYMMDD``."""

    if not isinstance(value, str):
        raise SnapshotSelectionError(
            f"Snapshot identifier must be a string, got {type(value).__name__}"
        )

    candidate = value.strip()

    dashed = _DATE_WITH_DASHES.match(candidate)

    if dashed is not None:
        candidate = "".join(dashed.groups())

    if not _COMPACT_DATE.match(candidate):
        raise SnapshotSelectionError(
            f"Snapshot identifier must be YYYY-MM-DD or YYYYMMDD, "
            f"got {value!r}"
        )

    return candidate


def format_snapshot_id(snapshot_id: str) -> str:
    """Return ``YYYYMMDD`` formatted for humans as ``YYYY-MM-DD``."""

    if not _COMPACT_DATE.match(snapshot_id or ""):
        return snapshot_id

    return f"{snapshot_id[0:4]}-{snapshot_id[4:6]}-{snapshot_id[6:8]}"


@dataclass(frozen=True)
class SnapshotChoice:
    """One entry in the snapshot picker."""

    index: int
    snapshot_id: str
    display: str
    is_latest: bool


def list_available_snapshots(
    *,
    bucket_url: str,
    limit: int | None = None,
    timeout_seconds: int = 60,
    discover: Callable[..., list[str]] | None = None,
) -> list[SnapshotChoice]:
    """Discover snapshots in the bucket, newest first."""

    discover_fn = discover if discover is not None else discover_snapshot_ids

    try:
        snapshot_ids = discover_fn(
            bucket_url=bucket_url,
            timeout_seconds=timeout_seconds,
        )
    except SnapshotError as exc:
        raise SnapshotSelectionError(
            f"Could not list snapshots from {bucket_url}: {exc}"
        ) from exc

    ordered = sorted(
        {s for s in snapshot_ids if _COMPACT_DATE.match(s or "")},
        reverse=True,
    )

    if limit is not None and limit > 0:
        ordered = ordered[:limit]

    return [
        SnapshotChoice(
            index=position + 1,
            snapshot_id=snapshot_id,
            display=format_snapshot_id(snapshot_id),
            is_latest=(position == 0),
        )
        for position, snapshot_id in enumerate(ordered)
    ]


def interpret_selection(
    selection: str | None,
    choices: list[SnapshotChoice],
) -> str | None:
    """Interpret a picker response.

    Returns ``None`` for "use the latest complete snapshot", otherwise the
    concrete ``YYYYMMDD`` identifier.
    """

    if selection is None:
        return None

    candidate = selection.strip()

    if candidate.lower() in LATEST_SENTINELS:
        return None

    if candidate.isdigit() and len(candidate) != 8:
        position = int(candidate)

        for choice in choices:
            if choice.index == position:
                return choice.snapshot_id

        raise SnapshotSelectionError(
            f"No snapshot at list position {position}"
        )

    return normalise_snapshot_id(candidate)


def render_snapshot_menu(choices: list[SnapshotChoice]) -> str:
    """Render the interactive picker exactly as the CLI displays it."""

    lines = ["Available PDB snapshots:", ""]

    if not choices:
        lines.append("  (no dated snapshots discovered)")
        return "\n".join(lines)

    lines.append("  1. latest complete snapshot [default]")

    for choice in choices:
        marker = "  <- latest discovered" if choice.is_latest else ""
        lines.append(
            f"  {choice.index + 1}. {choice.display}{marker}"
        )

    return "\n".join(lines)


def interpret_menu_response(
    response: str | None,
    choices: list[SnapshotChoice],
) -> str | None:
    """Interpret a response to :func:`render_snapshot_menu`.

    Position 1 is always "latest complete snapshot"; positions 2..n map onto
    the discovered list.  An explicit date is also accepted.
    """

    if response is None:
        return None

    candidate = response.strip()

    if candidate.lower() in LATEST_SENTINELS:
        return None

    if candidate.isdigit() and len(candidate) != 8:
        position = int(candidate)

        if position == 1:
            return None

        for choice in choices:
            if choice.index == position - 1:
                return choice.snapshot_id

        raise SnapshotSelectionError(
            f"No snapshot at menu position {position}"
        )

    return normalise_snapshot_id(candidate)


def resolve_snapshot_for_run(
    resolved: ResolvedRunConfig,
    *,
    requested: str | None = None,
    timeout_seconds: int = 60,
    resolver: Callable[..., ResolvedSnapshot] | None = None,
) -> tuple[ResolvedRunConfig, ResolvedSnapshot]:
    """Pin a run configuration to a concrete, verified snapshot identity.

    ``requested`` (CLI ``--snapshot`` or the UI picker) takes precedence over
    ``snapshot.snapshot_id`` already present in the configuration, which in turn
    takes precedence over the ``latest_complete`` default.
    """

    resolve_fn = resolver if resolver is not None else resolve_snapshot

    bucket_url = resolved.get("snapshot.bucket_url")

    if not bucket_url:
        raise RunConfigError("snapshot.bucket_url must be configured")

    snapshot_id: str | None

    if requested is not None:
        snapshot_id = normalise_snapshot_id(requested)
    else:
        configured = resolved.get("snapshot.snapshot_id")
        mode = resolved.get("snapshot.mode")

        if configured:
            snapshot_id = normalise_snapshot_id(str(configured))
        elif mode in {"fixed", "explicit"}:
            raise SnapshotSelectionError(
                "snapshot.mode requests an explicit snapshot but "
                "snapshot.snapshot_id is not set"
            )
        else:
            snapshot_id = None

    if snapshot_id is None:
        snapshot_config = {
            "mode": "latest_complete",
            "bucket_url": bucket_url,
        }
    else:
        snapshot_config = {
            "mode": "fixed",
            "snapshot_id": snapshot_id,
            "bucket_url": bucket_url,
        }

    try:
        snapshot = resolve_fn(
            snapshot_config,
            timeout_seconds=timeout_seconds,
        )
    except SnapshotError as exc:
        raise SnapshotSelectionError(
            f"Could not resolve snapshot: {exc}"
        ) from exc

    pinned = with_resolved_snapshot(
        resolved,
        snapshot_id=snapshot.snapshot_id,
        layout=snapshot.layout,
        source_prefix=snapshot.source_prefix,
        sample_mmcif_key=snapshot.sample_mmcif_key,
        selection_mode=snapshot.selection_mode,
    )

    return pinned, snapshot


def snapshot_provenance(snapshot: ResolvedSnapshot) -> dict[str, str]:
    """Return the snapshot identity block written into run provenance."""

    payload = asdict(snapshot)
    payload["display"] = format_snapshot_id(snapshot.snapshot_id)

    return payload
