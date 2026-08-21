"""Durable snapshot preservation and disposable hot materialisation.

Two layers, deliberately separated:

**A. Durable snapshot identity.**  What makes a completed run reproducible
after a temporary cache expires.  Content-addressed: an immutable source object
is stored once under its content identity, and each snapshot manifest
*references* it.  An object that is unchanged between two PDB snapshots is
preserved once and referenced twice.

**B. Hot working materialisation.**  A fast, node-local or scratch copy
optimised for parsing throughput.  Disposable by design: it may be deleted at
any time and rebuilt from layer A.

::

    PDB snapshot source (S3)
            |
            v  resolve to a concrete snapshot identity
    immutable Bronze manifest
            |
            v  preserve objects by content identity
    durable object store        objects/<content-id>
                                snapshots/<YYYYMMDD>.manifest.json
            |
            v  materialise for computation (disposable)
    hot cache                   cache/<YYYYMMDD>/...
            |
            v
    PDBClean stages

Nothing in this module moves data on import, and nothing here deletes anything.
It describes, plans and records; promotion is an explicit, separately approved
operation.

Object identity follows the verified provenance the pipeline already uses --
the S3 key, byte size and ETag recorded in the Bronze manifest, plus a content
hash where one has actually been computed.  Filename-only identity is never
sufficient and is never used.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable


SNAPSHOT_STORE_SCHEMA_NAME = "pdbclean_snapshot_store_manifest"
SNAPSHOT_STORE_SCHEMA_VERSION = "1.0"


# ---------------------------------------------------------------------------
# Availability states
# ---------------------------------------------------------------------------
#
# A snapshot is not simply "present" or "absent".  These are the states the
# pipeline and the UI distinguish.

#: Listed in the upstream archive, nothing held locally.
REMOTE_AVAILABLE = "remote_available"

#: Objects preserved in the durable content-addressed store.  A run that
#: reached this state stays reproducible after every cache is deleted.
PRESERVED = "preserved"

#: Materialised in the disposable hot cache, ready for fast parsing.
HOT = "hot"

#: Preserved *and* hot.
MATERIALISED = "materialised"

#: Preserved content verified against the recorded identity.
VERIFIED = "verified"

#: Nothing known about this snapshot locally.
UNKNOWN = "unknown"

AVAILABILITY_STATES: tuple[str, ...] = (
    UNKNOWN,
    REMOTE_AVAILABLE,
    HOT,
    PRESERVED,
    MATERIALISED,
    VERIFIED,
)


class SnapshotStoreError(RuntimeError):
    """Raised when the snapshot store cannot answer a question."""


# ---------------------------------------------------------------------------
# Object identity
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ObjectIdentity:
    """Verified identity of one immutable source object.

    ``content_id`` is what the durable store keys on.  It is derived from a
    real content hash when one has been computed, and otherwise from the
    archive's own verified metadata (key, size, ETag) -- never from the
    filename alone.
    """

    pdb_id: str
    source_key: str
    size_bytes: int
    etag: str | None = None
    content_sha256: str | None = None

    @property
    def content_id(self) -> str:
        """Stable content identity used as the durable store key."""

        if self.content_sha256:
            return f"sha256:{self.content_sha256}"

        if self.etag:
            # S3 ETags are content-derived for single-part uploads and remain
            # the archive's own verified identity for multipart ones.
            etag = self.etag.strip('"')

            return f"etag:{etag}:{self.size_bytes}"

        raise SnapshotStoreError(
            f"Object {self.source_key!r} has neither a content hash nor an "
            "ETag; filename-only identity is not sufficient for preservation."
        )

    @property
    def verified(self) -> bool:
        """Whether a real content hash backs this identity."""

        return self.content_sha256 is not None

    def to_dict(self) -> dict[str, Any]:
        return {
            "pdb_id": self.pdb_id,
            "source_key": self.source_key,
            "size_bytes": self.size_bytes,
            "etag": self.etag,
            "content_sha256": self.content_sha256,
            "content_id": self.content_id,
            "identity_verified": self.verified,
        }


# ---------------------------------------------------------------------------
# Store layout
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SnapshotStoreLayout:
    """Where the durable and hot layers live.

    Both roots are configuration, never hard-coded: a clone on another cluster,
    or a compute node with a different scratch, must be able to point them
    somewhere writable.  See ``storage.durable_snapshot_root`` and
    ``storage.hot_cache_root``.
    """

    durable_root: Path
    hot_root: Path

    @classmethod
    def from_config(
        cls,
        resolved: Any,
        *,
        repo_root: str | Path,
    ) -> "SnapshotStoreLayout":
        root = Path(repo_root)

        def _resolve(dotted: str, fallback: str) -> Path:
            value = resolved.get(dotted) or fallback
            candidate = Path(value)

            return candidate if candidate.is_absolute() else root / candidate

        return cls(
            durable_root=_resolve(
                "storage.durable_snapshot_root", "outputs/snapshot_store"
            ),
            hot_root=_resolve("storage.hot_cache_root", "outputs/snapshot_cache"),
        )

    @property
    def objects_root(self) -> Path:
        return self.durable_root / "objects"

    @property
    def manifests_root(self) -> Path:
        return self.durable_root / "snapshots"

    def object_path(self, content_id: str) -> Path:
        """Return the durable path for one content identity.

        Sharded two levels so a directory never holds hundreds of thousands of
        entries.
        """

        scheme, _, digest = content_id.partition(":")
        flat = digest.replace(":", "_")

        return self.objects_root / scheme / flat[:2] / flat[2:4] / flat

    def manifest_path(self, snapshot_id: str) -> Path:
        return self.manifests_root / f"{snapshot_id}.manifest.json"

    def hot_path(self, snapshot_id: str) -> Path:
        return self.hot_root / snapshot_id


# ---------------------------------------------------------------------------
# Snapshot manifests
# ---------------------------------------------------------------------------


@dataclass
class SnapshotManifest:
    """An immutable record of what one snapshot consisted of.

    Written once and never mutated.  Re-preserving a snapshot writes a new
    manifest only if none exists; an existing manifest is authoritative.
    """

    snapshot_id: str
    bucket_url: str
    source_prefix: str
    objects: list[ObjectIdentity] = field(default_factory=list)
    created_at: str | None = None
    preservation_status: str = REMOTE_AVAILABLE

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_name": SNAPSHOT_STORE_SCHEMA_NAME,
            "schema_version": SNAPSHOT_STORE_SCHEMA_VERSION,
            "snapshot_id": self.snapshot_id,
            "bucket_url": self.bucket_url,
            "source_prefix": self.source_prefix,
            "created_at": self.created_at,
            "preservation_status": self.preservation_status,
            "object_count": len(self.objects),
            "total_bytes": sum(obj.size_bytes for obj in self.objects),
            "objects": [obj.to_dict() for obj in self.objects],
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "SnapshotManifest":
        return cls(
            snapshot_id=payload["snapshot_id"],
            bucket_url=payload.get("bucket_url", ""),
            source_prefix=payload.get("source_prefix", ""),
            created_at=payload.get("created_at"),
            preservation_status=payload.get(
                "preservation_status", REMOTE_AVAILABLE
            ),
            objects=[
                ObjectIdentity(
                    pdb_id=entry["pdb_id"],
                    source_key=entry["source_key"],
                    size_bytes=entry["size_bytes"],
                    etag=entry.get("etag"),
                    content_sha256=entry.get("content_sha256"),
                )
                for entry in payload.get("objects", [])
            ],
        )


def write_snapshot_manifest(
    manifest: SnapshotManifest,
    layout: SnapshotStoreLayout,
    *,
    overwrite: bool = False,
) -> Path:
    """Write a snapshot manifest. Refuses to mutate an existing one."""

    target = layout.manifest_path(manifest.snapshot_id)

    if target.exists() and not overwrite:
        raise SnapshotStoreError(
            f"Snapshot manifest already exists and is immutable: {target}"
        )

    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(manifest.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return target


def read_snapshot_manifest(
    snapshot_id: str,
    layout: SnapshotStoreLayout,
) -> SnapshotManifest | None:
    """Read a preserved snapshot manifest, or ``None`` if not preserved."""

    target = layout.manifest_path(snapshot_id)

    if not target.is_file():
        return None

    return SnapshotManifest.from_dict(
        json.loads(target.read_text(encoding="utf-8"))
    )


# ---------------------------------------------------------------------------
# Status and planning
# ---------------------------------------------------------------------------


def snapshot_status(
    snapshot_id: str,
    layout: SnapshotStoreLayout,
    *,
    remote_available: bool = False,
) -> dict[str, Any]:
    """Report what is held for one snapshot, without touching any data."""

    manifest = read_snapshot_manifest(snapshot_id, layout)
    hot = layout.hot_path(snapshot_id)
    hot_present = hot.is_dir()

    if manifest is not None and hot_present:
        state = MATERIALISED
    elif manifest is not None:
        state = PRESERVED
    elif hot_present:
        state = HOT
    elif remote_available:
        state = REMOTE_AVAILABLE
    else:
        state = UNKNOWN

    payload: dict[str, Any] = {
        "snapshot_id": snapshot_id,
        "availability": state,
        "durable_root": str(layout.durable_root),
        "hot_root": str(layout.hot_root),
        "manifest_path": (
            str(layout.manifest_path(snapshot_id))
            if manifest is not None
            else None
        ),
        "hot_path": str(hot) if hot_present else None,
        "preserved_object_count": (
            len(manifest.objects) if manifest is not None else 0
        ),
        "preserved_bytes": (
            sum(obj.size_bytes for obj in manifest.objects)
            if manifest is not None
            else 0
        ),
        "reproducible_without_cache": manifest is not None,
    }

    return payload


def plan_preservation(
    objects: Iterable[ObjectIdentity],
    layout: SnapshotStoreLayout,
) -> dict[str, Any]:
    """Cost a preservation run without moving a single byte.

    Deduplication is the point: an object already preserved under the same
    content identity -- typically because an earlier snapshot referenced it --
    is referenced again, not copied again.
    """

    listed = list(objects)

    seen: set[str] = set()
    to_transfer: list[ObjectIdentity] = []
    already: list[ObjectIdentity] = []
    duplicates_within = 0

    for obj in listed:
        content_id = obj.content_id

        if content_id in seen:
            duplicates_within += 1
            continue

        seen.add(content_id)

        if layout.object_path(content_id).exists():
            already.append(obj)
        else:
            to_transfer.append(obj)

    return {
        "object_count": len(listed),
        "distinct_object_count": len(seen),
        "duplicate_references_within_snapshot": duplicates_within,
        "already_preserved_count": len(already),
        "already_preserved_bytes": sum(o.size_bytes for o in already),
        "to_transfer_count": len(to_transfer),
        "to_transfer_bytes": sum(o.size_bytes for o in to_transfer),
        "total_bytes_if_copied_naively": sum(o.size_bytes for o in listed),
        "bytes_saved_by_deduplication": (
            sum(o.size_bytes for o in listed)
            - sum(o.size_bytes for o in to_transfer)
        ),
        "destination": str(layout.objects_root),
        "unverified_identity_count": sum(
            1 for o in listed if not o.verified
        ),
    }


def provenance_block(
    snapshot_id: str,
    layout: SnapshotStoreLayout,
    *,
    selection_mode: str,
    remote_available: bool = False,
) -> dict[str, Any]:
    """The snapshot-durability block recorded in run provenance.

    Preservation and cache state are *operational* facts.  They record how
    reproducible a run is, and never alter the scientific snapshot identity.
    """

    status = snapshot_status(
        snapshot_id, layout, remote_available=remote_available
    )

    return {
        "snapshot_id": snapshot_id,
        "selection_mode": selection_mode,
        "availability": status["availability"],
        "reproducible_without_cache": status["reproducible_without_cache"],
        "durable_manifest": status["manifest_path"],
        "hot_cache_path": status["hot_path"],
        "note": (
            "Preservation and cache state describe availability only. The "
            "scientific identity of this run is the resolved snapshot id."
        ),
    }
