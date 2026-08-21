"""Bronze source-object resolution for one run's snapshot.

This is the layer that answers "which immutable source object did *this run*
parse for this PDB entry?" -- independently of anything that happened to the
chain afterwards.

Why it exists
-------------

Stage-14 removal means "this chain is not retained in the geometrically
deduplicated Gold training population". It does **not** mean the deposited
structure was deleted or lost. Visual inspection of a duplicate pair must
therefore resolve structures from the run's **source** layer, never from the
Gold retained set:

    resolved snapshot
        -> Bronze source manifest        <- this module
            -> immutable source object
                -> deposited mmCIF
                    -> model / chain selection

    (Gold retained/removed status is metadata shown alongside, and never
     controls whether a structure can be loaded.)

Snapshot fidelity
-----------------

Every key in the Bronze manifest is snapshot-scoped, e.g.::

    20260101/pub/pdb/data/structures/divided/mmCIF/ac/7acj.cif.gz

Fetching that key returns the object *as the snapshot froze it*, which is not
the same thing as fetching whatever the PDB serves today. The manifest also
records the object's ``etag`` and ``size_bytes``, so a materialised copy can be
verified byte-identical to the object the pipeline actually parsed. A
downloaded object whose ETag disagrees is rejected rather than displayed.
"""

from __future__ import annotations

import gzip
import shutil
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any


#: Refuse to materialise a single object larger than this. Structure files are
#: megabytes; anything far larger is not something to stream into a cache.
MAX_OBJECT_BYTES = 256 * 1024 * 1024

DEFAULT_TIMEOUT_SECONDS = 60


class SourceIndexError(RuntimeError):
    """Raised when a source object cannot be resolved or verified."""


@dataclass(frozen=True)
class SourceObject:
    """The immutable identity of one deposited entry in one snapshot."""

    pdb_id: str
    snapshot_id: str
    s3_key: str
    etag: str | None = None
    size_bytes: int | None = None
    source_layout: str | None = None

    def url(self, bucket_url: str) -> str:
        return f"{bucket_url.rstrip('/')}/{self.s3_key.lstrip('/')}"

    def to_dict(self) -> dict[str, Any]:
        return {
            "pdb_id": self.pdb_id,
            "snapshot_id": self.snapshot_id,
            "s3_key": self.s3_key,
            "etag": self.etag,
            "size_bytes": self.size_bytes,
            "source_layout": self.source_layout,
        }


class BronzeSourceIndex:
    """Look up source objects in a run's own Bronze manifest.

    Lookups use Arrow filter pushdown, so a single entry is resolved without
    loading the 246k-row manifest into memory.
    """

    def __init__(self, manifest_path: str | Path) -> None:
        self.manifest_path = Path(manifest_path)
        self._cache: dict[str, SourceObject | None] = {}

    @property
    def available(self) -> bool:
        return self.manifest_path.is_file()

    @classmethod
    def for_snapshot(
        cls,
        *,
        output_root: str | Path,
        snapshot_id: str,
    ) -> "BronzeSourceIndex":
        """The Bronze manifest belonging to one snapshot of one output root."""

        return cls(
            Path(output_root)
            / str(snapshot_id)
            / "bronze"
            / "source_manifest.parquet"
        )

    def lookup(self, pdb_id: str) -> SourceObject | None:
        """Return the source object for one PDB id, or ``None``."""

        key = pdb_id.lower()

        if key in self._cache:
            return self._cache[key]

        if not self.available:
            self._cache[key] = None

            return None

        import pyarrow.compute as pc
        import pyarrow.dataset as ds

        dataset = ds.dataset(str(self.manifest_path), format="parquet")

        table = dataset.to_table(
            filter=pc.equal(pc.utf8_lower(ds.field("pdb_id")), key),
            columns=[
                "snapshot",
                "pdb_id",
                "s3_key",
                "etag",
                "size_bytes",
                "source_layout",
            ],
        )

        if table.num_rows == 0:
            self._cache[key] = None

            return None

        row = table.to_pylist()[0]

        found = SourceObject(
            pdb_id=str(row["pdb_id"]).lower(),
            snapshot_id=str(row["snapshot"]),
            s3_key=str(row["s3_key"]),
            etag=(str(row["etag"]) if row.get("etag") else None),
            size_bytes=(
                int(row["size_bytes"])
                if row.get("size_bytes") is not None
                else None
            ),
            source_layout=row.get("source_layout"),
        )

        self._cache[key] = found

        return found


def _normalise_etag(value: str | None) -> str | None:
    return value.strip('"').strip() if value else None


def materialise(
    source: SourceObject,
    *,
    bucket_url: str,
    destination: str | Path,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
    verify: bool = True,
    opener: Any = None,
) -> Path:
    """Fetch one snapshot object into the disposable cache, verified.

    Fetches the **snapshot-scoped** key, so the bytes are the ones the run
    parsed, not whatever the PDB serves today. The response ETag is checked
    against the value the Bronze manifest recorded; a mismatch is an error, not
    a warning, because showing a different revision of a structure than the one
    the science used would be misleading.

    The object is decompressed into the cache so the viewer can serve plain
    mmCIF. This is derived data: deleting it loses nothing.
    """

    target = Path(destination)

    if target.is_file():
        return target

    if source.size_bytes is not None and source.size_bytes > MAX_OBJECT_BYTES:
        raise SourceIndexError(
            f"Source object {source.s3_key} is {source.size_bytes} bytes, "
            f"above the {MAX_OBJECT_BYTES}-byte materialisation limit."
        )

    url = source.url(bucket_url)
    fetch = opener or urllib.request.urlopen

    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")

    try:
        with fetch(url, timeout=timeout_seconds) as response:
            observed = _normalise_etag(
                response.headers.get("ETag")
                if hasattr(response, "headers")
                else None
            )
            expected = _normalise_etag(source.etag)

            if verify and expected and observed and observed != expected:
                raise SourceIndexError(
                    f"Source object {source.s3_key} does not match the ETag "
                    f"recorded for snapshot {source.snapshot_id} "
                    f"(expected {expected}, got {observed}). Refusing to show "
                    "a different revision than the one this run used."
                )

            payload = response.read()
    except urllib.error.URLError as exc:
        raise SourceIndexError(
            f"Could not retrieve {url}: {exc}"
        ) from exc

    try:
        if source.s3_key.endswith(".gz"):
            payload = gzip.decompress(payload)

        partial.write_bytes(payload)
        partial.replace(target)
    finally:
        if partial.exists():
            partial.unlink()

    return target


def cache_path_for(
    source: SourceObject,
    *,
    hot_root: str | Path,
) -> Path:
    """Where a materialised copy of one source object lives.

    Keyed by snapshot, so two snapshots' revisions of the same entry never
    collide.
    """

    return (
        Path(hot_root)
        / str(source.snapshot_id)
        / source.pdb_id[1:3]
        / f"{source.pdb_id}.cif"
    )


# ---------------------------------------------------------------------------
# Chain namespace resolution
# ---------------------------------------------------------------------------
#
# PDBClean's canonical chain identity is ``label_asym_id``. Deposited mmCIF
# files also carry ``auth_asym_id``, and the two frequently differ -- in the
# frozen 2026-01-01 release they differ for 73.5% of removed chains.
#
# The viewer therefore must not assume ``label_asym_id == auth_asym_id``. It
# resolves the pair from the run's own cleaning output, which records both for
# the whole eligible population, and states which identifier it selects on.


@dataclass(frozen=True)
class ChainNames:
    """Both chain identifiers for one chain, as the run recorded them."""

    pdb_id: str
    label_asym_id: str
    auth_asym_id: str | None = None
    model_id: int = 1

    @property
    def diverges(self) -> bool:
        return (
            self.auth_asym_id is not None
            and self.auth_asym_id != self.label_asym_id
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "pdb_id": self.pdb_id,
            "label_asym_id": self.label_asym_id,
            "auth_asym_id": self.auth_asym_id,
            "model_id": self.model_id,
            "diverges": self.diverges,
            # The viewer selects on the canonical identifier, which is what
            # the pipeline computed BRI for.
            "selected_on": "label_asym_id",
        }


class ChainNameIndex:
    """Resolve label/auth chain identifiers from a run's cleaning output.

    Reads the accepted-chain table, which covers the whole eligible
    population, so a chain that Stage 14 later removed is resolved exactly as
    readily as one that was retained.
    """

    def __init__(self, accepted_path: str | Path) -> None:
        self.accepted_path = Path(accepted_path)
        self._cache: dict[tuple[str, str], ChainNames | None] = {}

    @property
    def available(self) -> bool:
        return self.accepted_path.is_file()

    @classmethod
    def for_protocol(
        cls,
        *,
        output_root: str | Path,
        snapshot_id: str,
        protocol: str,
    ) -> "ChainNameIndex":
        return cls(
            Path(output_root)
            / str(snapshot_id)
            / str(protocol)
            / "quality"
            / "merged"
            / "accepted.parquet"
        )

    def lookup(self, pdb_id: str, label_asym_id: str) -> ChainNames | None:
        key = (pdb_id.lower(), label_asym_id)

        if key in self._cache:
            return self._cache[key]

        if not self.available:
            self._cache[key] = None

            return None

        import pyarrow.compute as pc
        import pyarrow.dataset as ds

        dataset = ds.dataset(str(self.accepted_path), format="parquet")

        table = dataset.to_table(
            filter=(
                pc.equal(pc.utf8_lower(ds.field("pdb_id")), key[0])
                & pc.equal(ds.field("label_chain_id"), label_asym_id)
            ),
            columns=["pdb_id", "label_chain_id", "auth_chain_id", "model_id"],
        )

        if table.num_rows == 0:
            self._cache[key] = None

            return None

        row = table.to_pylist()[0]

        found = ChainNames(
            pdb_id=str(row["pdb_id"]).lower(),
            label_asym_id=str(row["label_chain_id"]),
            auth_asym_id=(
                str(row["auth_chain_id"])
                if row.get("auth_chain_id") is not None
                else None
            ),
            model_id=int(row.get("model_id") or 1),
        )

        self._cache[key] = found

        return found
