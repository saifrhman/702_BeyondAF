"""Discovery of mmCIF objects from versioned PDB S3 snapshots."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import re
from datetime import datetime, timezone
from typing import Iterator
from urllib.parse import quote, urlencode
from urllib.request import urlopen
from xml.etree import ElementTree as ET

import gemmi


S3_NAMESPACE = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}


class SnapshotError(RuntimeError):
    """Raised when a PDB snapshot cannot be discovered safely."""


@dataclass(frozen=True)
class SnapshotObject:
    """Metadata for one compressed mmCIF object."""

    snapshot: str
    pdb_id: str
    s3_key: str
    size_bytes: int
    etag: str
    last_modified_utc: datetime


def _parse_page(
    xml_bytes: bytes,
    *,
    snapshot: str,
    source_prefix: str,
) -> tuple[list[SnapshotObject], bool, str | None]:
    """Parse one S3 ListObjectsV2 response."""

    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError as exc:
        raise SnapshotError(f"Invalid S3 XML response: {exc}") from exc

    objects: list[SnapshotObject] = []

    for item in root.findall("s3:Contents", S3_NAMESPACE):
        key = item.findtext("s3:Key", namespaces=S3_NAMESPACE)

        if not key or not key.endswith(".cif.gz"):
            continue

        if not key.startswith(source_prefix):
            raise SnapshotError(
                f"S3 returned object outside requested prefix: {key}"
            )

        filename = key.rsplit("/", 1)[-1]
        pdb_id = filename.removesuffix(".cif.gz").lower()

        size_text = item.findtext(
            "s3:Size",
            namespaces=S3_NAMESPACE,
        )
        etag_text = item.findtext(
            "s3:ETag",
            namespaces=S3_NAMESPACE,
        )
        modified_text = item.findtext(
            "s3:LastModified",
            namespaces=S3_NAMESPACE,
        )

        if size_text is None or etag_text is None or modified_text is None:
            raise SnapshotError(
                f"Incomplete metadata returned for {key}"
            )

        try:
            size_bytes = int(size_text)
        except ValueError as exc:
            raise SnapshotError(
                f"Invalid object size for {key}: {size_text!r}"
            ) from exc

        try:
            last_modified = datetime.fromisoformat(
                modified_text.replace("Z", "+00:00")
            ).astimezone(timezone.utc)
        except ValueError as exc:
            raise SnapshotError(
                f"Invalid LastModified timestamp for {key}: "
                f"{modified_text!r}"
            ) from exc

        objects.append(
            SnapshotObject(
                snapshot=snapshot,
                pdb_id=pdb_id,
                s3_key=key,
                size_bytes=size_bytes,
                etag=etag_text.strip('"'),
                last_modified_utc=last_modified,
            )
        )

    is_truncated = (
        root.findtext(
            "s3:IsTruncated",
            default="false",
            namespaces=S3_NAMESPACE,
        ).lower()
        == "true"
    )

    next_token = root.findtext(
        "s3:NextContinuationToken",
        namespaces=S3_NAMESPACE,
    )

    if is_truncated and not next_token:
        raise SnapshotError(
            "S3 response is truncated but contains no continuation token"
        )

    return objects, is_truncated, next_token


def iter_snapshot_objects(
    *,
    bucket_url: str,
    snapshot: str,
    source_prefix: str,
    page_size: int = 1000,
    timeout_seconds: int = 60,
) -> Iterator[SnapshotObject]:
    """Yield all compressed mmCIF objects from a fixed snapshot."""

    if not 1 <= page_size <= 1000:
        raise ValueError("page_size must be between 1 and 1000")

    expected_prefix = (
        f"{snapshot}/pub/pdb/data/structures/divided/mmCIF/"
    )

    if source_prefix != expected_prefix:
        raise SnapshotError(
            "source_prefix does not match the canonical divided mmCIF "
            f"path for snapshot {snapshot}"
        )

    base_url = bucket_url.rstrip("/") + "/"
    continuation_token: str | None = None

    while True:
        params = {
            "list-type": "2",
            "prefix": source_prefix,
            "max-keys": page_size,
        }

        if continuation_token is not None:
            params["continuation-token"] = continuation_token

        url = base_url + "?" + urlencode(params)

        try:
            with urlopen(url, timeout=timeout_seconds) as response:
                xml_bytes = response.read()
        except OSError as exc:
            raise SnapshotError(
                f"Failed to list snapshot objects: {exc}"
            ) from exc

        objects, truncated, continuation_token = _parse_page(
            xml_bytes,
            snapshot=snapshot,
            source_prefix=source_prefix,
        )

        yield from objects

        if not truncated:
            break


@dataclass(frozen=True)
class ResolvedSnapshot:
    """A snapshot that has been resolved and verified against S3."""

    selection_mode: str
    snapshot_id: str
    layout: str
    source_prefix: str
    sample_mmcif_key: str


def canonical_mmcif_prefix(snapshot_id: str) -> str:
    """Return the canonical divided-mmCIF prefix for a snapshot."""

    if (
        not isinstance(snapshot_id, str)
        or len(snapshot_id) != 8
        or not snapshot_id.isdigit()
    ):
        raise SnapshotError(
            f"Invalid snapshot ID {snapshot_id!r}; expected YYYYMMDD"
        )

    return (
        f"{snapshot_id}/pub/pdb/data/structures/"
        "divided/mmCIF/"
    )


def discover_snapshot_ids(
    *,
    bucket_url: str,
    timeout_seconds: int = 60,
) -> list[str]:
    """Discover dated top-level prefixes in the snapshot bucket."""

    base_url = bucket_url.rstrip("/") + "/"
    continuation_token: str | None = None
    snapshot_ids: set[str] = set()

    while True:
        params = {
            "list-type": "2",
            "delimiter": "/",
            "max-keys": 1000,
        }

        if continuation_token is not None:
            params["continuation-token"] = continuation_token

        url = base_url + "?" + urlencode(params)

        try:
            with urlopen(url, timeout=timeout_seconds) as response:
                xml_bytes = response.read()
        except OSError as exc:
            raise SnapshotError(
                f"Failed to discover snapshot prefixes: {exc}"
            ) from exc

        try:
            root = ET.fromstring(xml_bytes)
        except ET.ParseError as exc:
            raise SnapshotError(
                f"Invalid S3 XML response: {exc}"
            ) from exc

        for item in root.findall(
            "s3:CommonPrefixes",
            S3_NAMESPACE,
        ):
            prefix = item.findtext(
                "s3:Prefix",
                namespaces=S3_NAMESPACE,
            )

            if not prefix:
                continue

            candidate = prefix.rstrip("/")

            if len(candidate) == 8 and candidate.isdigit():
                snapshot_ids.add(candidate)

        truncated = (
            root.findtext(
                "s3:IsTruncated",
                default="false",
                namespaces=S3_NAMESPACE,
            ).lower()
            == "true"
        )

        if not truncated:
            break

        continuation_token = root.findtext(
            "s3:NextContinuationToken",
            namespaces=S3_NAMESPACE,
        )

        if not continuation_token:
            raise SnapshotError(
                "Snapshot-prefix listing was truncated but no "
                "continuation token was returned"
            )

    return sorted(snapshot_ids, reverse=True)



COORDINATE_ATOM_SITE_TAGS = (
    "_atom_site.group_PDB",
    "_atom_site.label_atom_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_asym_id",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
)


def is_coordinate_mmcif_bytes(
    compressed_bytes: bytes,
) -> bool:
    """Return True when a gzipped CIF contains usable atom coordinates."""

    try:
        raw_bytes = gzip.decompress(compressed_bytes)
    except (OSError, EOFError):
        return False

    try:
        document = gemmi.cif.read_string(
            raw_bytes.decode("utf-8")
        )
    except (UnicodeDecodeError, RuntimeError):
        return False

    if len(document) == 0:
        return False

    for block in document:
        for tag in COORDINATE_ATOM_SITE_TAGS:
            if len(block.find_values(tag)) == 0:
                break
        else:
            if len(block.find_values("_atom_site.Cartn_x")) > 0:
                return True

    return False


def download_s3_object_bytes(
    *,
    bucket_url: str,
    s3_key: str,
    timeout_seconds: int = 60,
) -> bytes:
    """Download one public S3 object by exact key."""

    object_url = (
        bucket_url.rstrip("/")
        + "/"
        + quote(s3_key, safe="/")
    )

    try:
        with urlopen(
            object_url,
            timeout=timeout_seconds,
        ) as response:
            return response.read()
    except OSError as exc:
        raise SnapshotError(
            f"Failed to download S3 object {s3_key}: {exc}"
        ) from exc


def download_verified_s3_object_bytes(
    *,
    bucket_url: str,
    s3_key: str,
    expected_size_bytes: int,
    expected_etag: str,
    timeout_seconds: int = 60,
) -> bytes:
    """Download one manifest object and verify immutable source metadata."""

    if expected_size_bytes <= 0:
        raise SnapshotError(
            "Expected S3 object size must be positive"
        )

    expected_etag = expected_etag.strip().strip('"')

    if not expected_etag:
        raise SnapshotError(
            "Expected S3 object ETag must be non-empty"
        )

    object_url = (
        bucket_url.rstrip("/")
        + "/"
        + quote(s3_key, safe="/")
    )

    try:
        with urlopen(
            object_url,
            timeout=timeout_seconds,
        ) as response:
            compressed_bytes = response.read()
            response_etag = response.headers.get("ETag")
    except OSError as exc:
        raise SnapshotError(
            f"Failed to download S3 object {s3_key}: {exc}"
        ) from exc

    actual_size_bytes = len(compressed_bytes)

    if actual_size_bytes != expected_size_bytes:
        raise SnapshotError(
            f"S3 object size mismatch for {s3_key}: "
            f"expected {expected_size_bytes}, "
            f"received {actual_size_bytes}"
        )

    if response_etag is None:
        raise SnapshotError(
            f"S3 object response has no ETag for {s3_key}"
        )

    actual_etag = response_etag.strip().strip('"')

    if actual_etag != expected_etag:
        raise SnapshotError(
            f"S3 object ETag mismatch for {s3_key}: "
            f"expected {expected_etag!r}, "
            f"received {actual_etag!r}"
        )

    return compressed_bytes


def object_is_coordinate_mmcif(
    *,
    bucket_url: str,
    s3_key: str,
    timeout_seconds: int = 60,
) -> bool:
    """Classify a `.cif.gz` object by its actual CIF contents."""

    if not s3_key.lower().endswith(".cif.gz"):
        return False

    compressed_bytes = download_s3_object_bytes(
        bucket_url=bucket_url,
        s3_key=s3_key,
        timeout_seconds=timeout_seconds,
    )

    return is_coordinate_mmcif_bytes(compressed_bytes)


def find_sample_mmcif(
    *,
    bucket_url: str,
    snapshot_id: str,
    timeout_seconds: int = 60,
    max_candidates: int = 20,
) -> SnapshotObject | None:
    """Find a genuine coordinate mmCIF in a snapshot.

    Object listing is recursive beneath the canonical divided-mmCIF prefix.
    Candidate `.cif.gz` files are inspected by content so validation,
    reflection, and map-coefficient CIF files cannot qualify accidentally.
    """

    source_prefix = canonical_mmcif_prefix(snapshot_id)

    candidates = iter_snapshot_objects(
        bucket_url=bucket_url,
        snapshot=snapshot_id,
        source_prefix=source_prefix,
        page_size=1000,
        timeout_seconds=timeout_seconds,
    )

    inspected = 0

    for candidate in candidates:
        if inspected >= max_candidates:
            break

        inspected += 1

        if object_is_coordinate_mmcif(
            bucket_url=bucket_url,
            s3_key=candidate.s3_key,
            timeout_seconds=timeout_seconds,
        ):
            return candidate

    return None


def resolve_snapshot(
    snapshot_config: dict,
    *,
    timeout_seconds: int = 60,
) -> ResolvedSnapshot:
    """Resolve and verify a fixed or automatically selected snapshot.

    Resolution first checks the standard wwPDB divided-mmCIF archive.
    If that is absent, the dated snapshot is searched recursively for
    genuine coordinate mmCIF files.

    This function verifies coordinate availability. Full-snapshot
    completeness is validated later when the complete manifest is built.
    """

    mode = snapshot_config.get("mode")
    bucket_url = snapshot_config.get("bucket_url")

    if not isinstance(bucket_url, str) or not bucket_url:
        raise SnapshotError(
            "snapshot.bucket_url must be configured"
        )

    def resolve_one(
        snapshot_id: str,
        selection_mode: str,
    ) -> ResolvedSnapshot | None:
        canonical_prefix = canonical_mmcif_prefix(snapshot_id)

        # Fast path: standard wwPDB coordinate archive.
        sample = find_sample_mmcif(
            bucket_url=bucket_url,
            snapshot_id=snapshot_id,
            timeout_seconds=timeout_seconds,
        )

        if sample is not None:
            return ResolvedSnapshot(
                selection_mode=selection_mode,
                snapshot_id=snapshot_id,
                layout="canonical_divided_mmcif",
                source_prefix=canonical_prefix,
                sample_mmcif_key=sample.s3_key,
            )

        # Fallback: recursively inspect the entire dated snapshot.
        recursive_sample = find_coordinate_mmcif_recursive(
            bucket_url=bucket_url,
            snapshot_id=snapshot_id,
            timeout_seconds=timeout_seconds,
        )

        if recursive_sample is not None:
            return ResolvedSnapshot(
                selection_mode=selection_mode,
                snapshot_id=snapshot_id,
                layout="recursive_coordinate_files",
                source_prefix=f"{snapshot_id}/",
                sample_mmcif_key=recursive_sample,
            )

        return None

    if mode == "fixed":
        snapshot_id = snapshot_config.get("snapshot_id")

        # Also validates YYYYMMDD format.
        canonical_mmcif_prefix(snapshot_id)

        resolved = resolve_one(
            snapshot_id,
            "fixed",
        )

        if resolved is None:
            raise SnapshotError(
                f"Snapshot {snapshot_id} contains no verified "
                "coordinate mmCIF files"
            )

        return resolved

    if mode == "latest_complete":
        snapshot_ids = discover_snapshot_ids(
            bucket_url=bucket_url,
            timeout_seconds=timeout_seconds,
        )

        if not snapshot_ids:
            raise SnapshotError(
                "No dated snapshot prefixes were discovered"
            )

        for snapshot_id in snapshot_ids:
            resolved = resolve_one(
                snapshot_id,
                "latest_complete",
            )

            if resolved is not None:
                return resolved

        raise SnapshotError(
            "Dated prefixes were found, but none contained "
            "verified coordinate mmCIF files"
        )

    raise SnapshotError(
        "snapshot.mode must be 'fixed' or 'latest_complete'"
    )


COORDINATE_FILENAME_PATTERN = re.compile(
    r"^[a-z0-9]{4}\.cif\.gz$",
    re.IGNORECASE,
)


def iter_s3_objects_recursive(
    *,
    bucket_url: str,
    prefix: str,
    timeout_seconds: int = 60,
    page_size: int = 1000,
) -> Iterator[tuple[str, int]]:
    """Recursively yield all S3 objects beneath an arbitrary prefix.

    No delimiter is used, so objects inside any depth of nested
    pseudo-directories are returned.
    """

    if not 1 <= page_size <= 1000:
        raise ValueError("page_size must be between 1 and 1000")

    base_url = bucket_url.rstrip("/") + "/"
    continuation_token: str | None = None

    while True:
        params = {
            "list-type": "2",
            "prefix": prefix,
            "max-keys": page_size,
        }

        if continuation_token is not None:
            params["continuation-token"] = continuation_token

        url = base_url + "?" + urlencode(params)

        try:
            with urlopen(url, timeout=timeout_seconds) as response:
                xml_bytes = response.read()
        except OSError as exc:
            raise SnapshotError(
                f"Failed to recursively list {prefix}: {exc}"
            ) from exc

        try:
            root = ET.fromstring(xml_bytes)
        except ET.ParseError as exc:
            raise SnapshotError(
                f"Invalid S3 XML response: {exc}"
            ) from exc

        for item in root.findall("s3:Contents", S3_NAMESPACE):
            key = item.findtext(
                "s3:Key",
                namespaces=S3_NAMESPACE,
            )
            size_text = item.findtext(
                "s3:Size",
                namespaces=S3_NAMESPACE,
            )

            if key is None or size_text is None:
                continue

            try:
                size_bytes = int(size_text)
            except ValueError as exc:
                raise SnapshotError(
                    f"Invalid object size for {key}: {size_text!r}"
                ) from exc

            yield key, size_bytes

        truncated = (
            root.findtext(
                "s3:IsTruncated",
                default="false",
                namespaces=S3_NAMESPACE,
            ).lower()
            == "true"
        )

        if not truncated:
            break

        continuation_token = root.findtext(
            "s3:NextContinuationToken",
            namespaces=S3_NAMESPACE,
        )

        if not continuation_token:
            raise SnapshotError(
                "Recursive S3 listing was truncated but no "
                "continuation token was returned"
            )


def find_coordinate_mmcif_recursive(
    *,
    bucket_url: str,
    snapshot_id: str,
    timeout_seconds: int = 60,
    max_coordinate_candidates: int = 25,
) -> str | None:
    """Search an entire dated snapshot recursively for coordinate mmCIF.

    This is a fallback for snapshots whose coordinate archive is not stored
    beneath the usual divided/mmCIF path.

    Files such as validation map-coefficient CIFs are ignored because their
    filenames do not match the canonical four-character PDB coordinate-file
    pattern. Matching candidates are then verified by CIF contents.
    """

    snapshot_prefix = f"{snapshot_id}/"
    inspected_candidates = 0

    for key, _ in iter_s3_objects_recursive(
        bucket_url=bucket_url,
        prefix=snapshot_prefix,
        timeout_seconds=timeout_seconds,
    ):
        filename = key.rsplit("/", 1)[-1]

        if not COORDINATE_FILENAME_PATTERN.fullmatch(filename):
            continue

        inspected_candidates += 1

        if object_is_coordinate_mmcif(
            bucket_url=bucket_url,
            s3_key=key,
            timeout_seconds=timeout_seconds,
        ):
            return key

        if inspected_candidates >= max_coordinate_candidates:
            break

    return None


def iter_resolved_snapshot_objects(
    *,
    bucket_url: str,
    resolved: ResolvedSnapshot,
    timeout_seconds: int = 60,
) -> Iterator[SnapshotObject]:
    """Yield coordinate objects according to the resolved snapshot layout."""

    if resolved.layout == "canonical_divided_mmcif":
        yield from iter_snapshot_objects(
            bucket_url=bucket_url,
            snapshot=resolved.snapshot_id,
            source_prefix=resolved.source_prefix,
            page_size=1000,
            timeout_seconds=timeout_seconds,
        )
        return

    if resolved.layout != "recursive_coordinate_files":
        raise SnapshotError(
            f"Unsupported resolved snapshot layout: {resolved.layout!r}"
        )

    base_url = bucket_url.rstrip("/") + "/"
    continuation_token: str | None = None
    prefix = f"{resolved.snapshot_id}/"

    while True:
        params = {
            "list-type": "2",
            "prefix": prefix,
            "max-keys": 1000,
        }

        if continuation_token is not None:
            params["continuation-token"] = continuation_token

        url = base_url + "?" + urlencode(params)

        try:
            with urlopen(url, timeout=timeout_seconds) as response:
                root = ET.fromstring(response.read())
        except (OSError, ET.ParseError) as exc:
            raise SnapshotError(
                f"Failed to enumerate recursive snapshot "
                f"{resolved.snapshot_id}: {exc}"
            ) from exc

        for item in root.findall("s3:Contents", S3_NAMESPACE):
            key = item.findtext("s3:Key", namespaces=S3_NAMESPACE)
            size_text = item.findtext("s3:Size", namespaces=S3_NAMESPACE)
            etag_text = item.findtext("s3:ETag", namespaces=S3_NAMESPACE)
            modified_text = item.findtext(
                "s3:LastModified",
                namespaces=S3_NAMESPACE,
            )

            if (
                key is None
                or size_text is None
                or etag_text is None
                or modified_text is None
            ):
                continue

            filename = key.rsplit("/", 1)[-1]

            if not COORDINATE_FILENAME_PATTERN.fullmatch(filename):
                continue

            if not object_is_coordinate_mmcif(
                bucket_url=bucket_url,
                s3_key=key,
                timeout_seconds=timeout_seconds,
            ):
                continue

            pdb_id = filename.removesuffix(".cif.gz").lower()

            try:
                size_bytes = int(size_text)
                last_modified = datetime.fromisoformat(
                    modified_text.replace("Z", "+00:00")
                ).astimezone(timezone.utc)
            except ValueError as exc:
                raise SnapshotError(
                    f"Invalid metadata for {key}: {exc}"
                ) from exc

            yield SnapshotObject(
                snapshot=resolved.snapshot_id,
                pdb_id=pdb_id,
                s3_key=key,
                size_bytes=size_bytes,
                etag=etag_text.strip('"'),
                last_modified_utc=last_modified,
            )

        truncated = (
            root.findtext(
                "s3:IsTruncated",
                default="false",
                namespaces=S3_NAMESPACE,
            ).lower()
            == "true"
        )

        if not truncated:
            break

        continuation_token = root.findtext(
            "s3:NextContinuationToken",
            namespaces=S3_NAMESPACE,
        )

        if not continuation_token:
            raise SnapshotError(
                "Recursive snapshot listing was truncated but "
                "returned no continuation token"
            )
