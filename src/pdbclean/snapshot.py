"""Discovery of mmCIF objects from versioned PDB S3 snapshots."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Iterator
from urllib.parse import urlencode
from urllib.request import urlopen
from xml.etree import ElementTree as ET


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
