"""Tests for fixed PDB snapshot discovery."""

from datetime import timezone

import pytest

from pdbclean.snapshot import (
    SnapshotError,
    _parse_page,
)


SNAPSHOT = "20260101"
PREFIX = "20260101/pub/pdb/data/structures/divided/mmCIF/"


def test_parse_valid_s3_page() -> None:
    xml = b"""<?xml version="1.0" encoding="UTF-8"?>
<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <IsTruncated>false</IsTruncated>
  <Contents>
    <Key>20260101/pub/pdb/data/structures/divided/mmCIF/00/100d.cif.gz</Key>
    <LastModified>2026-01-01T12:00:00.000Z</LastModified>
    <ETag>"abc123"</ETag>
    <Size>23548</Size>
  </Contents>
</ListBucketResult>
"""

    objects, truncated, token = _parse_page(
        xml,
        snapshot=SNAPSHOT,
        source_prefix=PREFIX,
    )

    assert len(objects) == 1

    obj = objects[0]

    assert obj.snapshot == SNAPSHOT
    assert obj.pdb_id == "100d"
    assert obj.size_bytes == 23548
    assert obj.etag == "abc123"
    assert obj.last_modified_utc.tzinfo == timezone.utc
    assert truncated is False
    assert token is None


def test_non_mmcif_objects_are_ignored() -> None:
    xml = b"""<?xml version="1.0" encoding="UTF-8"?>
<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <IsTruncated>false</IsTruncated>
  <Contents>
    <Key>20260101/pub/pdb/data/structures/divided/mmCIF/00/readme.txt</Key>
    <LastModified>2026-01-01T12:00:00.000Z</LastModified>
    <ETag>"abc123"</ETag>
    <Size>100</Size>
  </Contents>
</ListBucketResult>
"""

    objects, _, _ = _parse_page(
        xml,
        snapshot=SNAPSHOT,
        source_prefix=PREFIX,
    )

    assert objects == []


def test_object_outside_prefix_is_rejected() -> None:
    xml = b"""<?xml version="1.0" encoding="UTF-8"?>
<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <IsTruncated>false</IsTruncated>
  <Contents>
    <Key>20250101/pub/pdb/data/structures/divided/mmCIF/00/100d.cif.gz</Key>
    <LastModified>2026-01-01T12:00:00.000Z</LastModified>
    <ETag>"abc123"</ETag>
    <Size>23548</Size>
  </Contents>
</ListBucketResult>
"""

    with pytest.raises(
        SnapshotError,
        match="outside requested prefix",
    ):
        _parse_page(
            xml,
            snapshot=SNAPSHOT,
            source_prefix=PREFIX,
        )


def test_truncated_page_requires_continuation_token() -> None:
    xml = b"""<?xml version="1.0" encoding="UTF-8"?>
<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <IsTruncated>true</IsTruncated>
</ListBucketResult>
"""

    with pytest.raises(
        SnapshotError,
        match="no continuation token",
    ):
        _parse_page(
            xml,
            snapshot=SNAPSHOT,
            source_prefix=PREFIX,
        )


def test_resolve_fixed_prefers_canonical_layout(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from datetime import datetime, timezone

    from pdbclean.snapshot import (
        SnapshotObject,
        resolve_snapshot,
    )

    sample = SnapshotObject(
        snapshot="20260101",
        pdb_id="100d",
        s3_key=(
            "20260101/pub/pdb/data/structures/"
            "divided/mmCIF/00/100d.cif.gz"
        ),
        size_bytes=23548,
        etag="abc",
        last_modified_utc=datetime.now(timezone.utc),
    )

    monkeypatch.setattr(
        "pdbclean.snapshot.find_sample_mmcif",
        lambda **kwargs: sample,
    )

    resolved = resolve_snapshot(
        {
            "mode": "fixed",
            "snapshot_id": "20260101",
            "bucket_url": "https://example.invalid",
        }
    )

    assert resolved.snapshot_id == "20260101"
    assert resolved.layout == "canonical_divided_mmcif"
    assert resolved.sample_mmcif_key.endswith("100d.cif.gz")


def test_resolve_fixed_uses_recursive_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from pdbclean.snapshot import resolve_snapshot

    monkeypatch.setattr(
        "pdbclean.snapshot.find_sample_mmcif",
        lambda **kwargs: None,
    )

    monkeypatch.setattr(
        "pdbclean.snapshot.find_coordinate_mmcif_recursive",
        lambda **kwargs: (
            "20260101/some/new/layout/100d.cif.gz"
        ),
    )

    resolved = resolve_snapshot(
        {
            "mode": "fixed",
            "snapshot_id": "20260101",
            "bucket_url": "https://example.invalid",
        }
    )

    assert resolved.snapshot_id == "20260101"
    assert resolved.layout == "recursive_coordinate_files"
    assert resolved.source_prefix == "20260101/"


def test_latest_complete_skips_newer_invalid_snapshot(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from datetime import datetime, timezone

    from pdbclean.snapshot import (
        SnapshotObject,
        resolve_snapshot,
    )

    monkeypatch.setattr(
        "pdbclean.snapshot.discover_snapshot_ids",
        lambda **kwargs: ["20260415", "20260101"],
    )

    def fake_canonical(**kwargs):
        snapshot_id = kwargs["snapshot_id"]

        if snapshot_id == "20260415":
            return None

        return SnapshotObject(
            snapshot="20260101",
            pdb_id="100d",
            s3_key=(
                "20260101/pub/pdb/data/structures/"
                "divided/mmCIF/00/100d.cif.gz"
            ),
            size_bytes=23548,
            etag="abc",
            last_modified_utc=datetime.now(timezone.utc),
        )

    monkeypatch.setattr(
        "pdbclean.snapshot.find_sample_mmcif",
        fake_canonical,
    )

    monkeypatch.setattr(
        "pdbclean.snapshot.find_coordinate_mmcif_recursive",
        lambda **kwargs: None,
    )

    resolved = resolve_snapshot(
        {
            "mode": "latest_complete",
            "bucket_url": "https://example.invalid",
        }
    )

    assert resolved.snapshot_id == "20260101"
    assert resolved.layout == "canonical_divided_mmcif"


class _FakeS3Response:
    """Minimal urlopen-compatible response for verified-download tests."""

    def __init__(
        self,
        payload: bytes,
        *,
        etag: str | None,
    ) -> None:
        self._payload = payload
        self.headers = {}

        if etag is not None:
            self.headers["ETag"] = etag

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback) -> None:
        return None

    def read(self) -> bytes:
        return self._payload


def test_download_verified_s3_object_bytes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from pdbclean.snapshot import download_verified_s3_object_bytes

    payload = b"synthetic-cif-gzip"

    monkeypatch.setattr(
        "pdbclean.snapshot.urlopen",
        lambda *args, **kwargs: _FakeS3Response(
            payload,
            etag='"abc123"',
        ),
    )

    result = download_verified_s3_object_bytes(
        bucket_url="https://example.invalid",
        s3_key="20260101/path/1abc.cif.gz",
        expected_size_bytes=len(payload),
        expected_etag="abc123",
        timeout_seconds=5,
    )

    assert result == payload


def test_download_verified_s3_object_rejects_size_mismatch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from pdbclean.snapshot import download_verified_s3_object_bytes

    payload = b"abc"

    monkeypatch.setattr(
        "pdbclean.snapshot.urlopen",
        lambda *args, **kwargs: _FakeS3Response(
            payload,
            etag='"etag123"',
        ),
    )

    with pytest.raises(
        SnapshotError,
        match="S3 object size mismatch",
    ):
        download_verified_s3_object_bytes(
            bucket_url="https://example.invalid",
            s3_key="20260101/path/1abc.cif.gz",
            expected_size_bytes=999,
            expected_etag="etag123",
        )


def test_download_verified_s3_object_rejects_etag_mismatch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from pdbclean.snapshot import download_verified_s3_object_bytes

    payload = b"abc"

    monkeypatch.setattr(
        "pdbclean.snapshot.urlopen",
        lambda *args, **kwargs: _FakeS3Response(
            payload,
            etag='"actual-etag"',
        ),
    )

    with pytest.raises(
        SnapshotError,
        match="S3 object ETag mismatch",
    ):
        download_verified_s3_object_bytes(
            bucket_url="https://example.invalid",
            s3_key="20260101/path/1abc.cif.gz",
            expected_size_bytes=len(payload),
            expected_etag='"expected-etag"',
        )


def test_download_verified_s3_object_requires_response_etag(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from pdbclean.snapshot import download_verified_s3_object_bytes

    payload = b"abc"

    monkeypatch.setattr(
        "pdbclean.snapshot.urlopen",
        lambda *args, **kwargs: _FakeS3Response(
            payload,
            etag=None,
        ),
    )

    with pytest.raises(
        SnapshotError,
        match="response has no ETag",
    ):
        download_verified_s3_object_bytes(
            bucket_url="https://example.invalid",
            s3_key="20260101/path/1abc.cif.gz",
            expected_size_bytes=len(payload),
            expected_etag="etag123",
        )
