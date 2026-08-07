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
