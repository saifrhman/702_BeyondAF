"""Tests for explicit PDBClean Arrow schemas."""

import pyarrow as pa

from pdbclean.schemas import (
    SOURCE_MANIFEST_SCHEMA,
    SOURCE_MANIFEST_SCHEMA_VERSION,
)


def test_source_manifest_schema_identity() -> None:
    assert SOURCE_MANIFEST_SCHEMA_VERSION == "1.1"
    assert (
        SOURCE_MANIFEST_SCHEMA.metadata[b"schema_name"]
        == b"pdbclean_source_manifest"
    )
    assert (
        SOURCE_MANIFEST_SCHEMA.metadata[b"schema_version"]
        == b"1.1"
    )


def test_source_manifest_schema_fields() -> None:
    expected = {
        "snapshot": pa.string(),
        "source_layout": pa.string(),
        "pdb_id": pa.string(),
        "s3_key": pa.string(),
        "size_bytes": pa.int64(),
        "etag": pa.string(),
        "last_modified_utc": pa.timestamp("us", tz="UTC"),
        "manifest_generated_at_utc": pa.timestamp("us", tz="UTC"),
    }

    assert SOURCE_MANIFEST_SCHEMA.names == list(expected)

    for field in SOURCE_MANIFEST_SCHEMA:
        assert field.type == expected[field.name]
        assert field.nullable is False
