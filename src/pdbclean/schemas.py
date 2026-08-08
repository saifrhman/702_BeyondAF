"""Explicit Apache Arrow schemas used by the PDBClean pipeline."""

import pyarrow as pa


SOURCE_MANIFEST_SCHEMA_VERSION = "1.1"

SOURCE_MANIFEST_SCHEMA = pa.schema(
    [
        pa.field("snapshot", pa.string(), nullable=False),
        pa.field("source_layout", pa.string(), nullable=False),
        pa.field("pdb_id", pa.string(), nullable=False),
        pa.field("s3_key", pa.string(), nullable=False),
        pa.field("size_bytes", pa.int64(), nullable=False),
        pa.field("etag", pa.string(), nullable=False),
        pa.field(
            "last_modified_utc",
            pa.timestamp("us", tz="UTC"),
            nullable=False,
        ),
        pa.field(
            "manifest_generated_at_utc",
            pa.timestamp("us", tz="UTC"),
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": b"pdbclean_source_manifest",
        b"schema_version": SOURCE_MANIFEST_SCHEMA_VERSION.encode(),
    },
)
