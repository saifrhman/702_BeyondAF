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


def test_silver_chain_schema_identity() -> None:
    from pdbclean.schemas import (
        SILVER_CHAIN_SCHEMA,
        SILVER_CHAIN_SCHEMA_VERSION,
    )

    assert SILVER_CHAIN_SCHEMA_VERSION == "1.0"
    assert (
        SILVER_CHAIN_SCHEMA.metadata[b"schema_name"]
        == b"pdbclean_silver_chain"
    )
    assert (
        SILVER_CHAIN_SCHEMA.metadata[b"schema_version"]
        == b"1.0"
    )


def test_silver_chain_schema_required_fields() -> None:
    from pdbclean.schemas import SILVER_CHAIN_SCHEMA

    required = {
        "snapshot",
        "pdb_id",
        "model_id",
        "label_chain_id",
        "source_mmcif_key",
        "source_etag",
        "observed_residue_count",
        "observed_label_seq_ids",
        "observed_residue_names",
        "partial_occupancy_atom_count",
        "alternate_location_atom_count",
        "missing_backbone_atom_count",
        "duplicate_backbone_atom_count",
        "missing_internal_label_seq_ids",
        "nonstandard_residue_names",
        "atom_count",
    }

    assert required.issubset(set(SILVER_CHAIN_SCHEMA.names))


def test_silver_chain_canonical_identity_is_non_nullable() -> None:
    from pdbclean.schemas import SILVER_CHAIN_SCHEMA

    for name in (
        "snapshot",
        "pdb_id",
        "model_id",
        "label_chain_id",
    ):
        assert SILVER_CHAIN_SCHEMA.field(name).nullable is False
