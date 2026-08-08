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


SILVER_CHAIN_SCHEMA_VERSION = "1.0"

SILVER_CHAIN_SCHEMA = pa.schema(
    [
        # Provenance and canonical identity
        pa.field("snapshot", pa.string(), nullable=False),
        pa.field("pdb_id", pa.string(), nullable=False),
        pa.field("model_id", pa.int32(), nullable=False),
        pa.field("label_chain_id", pa.string(), nullable=False),
        pa.field("auth_chain_id", pa.string(), nullable=True),
        pa.field("entity_id", pa.string(), nullable=True),
        pa.field("polymer_type", pa.string(), nullable=True),

        # Source lineage
        pa.field("source_layout", pa.string(), nullable=False),
        pa.field("source_mmcif_key", pa.string(), nullable=False),
        pa.field("source_etag", pa.string(), nullable=False),

        # Observed residue information
        pa.field("observed_start_label_seq_id", pa.int32(), nullable=True),
        pa.field("observed_end_label_seq_id", pa.int32(), nullable=True),
        pa.field("observed_residue_count", pa.int32(), nullable=False),
        pa.field(
            "observed_label_seq_ids",
            pa.list_(pa.int32()),
            nullable=False,
        ),
        pa.field(
            "observed_residue_names",
            pa.list_(pa.string()),
            nullable=False,
        ),

        # Sequence derived from observed residues
        pa.field("observed_sequence", pa.string(), nullable=True),

        # Occupancy / alternate-location observations
        pa.field("minimum_occupancy", pa.float64(), nullable=True),
        pa.field("partial_occupancy_atom_count", pa.int32(), nullable=False),
        pa.field("alternate_location_atom_count", pa.int32(), nullable=False),

        # Backbone observations
        pa.field("missing_backbone_atom_count", pa.int32(), nullable=False),
        pa.field("duplicate_backbone_atom_count", pa.int32(), nullable=False),
        pa.field(
            "minimum_consecutive_backbone_distance",
            pa.float64(),
            nullable=True,
        ),

        # Residue continuity
        pa.field(
            "missing_internal_label_seq_ids",
            pa.list_(pa.int32()),
            nullable=False,
        ),

        # Amino-acid observations
        pa.field(
            "nonstandard_residue_names",
            pa.list_(pa.string()),
            nullable=False,
        ),

        # Parsing/accounting
        pa.field("atom_count", pa.int32(), nullable=False),
    ],
    metadata={
        b"schema_name": b"pdbclean_silver_chain",
        b"schema_version": SILVER_CHAIN_SCHEMA_VERSION.encode(),
    },
)
