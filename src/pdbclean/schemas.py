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

# ---------------------------------------------------------------------------
# Gold: Protocol 3.2 quality-cleaning outputs
# ---------------------------------------------------------------------------

GOLD_ACCEPTED_CHAIN_SCHEMA_VERSION = "1.0"

GOLD_ACCEPTED_CHAIN_SCHEMA = pa.schema(
    [
        # Canonical identity / provenance
        pa.field("snapshot", pa.string(), nullable=False),
        pa.field("pdb_id", pa.string(), nullable=False),
        pa.field("model_id", pa.int32(), nullable=False),
        pa.field("label_chain_id", pa.string(), nullable=False),
        pa.field("auth_chain_id", pa.string(), nullable=True),
        pa.field("entity_id", pa.string(), nullable=True),

        # Original observed chain range
        pa.field(
            "original_start_label_seq_id",
            pa.int32(),
            nullable=True,
        ),
        pa.field(
            "original_end_label_seq_id",
            pa.int32(),
            nullable=True,
        ),

        # Clean retained chain
        pa.field(
            "retained_start_label_seq_id",
            pa.int32(),
            nullable=False,
        ),
        pa.field(
            "retained_end_label_seq_id",
            pa.int32(),
            nullable=False,
        ),
        pa.field(
            "retained_residue_count",
            pa.int32(),
            nullable=False,
        ),
        pa.field(
            "retained_label_seq_ids",
            pa.list_(pa.int32()),
            nullable=False,
        ),
        pa.field(
            "retained_sequence",
            pa.string(),
            nullable=False,
        ),

        # Cleaning lineage
        pa.field(
            "terminal_trimmed",
            pa.bool_(),
            nullable=False,
        ),
        pa.field(
            "dirty_residue_count",
            pa.int32(),
            nullable=False,
        ),
        pa.field(
            "dirty_rule_ids",
            pa.list_(pa.string()),
            nullable=False,
        ),

        # Source / release lineage
        pa.field("source_mmcif_key", pa.string(), nullable=False),
        pa.field("source_etag", pa.string(), nullable=False),
        pa.field("cleaning_protocol", pa.string(), nullable=False),
        pa.field("pipeline_git_commit", pa.string(), nullable=False),
    ],
    metadata={
        b"schema_name": b"pdbclean_gold_accepted_chain",
        b"schema_version": GOLD_ACCEPTED_CHAIN_SCHEMA_VERSION.encode(),
    },
)


GOLD_REJECTED_CHAIN_SCHEMA_VERSION = "1.0"

GOLD_REJECTED_CHAIN_SCHEMA = pa.schema(
    [
        # Canonical identity / provenance
        pa.field("snapshot", pa.string(), nullable=False),
        pa.field("pdb_id", pa.string(), nullable=False),
        pa.field("model_id", pa.int32(), nullable=False),
        pa.field("label_chain_id", pa.string(), nullable=False),
        pa.field("auth_chain_id", pa.string(), nullable=True),
        pa.field("entity_id", pa.string(), nullable=True),

        # Terminal scientific outcome
        pa.field("terminal_status", pa.string(), nullable=False),
        pa.field("terminal_reason", pa.string(), nullable=False),
        pa.field("terminal_stage", pa.string(), nullable=False),
        pa.field(
            "missing_label_seq_ids",
            pa.list_(pa.int32()),
            nullable=False,
        ),

        # Cleaning lineage accumulated before rejection
        pa.field(
            "dirty_residue_count",
            pa.int32(),
            nullable=False,
        ),
        pa.field(
            "dirty_rule_ids",
            pa.list_(pa.string()),
            nullable=False,
        ),

        # Source / release lineage
        pa.field("source_mmcif_key", pa.string(), nullable=False),
        pa.field("source_etag", pa.string(), nullable=False),
        pa.field("cleaning_protocol", pa.string(), nullable=False),
        pa.field("pipeline_git_commit", pa.string(), nullable=False),
    ],
    metadata={
        b"schema_name": b"pdbclean_gold_rejected_chain",
        b"schema_version": GOLD_REJECTED_CHAIN_SCHEMA_VERSION.encode(),
    },
)


GOLD_DIRTY_RESIDUE_SCHEMA_VERSION = "1.0"

GOLD_DIRTY_RESIDUE_SCHEMA = pa.schema(
    [
        # Canonical source-chain identity
        pa.field("snapshot", pa.string(), nullable=False),
        pa.field("pdb_id", pa.string(), nullable=False),
        pa.field("model_id", pa.int32(), nullable=False),
        pa.field("label_chain_id", pa.string(), nullable=False),
        pa.field("auth_chain_id", pa.string(), nullable=True),
        pa.field("entity_id", pa.string(), nullable=True),

        # Dirty residue identity
        pa.field("label_seq_id", pa.int32(), nullable=False),
        pa.field(
            "deposited_residue_name",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "mapped_residue_code",
            pa.string(),
            nullable=True,
        ),

        # Protocol decision
        pa.field("rule_id", pa.string(), nullable=False),
        pa.field("dirty_type", pa.string(), nullable=False),
        pa.field("cleaning_stage", pa.string(), nullable=False),

        # Rule-specific evidence serialized as a JSON object.
        pa.field("details_json", pa.string(), nullable=False),

        # Source lineage
        pa.field("source_mmcif_key", pa.string(), nullable=False),
        pa.field("source_etag", pa.string(), nullable=False),
    ],
    metadata={
        b"schema_name": b"pdbclean_gold_dirty_residue",
        b"schema_version": GOLD_DIRTY_RESIDUE_SCHEMA_VERSION.encode(),
    },
)


GOLD_NON_CANDIDATE_CHAIN_SCHEMA_VERSION = "1.0"

GOLD_NON_CANDIDATE_CHAIN_SCHEMA = pa.schema(
    [
        # Canonical identity / provenance
        pa.field("snapshot", pa.string(), nullable=False),
        pa.field("pdb_id", pa.string(), nullable=False),
        pa.field("model_id", pa.int32(), nullable=False),
        pa.field("label_chain_id", pa.string(), nullable=False),
        pa.field("auth_chain_id", pa.string(), nullable=True),
        pa.field("entity_id", pa.string(), nullable=True),

        # Candidate-selection outcome
        pa.field("terminal_status", pa.string(), nullable=False),
        pa.field("terminal_reason", pa.string(), nullable=False),
        pa.field("terminal_stage", pa.string(), nullable=False),

        # Source / release lineage
        pa.field("source_mmcif_key", pa.string(), nullable=False),
        pa.field("source_etag", pa.string(), nullable=False),
        pa.field("cleaning_protocol", pa.string(), nullable=False),
        pa.field("pipeline_git_commit", pa.string(), nullable=False),
    ],
    metadata={
        b"schema_name": b"pdbclean_gold_non_candidate_chain",
        b"schema_version": GOLD_NON_CANDIDATE_CHAIN_SCHEMA_VERSION.encode(),
    },
)


GOLD_PROCESSING_ERROR_SCHEMA_VERSION = "1.0"

GOLD_PROCESSING_ERROR_SCHEMA = pa.schema(
    [
        pa.field("snapshot", pa.string(), nullable=False),
        pa.field("pdb_id", pa.string(), nullable=False),

        # Chain identity may be unavailable for entry/parser failures.
        pa.field("model_id", pa.int32(), nullable=True),
        pa.field("label_chain_id", pa.string(), nullable=True),

        pa.field("processing_stage", pa.string(), nullable=False),
        pa.field("error_type", pa.string(), nullable=False),
        pa.field("error_message", pa.string(), nullable=False),

        pa.field("source_mmcif_key", pa.string(), nullable=False),
        pa.field("source_etag", pa.string(), nullable=False),
        pa.field("pipeline_git_commit", pa.string(), nullable=False),
    ],
    metadata={
        b"schema_name": b"pdbclean_gold_processing_error",
        b"schema_version": GOLD_PROCESSING_ERROR_SCHEMA_VERSION.encode(),
    },
)

