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


# Silver is the deterministic in-memory parsed representation produced by
# mmcif_parser. This pipeline version does not persist a separate Silver
# dataset and therefore defines no persisted Arrow schema for Silver.

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


# ---------------------------------------------------------------------------
# Step 2: post-cleaning geometric-validation audit outputs
# ---------------------------------------------------------------------------

GEOMETRIC_VALIDATION_AUDIT_SCHEMA_VERSION = "1.0"

GEOMETRIC_VALIDATION_AUDIT_SCHEMA = pa.schema(
    [
        # Exact accepted-chain identity / retained lineage
        pa.field("snapshot", pa.string(), nullable=False),
        pa.field("pdb_id", pa.string(), nullable=False),
        pa.field("model_id", pa.int32(), nullable=False),
        pa.field("label_chain_id", pa.string(), nullable=False),
        pa.field(
            "retained_residue_count",
            pa.int32(),
            nullable=False,
        ),
        pa.field(
            "retained_label_seq_ids",
            pa.list_(pa.field("element", pa.int32())),
            nullable=False,
        ),

        # Immutable source / cleaning lineage
        pa.field("source_mmcif_key", pa.string(), nullable=False),
        pa.field("source_etag", pa.string(), nullable=False),
        pa.field("cleaning_protocol", pa.string(), nullable=False),
        pa.field(
            "quality_pipeline_git_commit",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "geometric_validation_pipeline_git_commit",
            pa.string(),
            nullable=False,
        ),

        # Scientific configuration used for this audit
        pa.field(
            "configured_minimum_backbone_distance_angstrom",
            pa.float64(),
            nullable=False,
        ),
        pa.field(
            "configured_minimum_triangle_angle_degrees",
            pa.float64(),
            nullable=False,
        ),

        # Scientific audit outcome
        pa.field("passed", pa.bool_(), nullable=False),
        pa.field(
            "minimum_observed_backbone_distance_angstrom",
            pa.float64(),
            nullable=True,
        ),
        pa.field(
            "minimum_observed_triangle_angle_degrees",
            pa.float64(),
            nullable=True,
        ),
        pa.field(
            "minimum_observed_basis_h_norm_angstrom",
            pa.float64(),
            nullable=True,
        ),
        pa.field("violation_count", pa.int32(), nullable=False),
        pa.field(
            "violation_types",
            pa.list_(pa.field("element", pa.string())),
            nullable=False,
        ),
        pa.field(
            "violation_residue_ids",
            pa.list_(pa.field("element", pa.int32())),
            nullable=False,
        ),
        pa.field(
            "violation_details",
            pa.list_(pa.field("element", pa.string())),
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": b"pdbclean_geometric_validation_audit",
        b"schema_version": (
            GEOMETRIC_VALIDATION_AUDIT_SCHEMA_VERSION.encode()
        ),
    },
)


GEOMETRIC_VALIDATION_PROCESSING_ERROR_SCHEMA_VERSION = "1.0"

GEOMETRIC_VALIDATION_PROCESSING_ERROR_SCHEMA = pa.schema(
    [
        # Every Step-2 error corresponds to one accepted Gold chain.
        pa.field("snapshot", pa.string(), nullable=False),
        pa.field("pdb_id", pa.string(), nullable=False),
        pa.field("model_id", pa.int32(), nullable=False),
        pa.field("label_chain_id", pa.string(), nullable=False),

        pa.field("processing_stage", pa.string(), nullable=False),
        pa.field("error_type", pa.string(), nullable=False),
        pa.field("error_message", pa.string(), nullable=False),

        # Source / upstream / Step-2 provenance
        pa.field("source_mmcif_key", pa.string(), nullable=False),
        pa.field("source_etag", pa.string(), nullable=False),
        pa.field("cleaning_protocol", pa.string(), nullable=False),
        pa.field(
            "quality_pipeline_git_commit",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "geometric_validation_pipeline_git_commit",
            pa.string(),
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": (
            b"pdbclean_geometric_validation_processing_error"
        ),
        b"schema_version": (
            GEOMETRIC_VALIDATION_PROCESSING_ERROR_SCHEMA_VERSION.encode()
        ),
    },
)


# ---------------------------------------------------------------------------
# Stage 2 finalization: canonical Stage-3 input population
# ---------------------------------------------------------------------------

STAGE3_ELIGIBLE_CHAIN_SCHEMA_VERSION = "1.0"

# Preserve the complete authoritative Stage-1 accepted-chain record.
STAGE3_ELIGIBLE_CHAIN_SCHEMA = (
    GOLD_ACCEPTED_CHAIN_SCHEMA.with_metadata(
        {
            b"schema_name": b"pdbclean_stage3_eligible_chain",
            b"schema_version": (
                STAGE3_ELIGIBLE_CHAIN_SCHEMA_VERSION.encode()
            ),
        }
    )
)


STAGE3_QUARANTINED_CHAIN_SCHEMA_VERSION = "1.0"

# Quarantined chains preserve the same authoritative Stage-1 lineage.
# The geometric evidence explaining quarantine remains in the Step-2
# audit dataset.
STAGE3_QUARANTINED_CHAIN_SCHEMA = (
    GOLD_ACCEPTED_CHAIN_SCHEMA.with_metadata(
        {
            b"schema_name": b"pdbclean_stage3_quarantined_chain",
            b"schema_version": (
                STAGE3_QUARANTINED_CHAIN_SCHEMA_VERSION.encode()
            ),
        }
    )
)



# ---------------------------------------------------------------------------
# Stage 3: Definition 3.4 BRI production artifacts
# ---------------------------------------------------------------------------

STAGE3_BRI_CHAIN_SCHEMA_VERSION = "1.0"

# One BRI residue row always contains the nine strong coordinates from
# MATCH Definition 3.4. The outer list has one row per retained residue.
_STAGE3_BRI_ROW_TYPE = pa.list_(
    pa.field(
        "coordinate",
        pa.float64(),
        nullable=False,
    ),
    9,
)

_STAGE3_BRI_MATRIX_TYPE = pa.list_(
    pa.field(
        "row",
        _STAGE3_BRI_ROW_TYPE,
        nullable=False,
    )
)

# Preserve all canonical Stage-3 input lineage. Rename the generic
# Stage-1 producer field to make its role explicit once additional
# producer commits are attached.
_STAGE3_BRI_UPSTREAM_FIELDS = [
    field
    for field in STAGE3_ELIGIBLE_CHAIN_SCHEMA
    if field.name != "pipeline_git_commit"
]

STAGE3_BRI_CHAIN_SCHEMA = pa.schema(
    [
        *_STAGE3_BRI_UPSTREAM_FIELDS,

        # Upstream and Stage-3 producer provenance
        pa.field(
            "quality_pipeline_git_commit",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "geometric_validation_pipeline_git_commit",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "geometric_validation_finalizer_git_commit",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "bri_pipeline_git_commit",
            pa.string(),
            nullable=False,
        ),

        # Canonical Definition 3.4 BRI: m rows x 9 coordinates
        pa.field(
            "bri",
            _STAGE3_BRI_MATRIX_TYPE,
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": b"pdbclean_stage3_bri_chain",
        b"schema_version": (
            STAGE3_BRI_CHAIN_SCHEMA_VERSION.encode()
        ),
        b"bri_definition": b"MATCH Definition 3.4",
        b"bri_columns": (
            b"x(N),y(N),z(N),x(A),y(A),z(A),x(C),y(C),z(C)"
        ),
        b"bri_decimal_places": b"3",
        b"bri_canonicalization": b"numpy.around(decimal=3)",
    },
)


STAGE3_BRI_PROCESSING_ERROR_SCHEMA_VERSION = "1.0"

STAGE3_BRI_PROCESSING_ERROR_SCHEMA = pa.schema(
    [
        # Every error corresponds to exactly one canonical Stage-3
        # eligible-chain identity.
        pa.field("snapshot", pa.string(), nullable=False),
        pa.field("pdb_id", pa.string(), nullable=False),
        pa.field("model_id", pa.int32(), nullable=False),
        pa.field("label_chain_id", pa.string(), nullable=False),
        pa.field(
            "retained_residue_count",
            pa.int32(),
            nullable=False,
        ),
        pa.field(
            "retained_label_seq_ids",
            pa.list_(pa.field("element", pa.int32())),
            nullable=False,
        ),

        # Terminal processing failure
        pa.field("processing_stage", pa.string(), nullable=False),
        pa.field("error_type", pa.string(), nullable=False),
        pa.field("error_message", pa.string(), nullable=False),

        # Immutable source and complete producer lineage
        pa.field("source_mmcif_key", pa.string(), nullable=False),
        pa.field("source_etag", pa.string(), nullable=False),
        pa.field("cleaning_protocol", pa.string(), nullable=False),
        pa.field(
            "quality_pipeline_git_commit",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "geometric_validation_pipeline_git_commit",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "geometric_validation_finalizer_git_commit",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "bri_pipeline_git_commit",
            pa.string(),
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": (
            b"pdbclean_stage3_bri_processing_error"
        ),
        b"schema_version": (
            STAGE3_BRI_PROCESSING_ERROR_SCHEMA_VERSION.encode()
        ),
    },
)


# ---------------------------------------------------------------------------
# Stage 5: MATCH Definition 5.1 Brain production artifacts
# ---------------------------------------------------------------------------

STAGE5_BRAIN_CHAIN_SCHEMA_VERSION = "1.0"
STAGE5_BRAIN_UNDEFINED_CHAIN_SCHEMA_VERSION = "1.0"
STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA_VERSION = "1.0"

# Definition 5.1 has exactly nine strong-coordinate averages.
_STAGE5_BRAIN_VECTOR_TYPE = pa.list_(
    pa.field(
        "coordinate",
        pa.float64(),
        nullable=False,
    ),
    9,
)

# Brain is derived entirely from the canonical Stage-3 BRI. Preserve all
# Stage-3 lineage and producer provenance but do not duplicate the full
# m x 9 BRI matrix into the Brain artifact.
_STAGE5_BRAIN_UPSTREAM_FIELDS = [
    field
    for field in STAGE3_BRI_CHAIN_SCHEMA
    if field.name != "bri"
]

STAGE5_BRAIN_CHAIN_SCHEMA = pa.schema(
    [
        *_STAGE5_BRAIN_UPSTREAM_FIELDS,

        pa.field(
            "brain_pipeline_git_commit",
            pa.string(),
            nullable=False,
        ),

        # MATCH Definition 5.1:
        # mean of each of the nine BRI columns over rows 2..m.
        pa.field(
            "brain",
            _STAGE5_BRAIN_VECTOR_TYPE,
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": b"pdbclean_stage5_brain_chain",
        b"schema_version": (
            STAGE5_BRAIN_CHAIN_SCHEMA_VERSION.encode()
        ),
        b"brain_definition": b"MATCH Definition 5.1",
        b"brain_columns": (
            b"x(N),y(N),z(N),x(A),y(A),z(A),x(C),y(C),z(C)"
        ),
        b"brain_input": b"canonical Stage-3 BRI",
        b"brain_rows": b"rows 2..m; first BRI row excluded",
        b"brain_result_rounding": b"none",
        b"brain_minimum_m": b"2",
    },
)


# m=1 is not an error. Definition 5.1 excludes the first BRI row, so
# there are no rows left to average. Preserve the chain explicitly so
# downstream exact comparison can bypass the Brain prefilter.
STAGE5_BRAIN_UNDEFINED_CHAIN_SCHEMA = pa.schema(
    [
        *_STAGE5_BRAIN_UPSTREAM_FIELDS,

        pa.field(
            "brain_pipeline_git_commit",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "undefined_reason",
            pa.string(),
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": b"pdbclean_stage5_brain_undefined_chain",
        b"schema_version": (
            STAGE5_BRAIN_UNDEFINED_CHAIN_SCHEMA_VERSION.encode()
        ),
        b"brain_definition": b"MATCH Definition 5.1",
        b"brain_defined": b"false",
        b"brain_undefined_condition": b"m=1",
    },
)


# Processing errors are distinct from mathematically undefined m=1
# chains. Preserve the exact identity/source/provenance fields required
# to trace each failed Stage-5 input back to its canonical BRI record.
_STAGE5_BRAIN_ERROR_UPSTREAM_FIELD_NAMES = (
    "snapshot",
    "pdb_id",
    "model_id",
    "label_chain_id",
    "retained_residue_count",
    "retained_label_seq_ids",
    "source_mmcif_key",
    "source_etag",
    "cleaning_protocol",
    "quality_pipeline_git_commit",
    "geometric_validation_pipeline_git_commit",
    "geometric_validation_finalizer_git_commit",
    "bri_pipeline_git_commit",
)

_STAGE5_BRAIN_ERROR_UPSTREAM_FIELDS = [
    STAGE3_BRI_CHAIN_SCHEMA.field(name)
    for name in _STAGE5_BRAIN_ERROR_UPSTREAM_FIELD_NAMES
]

STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA = pa.schema(
    [
        *_STAGE5_BRAIN_ERROR_UPSTREAM_FIELDS,

        pa.field(
            "brain_pipeline_git_commit",
            pa.string(),
            nullable=False,
        ),

        pa.field(
            "processing_stage",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "error_type",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "error_message",
            pa.string(),
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": b"pdbclean_stage5_brain_processing_error",
        b"schema_version": (
            STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA_VERSION.encode()
        ),
    },
)
