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


def test_gold_accepted_chain_schema_identity() -> None:
    from pdbclean.schemas import (
        GOLD_ACCEPTED_CHAIN_SCHEMA,
        GOLD_ACCEPTED_CHAIN_SCHEMA_VERSION,
    )

    assert GOLD_ACCEPTED_CHAIN_SCHEMA_VERSION == "1.0"
    assert (
        GOLD_ACCEPTED_CHAIN_SCHEMA.metadata[b"schema_name"]
        == b"pdbclean_gold_accepted_chain"
    )


def test_gold_accepted_chain_schema_required_fields() -> None:
    from pdbclean.schemas import GOLD_ACCEPTED_CHAIN_SCHEMA

    required = {
        "snapshot",
        "pdb_id",
        "model_id",
        "label_chain_id",
        "original_start_label_seq_id",
        "original_end_label_seq_id",
        "retained_start_label_seq_id",
        "retained_end_label_seq_id",
        "retained_residue_count",
        "retained_label_seq_ids",
        "retained_sequence",
        "terminal_trimmed",
        "dirty_residue_count",
        "dirty_rule_ids",
        "source_mmcif_key",
        "source_etag",
        "cleaning_protocol",
        "pipeline_git_commit",
    }

    assert required.issubset(set(GOLD_ACCEPTED_CHAIN_SCHEMA.names))

    for name in (
        "snapshot",
        "pdb_id",
        "model_id",
        "label_chain_id",
        "retained_start_label_seq_id",
        "retained_end_label_seq_id",
        "retained_residue_count",
        "retained_label_seq_ids",
        "retained_sequence",
    ):
        assert GOLD_ACCEPTED_CHAIN_SCHEMA.field(name).nullable is False


def test_gold_rejected_chain_schema_identity() -> None:
    from pdbclean.schemas import (
        GOLD_REJECTED_CHAIN_SCHEMA,
        GOLD_REJECTED_CHAIN_SCHEMA_VERSION,
    )

    assert GOLD_REJECTED_CHAIN_SCHEMA_VERSION == "1.0"
    assert (
        GOLD_REJECTED_CHAIN_SCHEMA.metadata[b"schema_name"]
        == b"pdbclean_gold_rejected_chain"
    )

    required = {
        "terminal_status",
        "terminal_reason",
        "terminal_stage",
        "missing_label_seq_ids",
        "dirty_residue_count",
        "dirty_rule_ids",
    }

    assert required.issubset(set(GOLD_REJECTED_CHAIN_SCHEMA.names))


def test_gold_dirty_residue_schema_identity() -> None:
    from pdbclean.schemas import (
        GOLD_DIRTY_RESIDUE_SCHEMA,
        GOLD_DIRTY_RESIDUE_SCHEMA_VERSION,
    )

    assert GOLD_DIRTY_RESIDUE_SCHEMA_VERSION == "1.0"
    assert (
        GOLD_DIRTY_RESIDUE_SCHEMA.metadata[b"schema_name"]
        == b"pdbclean_gold_dirty_residue"
    )

    required = {
        "pdb_id",
        "model_id",
        "label_chain_id",
        "label_seq_id",
        "deposited_residue_name",
        "mapped_residue_code",
        "rule_id",
        "dirty_type",
        "cleaning_stage",
        "details_json",
    }

    assert required.issubset(set(GOLD_DIRTY_RESIDUE_SCHEMA.names))


def test_gold_non_candidate_chain_schema_identity() -> None:
    from pdbclean.schemas import (
        GOLD_NON_CANDIDATE_CHAIN_SCHEMA,
        GOLD_NON_CANDIDATE_CHAIN_SCHEMA_VERSION,
    )

    assert GOLD_NON_CANDIDATE_CHAIN_SCHEMA_VERSION == "1.0"
    assert (
        GOLD_NON_CANDIDATE_CHAIN_SCHEMA.metadata[b"schema_name"]
        == b"pdbclean_gold_non_candidate_chain"
    )

    required = {
        "snapshot",
        "pdb_id",
        "model_id",
        "label_chain_id",
        "terminal_status",
        "terminal_reason",
        "terminal_stage",
        "source_mmcif_key",
        "source_etag",
        "cleaning_protocol",
        "pipeline_git_commit",
    }

    assert required.issubset(
        set(GOLD_NON_CANDIDATE_CHAIN_SCHEMA.names)
    )


def test_gold_processing_error_schema_allows_missing_chain_identity() -> None:
    from pdbclean.schemas import GOLD_PROCESSING_ERROR_SCHEMA

    assert (
        GOLD_PROCESSING_ERROR_SCHEMA.metadata[b"schema_name"]
        == b"pdbclean_gold_processing_error"
    )

    assert GOLD_PROCESSING_ERROR_SCHEMA.field("model_id").nullable is True
    assert (
        GOLD_PROCESSING_ERROR_SCHEMA.field("label_chain_id").nullable
        is True
    )


def test_gold_quality_schema_canonical_identity_types_match() -> None:
    from pdbclean.schemas import (
        GOLD_ACCEPTED_CHAIN_SCHEMA,
        GOLD_DIRTY_RESIDUE_SCHEMA,
        GOLD_NON_CANDIDATE_CHAIN_SCHEMA,
        GOLD_REJECTED_CHAIN_SCHEMA,
    )

    schemas = (
        GOLD_ACCEPTED_CHAIN_SCHEMA,
        GOLD_REJECTED_CHAIN_SCHEMA,
        GOLD_NON_CANDIDATE_CHAIN_SCHEMA,
        GOLD_DIRTY_RESIDUE_SCHEMA,
    )

    for schema in schemas:
        assert schema.field("snapshot").type == pa.string()
        assert schema.field("pdb_id").type == pa.string()
        assert schema.field("model_id").type == pa.int32()
        assert schema.field("label_chain_id").type == pa.string()



def test_stage3_bri_chain_schema_identity() -> None:
    from pdbclean.bri import BRI_COLUMNS
    from pdbclean.schemas import (
        STAGE3_BRI_CHAIN_SCHEMA,
        STAGE3_BRI_CHAIN_SCHEMA_VERSION,
    )

    assert STAGE3_BRI_CHAIN_SCHEMA_VERSION == "1.0"

    metadata = STAGE3_BRI_CHAIN_SCHEMA.metadata

    assert metadata[b"schema_name"] == b"pdbclean_stage3_bri_chain"
    assert metadata[b"schema_version"] == b"1.0"
    assert metadata[b"bri_definition"] == b"MATCH Definition 3.4"
    assert metadata[b"bri_decimal_places"] == b"3"
    assert (
        metadata[b"bri_canonicalization"]
        == b"numpy.around(decimal=3)"
    )
    assert (
        metadata[b"bri_columns"].decode()
        == ",".join(BRI_COLUMNS)
    )


def test_stage3_bri_chain_schema_preserves_stage3_input_lineage() -> None:
    from pdbclean.schemas import (
        STAGE3_BRI_CHAIN_SCHEMA,
        STAGE3_ELIGIBLE_CHAIN_SCHEMA,
    )

    upstream_fields = [
        field
        for field in STAGE3_ELIGIBLE_CHAIN_SCHEMA
        if field.name != "pipeline_git_commit"
    ]

    assert STAGE3_BRI_CHAIN_SCHEMA.names[
        : len(upstream_fields)
    ] == [field.name for field in upstream_fields]

    for expected in upstream_fields:
        observed = STAGE3_BRI_CHAIN_SCHEMA.field(expected.name)

        assert observed.type == expected.type
        assert observed.nullable == expected.nullable

    quality_commit = STAGE3_BRI_CHAIN_SCHEMA.field(
        "quality_pipeline_git_commit"
    )
    upstream_commit = STAGE3_ELIGIBLE_CHAIN_SCHEMA.field(
        "pipeline_git_commit"
    )

    assert quality_commit.type == upstream_commit.type
    assert quality_commit.nullable is False


def test_stage3_bri_chain_schema_has_exact_m_by_9_payload_type() -> None:
    from pdbclean.schemas import STAGE3_BRI_CHAIN_SCHEMA

    expected_row_type = pa.list_(
        pa.field(
            "coordinate",
            pa.float64(),
            nullable=False,
        ),
        9,
    )

    expected_bri_type = pa.list_(
        pa.field(
            "row",
            expected_row_type,
            nullable=False,
        )
    )

    bri_field = STAGE3_BRI_CHAIN_SCHEMA.field("bri")

    assert bri_field.type == expected_bri_type
    assert bri_field.nullable is False

    for name in (
        "quality_pipeline_git_commit",
        "geometric_validation_pipeline_git_commit",
        "geometric_validation_finalizer_git_commit",
        "bri_pipeline_git_commit",
    ):
        field = STAGE3_BRI_CHAIN_SCHEMA.field(name)

        assert field.type == pa.string()
        assert field.nullable is False


def test_stage3_bri_processing_error_schema_identity_and_lineage() -> None:
    from pdbclean.schemas import (
        STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
        STAGE3_BRI_PROCESSING_ERROR_SCHEMA_VERSION,
    )

    schema = STAGE3_BRI_PROCESSING_ERROR_SCHEMA

    assert STAGE3_BRI_PROCESSING_ERROR_SCHEMA_VERSION == "1.0"
    assert (
        schema.metadata[b"schema_name"]
        == b"pdbclean_stage3_bri_processing_error"
    )
    assert schema.metadata[b"schema_version"] == b"1.0"

    required = {
        "snapshot",
        "pdb_id",
        "model_id",
        "label_chain_id",
        "retained_residue_count",
        "retained_label_seq_ids",
        "processing_stage",
        "error_type",
        "error_message",
        "source_mmcif_key",
        "source_etag",
        "cleaning_protocol",
        "quality_pipeline_git_commit",
        "geometric_validation_pipeline_git_commit",
        "geometric_validation_finalizer_git_commit",
        "bri_pipeline_git_commit",
    }

    assert set(schema.names) == required

    for field in schema:
        assert field.nullable is False

    assert schema.field("retained_residue_count").type == pa.int32()
    assert schema.field("retained_label_seq_ids").type == pa.list_(
        pa.field("element", pa.int32())
    )


def test_stage5_brain_chain_schema_preserves_bri_lineage_without_matrix() -> None:
    from pdbclean.schemas import (
        STAGE3_BRI_CHAIN_SCHEMA,
        STAGE5_BRAIN_CHAIN_SCHEMA,
        STAGE5_BRAIN_CHAIN_SCHEMA_VERSION,
    )

    schema = STAGE5_BRAIN_CHAIN_SCHEMA

    assert STAGE5_BRAIN_CHAIN_SCHEMA_VERSION == "1.0"
    assert (
        schema.metadata[b"schema_name"]
        == b"pdbclean_stage5_brain_chain"
    )
    assert schema.metadata[b"schema_version"] == b"1.0"

    upstream_fields = [
        field
        for field in STAGE3_BRI_CHAIN_SCHEMA
        if field.name != "bri"
    ]

    assert schema.names[: len(upstream_fields)] == [
        field.name
        for field in upstream_fields
    ]

    for expected in upstream_fields:
        observed = schema.field(expected.name)

        assert observed.type == expected.type
        assert observed.nullable == expected.nullable

    assert "bri" not in schema.names

    bri_finalizer = schema.field(
        "bri_finalizer_git_commit"
    )
    brain_producer = schema.field(
        "brain_pipeline_git_commit"
    )

    assert bri_finalizer.type == pa.string()
    assert bri_finalizer.nullable is False
    assert brain_producer.type == pa.string()
    assert brain_producer.nullable is False


def test_stage5_brain_payload_is_exact_nonnullable_float64_vector9() -> None:
    from pdbclean.schemas import (
        STAGE5_BRAIN_CHAIN_SCHEMA,
    )

    expected_type = pa.list_(
        pa.field(
            "coordinate",
            pa.float64(),
            nullable=False,
        ),
        9,
    )

    field = STAGE5_BRAIN_CHAIN_SCHEMA.field(
        "brain"
    )

    assert field.type == expected_type
    assert field.nullable is False

    metadata = STAGE5_BRAIN_CHAIN_SCHEMA.metadata

    assert metadata[b"brain_definition"] == (
        b"MATCH Definition 5.1"
    )
    assert metadata[b"brain_rows"] == (
        b"rows 2..m; first BRI row excluded"
    )
    assert metadata[b"brain_result_rounding"] == b"none"
    assert metadata[b"brain_minimum_m"] == b"2"


def test_stage5_brain_undefined_schema_is_explicit_terminal_outcome() -> None:
    from pdbclean.schemas import (
        STAGE3_BRI_CHAIN_SCHEMA,
        STAGE5_BRAIN_UNDEFINED_CHAIN_SCHEMA,
        STAGE5_BRAIN_UNDEFINED_CHAIN_SCHEMA_VERSION,
    )

    schema = STAGE5_BRAIN_UNDEFINED_CHAIN_SCHEMA

    assert STAGE5_BRAIN_UNDEFINED_CHAIN_SCHEMA_VERSION == "1.0"
    assert (
        schema.metadata[b"schema_name"]
        == b"pdbclean_stage5_brain_undefined_chain"
    )
    assert schema.metadata[b"schema_version"] == b"1.0"
    assert schema.metadata[b"brain_defined"] == b"false"
    assert schema.metadata[b"brain_undefined_condition"] == b"m=1"

    upstream_fields = [
        field
        for field in STAGE3_BRI_CHAIN_SCHEMA
        if field.name != "bri"
    ]

    assert schema.names[: len(upstream_fields)] == [
        field.name
        for field in upstream_fields
    ]

    assert "bri" not in schema.names
    assert "brain" not in schema.names

    reason = schema.field(
        "undefined_reason"
    )

    assert reason.type == pa.string()
    assert reason.nullable is False

    bri_finalizer = schema.field(
        "bri_finalizer_git_commit"
    )
    brain_producer = schema.field(
        "brain_pipeline_git_commit"
    )

    assert bri_finalizer.type == pa.string()
    assert bri_finalizer.nullable is False
    assert brain_producer.type == pa.string()
    assert brain_producer.nullable is False


def test_stage5_brain_processing_error_schema_has_traceable_lineage() -> None:
    from pdbclean.schemas import (
        STAGE3_BRI_CHAIN_SCHEMA,
        STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA,
        STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA_VERSION,
    )

    schema = STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA

    assert STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA_VERSION == "1.0"
    assert (
        schema.metadata[b"schema_name"]
        == b"pdbclean_stage5_brain_processing_error"
    )
    assert schema.metadata[b"schema_version"] == b"1.0"

    upstream_names = (
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

    for name in upstream_names:
        observed = schema.field(name)
        expected = STAGE3_BRI_CHAIN_SCHEMA.field(name)

        assert observed.type == expected.type
        assert observed.nullable == expected.nullable

    for name in (
        "bri_finalizer_git_commit",
        "brain_pipeline_git_commit",
        "processing_stage",
        "error_type",
        "error_message",
    ):
        field = schema.field(name)

        assert field.type == pa.string()
        assert field.nullable is False

    for field in schema:
        assert field.nullable is False
