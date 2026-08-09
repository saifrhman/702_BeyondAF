"""Tests for Protocol 3.2 Gold materialization."""

import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.cleaning import clean_protocol32_chain
from pdbclean.gold import (
    GoldChainRecords,
    GoldProvenance,
    QualityTaskContext,
    build_quality_task_summary,
    gold_records_to_tables,
    materialize_gold_chain,
    write_gold_quality_shards,
    write_quality_task_summary_atomic,
)
from pdbclean.mmcif_parser import AtomObservation, ChainObservation
from pdbclean.schemas import (
    GOLD_ACCEPTED_CHAIN_SCHEMA,
    GOLD_DIRTY_RESIDUE_SCHEMA,
    GOLD_NON_CANDIDATE_CHAIN_SCHEMA,
    GOLD_PROCESSING_ERROR_SCHEMA,
    GOLD_REJECTED_CHAIN_SCHEMA,
)


def _provenance() -> GoldProvenance:
    return GoldProvenance(
        snapshot="20260101",
        source_mmcif_key="test/test.cif.gz",
        source_etag="etag123",
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
    )



def _atom(
    *,
    atom_name: str,
    residue_id: int,
) -> AtomObservation:
    offset = {
        "N": 0.0,
        "CA": 1.0,
        "C": 2.0,
    }[atom_name]

    return AtomObservation(
        model_id=1,
        label_chain_id="A",
        auth_chain_id="X",
        entity_id="1",
        label_seq_id=residue_id,
        auth_seq_id=str(residue_id),
        residue_name="ALA",
        atom_name=atom_name,
        alt_id=None,
        occupancy=1.0,
        x=(residue_id - 1) * 10.0 + offset,
        y=0.0,
        z=0.0,
        group_pdb="ATOM",
        occupancy_raw="1.00",
    )



def test_materialize_non_candidate_chain() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        auth_chain_id="X",
        entity_id="1",
        entry_has_polypeptide=True,
        atoms=[],
    )

    result = clean_protocol32_chain(chain)
    records = materialize_gold_chain(
        chain,
        result,
        _provenance(),
    )

    assert records.accepted_chain is None
    assert records.rejected_chain is None
    assert records.dirty_residues == ()

    record = records.non_candidate_chain
    assert record is not None
    assert record["terminal_status"] == "non_candidate"
    assert (
        record["terminal_reason"]
        == "empty_protocol32_backbone_projection"
    )
    assert record["terminal_stage"] == "candidate_selection"

    table = pa.Table.from_pylist(
        [record],
        schema=GOLD_NON_CANDIDATE_CHAIN_SCHEMA,
    )

    assert table.num_rows == 1


def test_materialize_rejected_chain_with_dirty_lineage() -> None:
    atoms = [
        _atom(atom_name=name, residue_id=residue_id)
        for residue_id in (1, 2, 3)
        for name in ("N", "CA", "C")
    ]

    index = next(
        i
        for i, atom in enumerate(atoms)
        if atom.label_seq_id == 2 and atom.atom_name == "CA"
    )
    atoms[index] = AtomObservation(
        **{
            **atoms[index].__dict__,
            "occupancy": 0.5,
            "occupancy_raw": "0.50",
        }
    )

    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        auth_chain_id="X",
        entity_id="1",
        entry_has_polypeptide=True,
        atoms=atoms,
    )

    result = clean_protocol32_chain(chain)
    records = materialize_gold_chain(
        chain,
        result,
        _provenance(),
    )

    assert records.accepted_chain is None
    assert records.non_candidate_chain is None

    rejected = records.rejected_chain
    assert rejected is not None
    assert rejected["terminal_status"] == "rejected"
    assert rejected["terminal_stage"] == "Q003_after_Q002"
    assert rejected["missing_label_seq_ids"] == [2]
    assert rejected["dirty_residue_count"] == 1
    assert rejected["dirty_rule_ids"] == ["Q002"]

    assert len(records.dirty_residues) == 1
    dirty = records.dirty_residues[0]
    assert dirty["label_seq_id"] == 2
    assert dirty["rule_id"] == "Q002"
    assert dirty["dirty_type"] == "disordered"
    assert '"occupancy_raw":"0.50"' in dirty["details_json"]

    rejected_table = pa.Table.from_pylist(
        [rejected],
        schema=GOLD_REJECTED_CHAIN_SCHEMA,
    )
    dirty_table = pa.Table.from_pylist(
        list(records.dirty_residues),
        schema=GOLD_DIRTY_RESIDUE_SCHEMA,
    )

    assert rejected_table.num_rows == 1
    assert dirty_table.num_rows == 1


def test_materialize_accepted_chain_after_terminal_trim() -> None:
    atoms = [
        _atom(atom_name=name, residue_id=residue_id)
        for residue_id in (1, 2, 3)
        for name in ("N", "CA", "C")
    ]

    index = next(
        i
        for i, atom in enumerate(atoms)
        if atom.label_seq_id == 3 and atom.atom_name == "CA"
    )
    atoms[index] = AtomObservation(
        **{
            **atoms[index].__dict__,
            "occupancy": 0.5,
            "occupancy_raw": "0.50",
        }
    )

    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        auth_chain_id="X",
        entity_id="1",
        entry_has_polypeptide=True,
        atoms=atoms,
    )

    result = clean_protocol32_chain(chain)
    records = materialize_gold_chain(
        chain,
        result,
        _provenance(),
    )

    assert result.status == "accepted"
    assert records.rejected_chain is None
    assert records.non_candidate_chain is None

    accepted = records.accepted_chain
    assert accepted is not None

    assert accepted["original_start_label_seq_id"] == 1
    assert accepted["original_end_label_seq_id"] == 3

    assert accepted["retained_start_label_seq_id"] == 1
    assert accepted["retained_end_label_seq_id"] == 2
    assert accepted["retained_residue_count"] == 2
    assert accepted["retained_label_seq_ids"] == [1, 2]
    assert accepted["retained_sequence"] == "AA"

    assert accepted["terminal_trimmed"] is True
    assert accepted["dirty_residue_count"] == 1
    assert accepted["dirty_rule_ids"] == ["Q002"]

    assert len(records.dirty_residues) == 1
    assert records.dirty_residues[0]["label_seq_id"] == 3
    assert records.dirty_residues[0]["rule_id"] == "Q002"

    accepted_table = pa.Table.from_pylist(
        [accepted],
        schema=GOLD_ACCEPTED_CHAIN_SCHEMA,
    )
    dirty_table = pa.Table.from_pylist(
        list(records.dirty_residues),
        schema=GOLD_DIRTY_RESIDUE_SCHEMA,
    )

    assert accepted_table.num_rows == 1
    assert dirty_table.num_rows == 1


def test_materialize_clean_accepted_chain_uses_pinned_bri_mapping() -> None:
    atoms = [
        _atom(atom_name=name, residue_id=residue_id)
        for residue_id in (1, 2)
        for name in ("N", "CA", "C")
    ]

    # Pinned BRI v1.2.2 CCD mapping accepts LLP as lysine ("K").
    atoms = [
        AtomObservation(
            **{
                **atom.__dict__,
                "residue_name": (
                    "LLP"
                    if atom.label_seq_id == 2
                    else atom.residue_name
                ),
            }
        )
        for atom in atoms
    ]

    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        auth_chain_id="X",
        entity_id="1",
        entry_has_polypeptide=True,
        atoms=atoms,
    )

    result = clean_protocol32_chain(chain)
    records = materialize_gold_chain(
        chain,
        result,
        _provenance(),
    )

    assert result.status == "accepted"

    accepted = records.accepted_chain
    assert accepted is not None
    assert accepted["retained_label_seq_ids"] == [1, 2]
    assert accepted["retained_sequence"] == "AK"
    assert accepted["terminal_trimmed"] is False
    assert accepted["dirty_residue_count"] == 0
    assert accepted["dirty_rule_ids"] == []
    assert records.dirty_residues == ()

    table = pa.Table.from_pylist(
        [accepted],
        schema=GOLD_ACCEPTED_CHAIN_SCHEMA,
    )

    assert table.num_rows == 1


def test_gold_records_to_tables_preserves_explicit_schemas() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        auth_chain_id="X",
        entity_id="1",
        entry_has_polypeptide=True,
        atoms=[],
    )

    result = clean_protocol32_chain(chain)
    record = materialize_gold_chain(
        chain,
        result,
        _provenance(),
    )

    tables = gold_records_to_tables([record])

    assert tables.accepted_chains.num_rows == 0
    assert tables.rejected_chains.num_rows == 0
    assert tables.non_candidate_chains.num_rows == 1
    assert tables.dirty_residues.num_rows == 0
    assert tables.processing_errors.num_rows == 0

    assert tables.accepted_chains.schema == GOLD_ACCEPTED_CHAIN_SCHEMA
    assert tables.rejected_chains.schema == GOLD_REJECTED_CHAIN_SCHEMA
    assert (
        tables.non_candidate_chains.schema
        == GOLD_NON_CANDIDATE_CHAIN_SCHEMA
    )
    assert tables.dirty_residues.schema == GOLD_DIRTY_RESIDUE_SCHEMA
    assert tables.processing_errors.schema == GOLD_PROCESSING_ERROR_SCHEMA


def test_gold_records_to_tables_requires_one_terminal_outcome() -> None:
    invalid = GoldChainRecords()

    try:
        gold_records_to_tables([invalid])
    except ValueError as exc:
        assert "exactly one terminal chain-level outcome" in str(exc)
    else:
        raise AssertionError(
            "Expected invalid GoldChainRecords to raise ValueError"
        )


def test_write_gold_quality_shards_preserves_schema_and_rows(
    tmp_path,
) -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        auth_chain_id="X",
        entity_id="1",
        entry_has_polypeptide=True,
        atoms=[],
    )

    result = clean_protocol32_chain(chain)
    record = materialize_gold_chain(
        chain,
        result,
        _provenance(),
    )
    tables = gold_records_to_tables([record])

    outputs = write_gold_quality_shards(
        tables,
        tmp_path,
        task_id=7,
    )

    assert set(outputs) == {
        "accepted",
        "rejected",
        "non_candidates",
        "dirty_residues",
        "errors",
    }

    accepted = pq.read_table(outputs["accepted"])
    rejected = pq.read_table(outputs["rejected"])
    non_candidates = pq.read_table(outputs["non_candidates"])
    dirty = pq.read_table(outputs["dirty_residues"])
    errors = pq.read_table(outputs["errors"])

    assert accepted.schema == GOLD_ACCEPTED_CHAIN_SCHEMA
    assert rejected.schema == GOLD_REJECTED_CHAIN_SCHEMA
    assert non_candidates.schema == GOLD_NON_CANDIDATE_CHAIN_SCHEMA
    assert dirty.schema == GOLD_DIRTY_RESIDUE_SCHEMA
    assert errors.schema == GOLD_PROCESSING_ERROR_SCHEMA

    assert accepted.num_rows == 0
    assert rejected.num_rows == 0
    assert non_candidates.num_rows == 1
    assert dirty.num_rows == 0
    assert errors.num_rows == 0

    assert not list(tmp_path.rglob("*.tmp"))


def test_write_gold_quality_shards_supports_all_empty_tables(
    tmp_path,
) -> None:
    tables = gold_records_to_tables([])

    outputs = write_gold_quality_shards(
        tables,
        tmp_path,
        task_id="empty",
    )

    expected = {
        "accepted": GOLD_ACCEPTED_CHAIN_SCHEMA,
        "rejected": GOLD_REJECTED_CHAIN_SCHEMA,
        "non_candidates": GOLD_NON_CANDIDATE_CHAIN_SCHEMA,
        "dirty_residues": GOLD_DIRTY_RESIDUE_SCHEMA,
        "errors": GOLD_PROCESSING_ERROR_SCHEMA,
    }

    for name, schema in expected.items():
        table = pq.read_table(outputs[name])
        assert table.num_rows == 0
        assert table.schema == schema

    assert not list(tmp_path.rglob("*.tmp"))


def test_gold_records_to_tables_materializes_processing_errors() -> None:
    error = {
        "snapshot": "20260101",
        "pdb_id": "1abc",
        "model_id": None,
        "label_chain_id": None,
        "processing_stage": "mmcif_parse",
        "error_type": "MMCIFParseError",
        "error_message": "synthetic parser failure",
        "source_mmcif_key": "ab/1abc.cif.gz",
        "source_etag": "etag-error",
        "pipeline_git_commit": "deadbeef",
    }

    tables = gold_records_to_tables(
        [],
        processing_errors=[error],
    )

    assert tables.accepted_chains.num_rows == 0
    assert tables.rejected_chains.num_rows == 0
    assert tables.non_candidate_chains.num_rows == 0
    assert tables.dirty_residues.num_rows == 0
    assert tables.processing_errors.num_rows == 1
    assert tables.processing_errors.schema == GOLD_PROCESSING_ERROR_SCHEMA

    row = tables.processing_errors.to_pylist()[0]
    assert row["pdb_id"] == "1abc"
    assert row["model_id"] is None
    assert row["label_chain_id"] is None
    assert row["processing_stage"] == "mmcif_parse"
    assert row["error_type"] == "MMCIFParseError"


def test_build_quality_task_summary_separates_error_levels() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        auth_chain_id="X",
        entity_id="1",
        entry_has_polypeptide=True,
        atoms=[],
    )

    non_candidate = materialize_gold_chain(
        chain,
        clean_protocol32_chain(chain),
        _provenance(),
    )

    chain_error = {
        "snapshot": "20260101",
        "pdb_id": "1abc",
        "model_id": 1,
        "label_chain_id": "A",
        "processing_stage": "quality_cleaning",
        "error_type": "RuntimeError",
        "error_message": "synthetic chain failure",
        "source_mmcif_key": "ab/1abc.cif.gz",
        "source_etag": "etag-chain",
        "pipeline_git_commit": "deadbeef",
    }

    source_error = {
        "snapshot": "20260101",
        "pdb_id": "2def",
        "model_id": None,
        "label_chain_id": None,
        "processing_stage": "mmcif_parse",
        "error_type": "MMCIFParseError",
        "error_message": "synthetic source failure",
        "source_mmcif_key": "de/2def.cif.gz",
        "source_etag": "etag-source",
        "pipeline_git_commit": "deadbeef",
    }

    tables = gold_records_to_tables(
        [non_candidate],
        processing_errors=[chain_error, source_error],
    )

    summary = build_quality_task_summary(
        tables,
        QualityTaskContext(
            task_id="7",
            snapshot="20260101",
            cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
            pipeline_git_commit="deadbeef",
            started_at_utc="2026-08-09T20:00:00Z",
            completed_at_utc="2026-08-09T20:01:00Z",
            runtime_seconds=60.0,
            input_source_object_count=2,
            successful_source_object_count=1,
            failed_source_object_count=1,
            parsed_silver_chain_count=2,
            candidate_entry_count=1,
            candidate_chain_count=1,
            slurm_job_id="12345",
            slurm_array_task_id="7",
            peak_memory_bytes=1024,
        ),
    )

    assert summary["non_candidate_chain_count"] == 1
    assert summary["processing_error_count"] == 2
    assert summary["chain_level_processing_error_count"] == 1
    assert summary["source_entry_processing_error_count"] == 1

    assert summary["source_object_accounting_valid"] is True
    assert summary["parsed_chain_accounting_valid"] is True

    assert summary["processing_errors_by_stage"] == {
        "mmcif_parse": 1,
        "quality_cleaning": 1,
    }
    assert summary["processing_errors_by_type"] == {
        "MMCIFParseError": 1,
        "RuntimeError": 1,
    }

    assert summary["total_gold_record_count"] == 3


def test_build_quality_task_summary_detects_bad_accounting() -> None:
    tables = gold_records_to_tables([])

    summary = build_quality_task_summary(
        tables,
        QualityTaskContext(
            task_id="bad",
            snapshot="20260101",
            cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
            pipeline_git_commit="deadbeef",
            started_at_utc="2026-08-09T20:00:00Z",
            completed_at_utc="2026-08-09T20:00:01Z",
            runtime_seconds=1.0,
            input_source_object_count=2,
            successful_source_object_count=1,
            failed_source_object_count=0,
            parsed_silver_chain_count=1,
            candidate_entry_count=0,
            candidate_chain_count=0,
        ),
    )

    assert summary["source_object_accounting_valid"] is False
    assert summary["parsed_chain_accounting_valid"] is False


def test_write_quality_task_summary_atomic(tmp_path) -> None:
    tables = gold_records_to_tables([])

    summary = build_quality_task_summary(
        tables,
        QualityTaskContext(
            task_id="12",
            snapshot="20260101",
            cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
            pipeline_git_commit="deadbeef",
            started_at_utc="2026-08-09T20:00:00Z",
            completed_at_utc="2026-08-09T20:00:01Z",
            runtime_seconds=1.0,
            input_source_object_count=0,
            successful_source_object_count=0,
            failed_source_object_count=0,
            parsed_silver_chain_count=0,
            candidate_entry_count=0,
            candidate_chain_count=0,
        ),
    )

    output = write_quality_task_summary_atomic(
        summary,
        tmp_path,
    )

    assert output == tmp_path / "summaries" / "task_12.json"
    assert output.exists()
    assert not list(tmp_path.rglob("*.tmp"))

    payload = output.read_text(encoding="utf-8")

    assert payload.endswith("\n")
    assert '"parsed_chain_accounting_valid": true' in payload
    assert '"source_object_accounting_valid": true' in payload

    # Rewriting the same logical summary is deterministic.
    first = payload
    write_quality_task_summary_atomic(summary, tmp_path)
    second = output.read_text(encoding="utf-8")

    assert second == first


def test_write_quality_task_summary_refuses_invalid_accounting(
    tmp_path,
) -> None:
    tables = gold_records_to_tables([])

    summary = build_quality_task_summary(
        tables,
        QualityTaskContext(
            task_id="bad",
            snapshot="20260101",
            cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
            pipeline_git_commit="deadbeef",
            started_at_utc="2026-08-09T20:00:00Z",
            completed_at_utc="2026-08-09T20:00:01Z",
            runtime_seconds=1.0,
            input_source_object_count=1,
            successful_source_object_count=0,
            failed_source_object_count=0,
            parsed_silver_chain_count=0,
            candidate_entry_count=0,
            candidate_chain_count=0,
        ),
    )

    try:
        write_quality_task_summary_atomic(summary, tmp_path)
    except ValueError as exc:
        assert "source-object accounting failed" in str(exc)
    else:
        raise AssertionError(
            "Expected invalid accounting to prevent summary publication"
        )

    assert not (tmp_path / "summaries" / "task_bad.json").exists()
    assert not list(tmp_path.rglob("*.tmp"))


def test_write_quality_task_summary_rejects_unsafe_task_id(
    tmp_path,
) -> None:
    tables = gold_records_to_tables([])

    summary = build_quality_task_summary(
        tables,
        QualityTaskContext(
            task_id="../escape",
            snapshot="20260101",
            cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
            pipeline_git_commit="deadbeef",
            started_at_utc="2026-08-09T20:00:00Z",
            completed_at_utc="2026-08-09T20:00:01Z",
            runtime_seconds=1.0,
            input_source_object_count=0,
            successful_source_object_count=0,
            failed_source_object_count=0,
            parsed_silver_chain_count=0,
            candidate_entry_count=0,
            candidate_chain_count=0,
        ),
    )

    try:
        write_quality_task_summary_atomic(summary, tmp_path)
    except ValueError as exc:
        assert "Unsafe quality-task task_id" in str(exc)
    else:
        raise AssertionError("Expected unsafe task_id to be rejected")
