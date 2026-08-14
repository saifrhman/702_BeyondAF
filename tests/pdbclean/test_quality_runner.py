"""Tests for quality-stage batch orchestration."""

import pytest

from pdbclean.mmcif_parser import AtomObservation, ChainObservation
from pdbclean.quality_runner import (
    QualityRunnerError,
    candidate_accounting,
    select_configured_model_chains,
)


def _chain(
    model_id: int,
    label_chain_id: str,
) -> ChainObservation:
    return ChainObservation(
        pdb_id="test",
        model_id=model_id,
        label_chain_id=label_chain_id,
        entry_has_polypeptide=True,
    )



def _candidate_chain(
    pdb_id: str,
    label_chain_id: str,
) -> ChainObservation:
    atom = AtomObservation(
        model_id=1,
        label_chain_id=label_chain_id,
        auth_chain_id=label_chain_id,
        entity_id="1",
        label_seq_id=1,
        auth_seq_id="1",
        residue_name="ALA",
        atom_name="CA",
        alt_id=None,
        occupancy=1.0,
        x=0.0,
        y=0.0,
        z=0.0,
        group_pdb="ATOM",
        occupancy_raw="1.00",
    )

    return ChainObservation(
        pdb_id=pdb_id,
        model_id=1,
        label_chain_id=label_chain_id,
        entry_has_polypeptide=True,
        atoms=[atom],
    )


def _non_candidate_chain(
    pdb_id: str,
    label_chain_id: str,
) -> ChainObservation:
    atom = AtomObservation(
        model_id=1,
        label_chain_id=label_chain_id,
        auth_chain_id=label_chain_id,
        entity_id="2",
        label_seq_id=None,
        auth_seq_id=None,
        residue_name="SO4",
        atom_name="S",
        alt_id=None,
        occupancy=1.0,
        x=0.0,
        y=0.0,
        z=0.0,
        group_pdb="HETATM",
        occupancy_raw="1.00",
    )

    return ChainObservation(
        pdb_id=pdb_id,
        model_id=1,
        label_chain_id=label_chain_id,
        entry_has_polypeptide=True,
        atoms=[atom],
    )

def test_select_configured_model_chains_keeps_only_model_one() -> None:
    chains = [
        _chain(1, "A"),
        _chain(2, "A"),
        _chain(1, "B"),
        _chain(3, "C"),
    ]

    selected = select_configured_model_chains(
        chains,
        {
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
    )

    assert [
        (chain.model_id, chain.label_chain_id)
        for chain in selected
    ] == [
        (1, "A"),
        (1, "B"),
    ]


def test_select_configured_model_chains_respects_configured_model_id() -> None:
    chains = [
        _chain(1, "A"),
        _chain(2, "A"),
    ]

    selected = select_configured_model_chains(
        chains,
        {
            "models": {
                "policy": "first_model",
                "model_id": 2,
            }
        },
    )

    assert len(selected) == 1
    assert selected[0].model_id == 2


def test_select_configured_model_chains_all_models_keeps_every_model() -> None:
    chains = [
        _chain(1, "A"),
        _chain(2, "A"),
        _chain(1, "B"),
        _chain(3, "C"),
    ]

    selected = select_configured_model_chains(
        chains,
        {
            "models": {
                "policy": "all_models",
            }
        },
    )

    assert [
        (chain.model_id, chain.label_chain_id)
        for chain in selected
    ] == [
        (1, "A"),
        (2, "A"),
        (1, "B"),
        (3, "C"),
    ]


def test_select_configured_model_chains_rejects_unsupported_policy() -> None:
    with pytest.raises(
        QualityRunnerError,
        match="Unsupported model-selection policy",
    ):
        select_configured_model_chains(
            [_chain(1, "A")],
            {
                "models": {
                    "policy": "unsupported",
                }
            },
        )



def test_candidate_accounting_counts_unique_entries_and_chains() -> None:
    chains = [
        _candidate_chain("1abc", "A"),
        _candidate_chain("1abc", "B"),
        _candidate_chain("2def", "A"),
    ]

    entry_count, chain_count = candidate_accounting(chains)

    assert entry_count == 2
    assert chain_count == 3


def test_candidate_accounting_excludes_non_candidates() -> None:
    chains = [
        _candidate_chain("1abc", "A"),
        _non_candidate_chain("1abc", "B"),
        _non_candidate_chain("2def", "A"),
    ]

    entry_count, chain_count = candidate_accounting(chains)

    assert entry_count == 1
    assert chain_count == 1


def test_candidate_accounting_empty_input() -> None:
    assert candidate_accounting([]) == (0, 0)


def _gzip_cif(text: str) -> bytes:
    import gzip

    return gzip.compress(text.encode("utf-8"))


def _gold_provenance():
    from pdbclean.gold import GoldProvenance

    return GoldProvenance(
        snapshot="20260101",
        source_mmcif_key="20260101/path/test.cif.gz",
        source_etag="etag123",
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
    )


def _multimodel_cif_bytes() -> bytes:
    cif = """data_test
#
loop_
_entity_poly.entity_id
_entity_poly.type
1 'polypeptide(L)'
#
loop_
_atom_site.group_PDB
_atom_site.pdbx_PDB_model_num
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.label_comp_id
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.occupancy
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 A X 1 1 1 ALA N  . 1.00 0.0 0.0 0.0
ATOM 1 A X 1 1 1 ALA CA . 1.00 1.0 0.0 0.0
ATOM 1 A X 1 1 1 ALA C  . 1.00 2.0 0.0 0.0
ATOM 2 A X 1 1 1 ALA N  . 1.00 0.0 1.0 0.0
ATOM 2 A X 1 1 1 ALA CA . 1.00 1.0 1.0 0.0
ATOM 2 A X 1 1 1 ALA C  . 1.00 2.0 1.0 0.0
#
"""
    return _gzip_cif(cif)


def test_process_verified_mmcif_selects_only_configured_model() -> None:
    from pdbclean.quality_runner import process_verified_mmcif_bytes

    result = process_verified_mmcif_bytes(
        _multimodel_cif_bytes(),
        pdb_id="TEST",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        provenance=_gold_provenance(),
    )

    assert result.pdb_id == "test"

    # Parser preserves both deposited models.
    assert result.parsed_silver_chain_count == 2

    # Only configured model 1 proceeds to Protocol 3.2.
    assert result.selected_silver_chain_count == 1
    assert result.candidate_entry_count == 1
    assert result.candidate_chain_count == 1

    assert result.source_failed is False
    assert result.processing_errors == ()
    assert len(result.gold_records) == 1

    gold = result.gold_records[0]
    assert gold.accepted_chain is not None
    assert gold.accepted_chain["model_id"] == 1
    assert gold.accepted_chain["retained_sequence"] == "A"


def test_process_verified_mmcif_all_models_processes_every_model() -> None:
    from pdbclean.quality_runner import process_verified_mmcif_bytes

    result = process_verified_mmcif_bytes(
        _multimodel_cif_bytes(),
        pdb_id="TEST",
        selection_config={
            "models": {
                "policy": "all_models",
            }
        },
        provenance=_gold_provenance(),
    )

    # Parser and selection both preserve the two deposited models.
    assert result.parsed_silver_chain_count == 2
    assert result.selected_silver_chain_count == 2

    assert result.candidate_entry_count == 1
    assert result.candidate_chain_count == 2

    assert result.source_failed is False
    assert result.processing_errors == ()
    assert len(result.gold_records) == 2

    assert [
        record.accepted_chain["model_id"]
        for record in result.gold_records
        if record.accepted_chain is not None
    ] == [1, 2]


def test_process_verified_mmcif_records_parse_failure() -> None:
    from pdbclean.quality_runner import process_verified_mmcif_bytes

    result = process_verified_mmcif_bytes(
        b"not-a-gzip-stream",
        pdb_id="BAD1",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        provenance=_gold_provenance(),
    )

    assert result.pdb_id == "bad1"
    assert result.source_failed is True
    assert result.parsed_silver_chain_count == 0
    assert result.selected_silver_chain_count == 0
    assert result.candidate_entry_count == 0
    assert result.candidate_chain_count == 0
    assert result.gold_records == ()

    assert len(result.processing_errors) == 1
    error = result.processing_errors[0]

    assert error["model_id"] is None
    assert error["label_chain_id"] is None
    assert error["processing_stage"] == "mmcif_parse"
    assert error["error_type"] == "MMCIFParseError"


def test_process_verified_mmcif_records_chain_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from pdbclean.quality_runner import process_verified_mmcif_bytes

    def fail_cleaning(chain, **kwargs):
        raise RuntimeError("synthetic cleaning failure")

    monkeypatch.setattr(
        "pdbclean.quality_runner.clean_protocol32_chain",
        fail_cleaning,
    )

    result = process_verified_mmcif_bytes(
        _multimodel_cif_bytes(),
        pdb_id="TEST",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        provenance=_gold_provenance(),
    )

    assert result.source_failed is False
    assert result.parsed_silver_chain_count == 2
    assert result.selected_silver_chain_count == 1
    assert result.candidate_entry_count == 1
    assert result.candidate_chain_count == 1

    assert result.gold_records == ()
    assert len(result.processing_errors) == 1

    error = result.processing_errors[0]

    assert error["pdb_id"] == "test"
    assert error["model_id"] == 1
    assert error["label_chain_id"] == "A"
    assert error["processing_stage"] == "quality_cleaning"
    assert error["error_type"] == "RuntimeError"
    assert error["error_message"] == "synthetic cleaning failure"



def test_process_verified_mmcif_does_not_swallow_unexpected_parser_bug(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from pdbclean.quality_runner import process_verified_mmcif_bytes

    def fail_parser(*args, **kwargs):
        raise RuntimeError("synthetic parser bug")

    monkeypatch.setattr(
        "pdbclean.quality_runner.parse_coordinate_mmcif_bytes",
        fail_parser,
    )

    with pytest.raises(RuntimeError, match="synthetic parser bug"):
        process_verified_mmcif_bytes(
            _multimodel_cif_bytes(),
            pdb_id="TEST",
            selection_config={
                "models": {
                    "policy": "first_model",
                    "model_id": 1,
                }
            },
            provenance=_gold_provenance(),
        )


def _manifest_row() -> dict:
    compressed = _multimodel_cif_bytes()

    return {
        "snapshot": "20260101",
        "source_layout": "divided_mmcif",
        "pdb_id": "TEST",
        "s3_key": "20260101/pub/pdb/data/structures/divided/mmCIF/te/test.cif.gz",
        "size_bytes": len(compressed),
        "etag": "etag123",
    }


def test_process_manifest_source_downloads_verifies_and_processes() -> None:
    from pdbclean.quality_runner import process_manifest_source

    compressed = _multimodel_cif_bytes()
    calls = []

    def downloader(**kwargs):
        calls.append(kwargs)
        return compressed

    result = process_manifest_source(
        _manifest_row(),
        bucket_url="https://example.invalid",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
        timeout_seconds=37,
        downloader=downloader,
    )

    assert calls == [
        {
            "bucket_url": "https://example.invalid",
            "s3_key": (
                "20260101/pub/pdb/data/structures/"
                "divided/mmCIF/te/test.cif.gz"
            ),
            "expected_size_bytes": len(compressed),
            "expected_etag": "etag123",
            "timeout_seconds": 37,
        }
    ]

    assert result.source_failed is False
    assert result.parsed_silver_chain_count == 2
    assert result.selected_silver_chain_count == 1
    assert result.candidate_entry_count == 1
    assert result.candidate_chain_count == 1
    assert len(result.gold_records) == 1
    assert result.processing_errors == ()

    accepted = result.gold_records[0].accepted_chain
    assert accepted is not None
    assert accepted["snapshot"] == "20260101"
    assert accepted["pdb_id"] == "test"
    assert accepted["source_etag"] == "etag123"
    assert accepted["pipeline_git_commit"] == "deadbeef"


def test_process_manifest_source_records_download_verification_failure() -> None:
    from pdbclean.quality_runner import process_manifest_source
    from pdbclean.snapshot import SnapshotError

    def downloader(**kwargs):
        raise SnapshotError("synthetic ETag mismatch")

    result = process_manifest_source(
        _manifest_row(),
        bucket_url="https://example.invalid",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
        downloader=downloader,
    )

    assert result.source_failed is True
    assert result.parsed_silver_chain_count == 0
    assert result.selected_silver_chain_count == 0
    assert result.gold_records == ()

    assert len(result.processing_errors) == 1
    error = result.processing_errors[0]

    assert error["pdb_id"] == "test"
    assert error["model_id"] is None
    assert error["label_chain_id"] is None
    assert error["processing_stage"] == "source_download_verify"
    assert error["error_type"] == "SnapshotError"
    assert error["error_message"] == "synthetic ETag mismatch"


def test_process_manifest_source_does_not_swallow_unexpected_downloader_bug() -> None:
    from pdbclean.quality_runner import process_manifest_source

    def downloader(**kwargs):
        raise RuntimeError("synthetic downloader bug")

    with pytest.raises(RuntimeError, match="synthetic downloader bug"):
        process_manifest_source(
            _manifest_row(),
            bucket_url="https://example.invalid",
            selection_config={
                "models": {
                    "policy": "first_model",
                    "model_id": 1,
                }
            },
            cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
            pipeline_git_commit="deadbeef",
            downloader=downloader,
        )


def test_process_manifest_source_rejects_malformed_manifest_row() -> None:
    from pdbclean.quality_runner import (
        QualityRunnerError,
        process_manifest_source,
    )

    row = _manifest_row()
    del row["etag"]

    with pytest.raises(
        QualityRunnerError,
        match=r"missing required field\(s\): etag",
    ):
        process_manifest_source(
            row,
            bucket_url="https://example.invalid",
            selection_config={
                "models": {
                    "policy": "first_model",
                    "model_id": 1,
                }
            },
            cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
            pipeline_git_commit="deadbeef",
        )


def test_process_manifest_batch_aggregates_sources_and_gold_tables() -> None:
    from pdbclean.gold import GoldProvenance
    from pdbclean.quality_runner import (
        SourceQualityResult,
        process_manifest_batch,
        process_verified_mmcif_bytes,
    )

    rows = [
        {"pdb_id": "GOOD"},
        {"pdb_id": "BAD"},
    ]

    def source_processor(manifest_row, **kwargs):
        if manifest_row["pdb_id"] == "GOOD":
            return process_verified_mmcif_bytes(
                _multimodel_cif_bytes(),
                pdb_id="GOOD",
                selection_config=kwargs["selection_config"],
                provenance=GoldProvenance(
                    snapshot="20260101",
                    source_mmcif_key="good.cif.gz",
                    source_etag="etag-good",
                    cleaning_protocol=kwargs["cleaning_protocol"],
                    pipeline_git_commit=kwargs["pipeline_git_commit"],
                ),
            )

        return SourceQualityResult(
            pdb_id="bad",
            parsed_silver_chain_count=0,
            selected_silver_chain_count=0,
            candidate_entry_count=0,
            candidate_chain_count=0,
            processing_errors=(
                {
                    "snapshot": "20260101",
                    "pdb_id": "bad",
                    "model_id": None,
                    "label_chain_id": None,
                    "processing_stage": "source_download_verify",
                    "error_type": "SnapshotError",
                    "error_message": "synthetic source failure",
                    "source_mmcif_key": "bad.cif.gz",
                    "source_etag": "etag-bad",
                    "pipeline_git_commit": "deadbeef",
                },
            ),
            source_failed=True,
        )

    result = process_manifest_batch(
        rows,
        bucket_url="https://example.invalid",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
        source_processor=source_processor,
    )

    assert result.input_source_object_count == 2
    assert result.successful_source_object_count == 1
    assert result.failed_source_object_count == 1

    assert result.parsed_silver_chain_count == 2
    assert result.selected_silver_chain_count == 1
    assert result.candidate_entry_count == 1
    assert result.candidate_chain_count == 1

    assert result.tables.accepted_chains.num_rows == 1
    assert result.tables.rejected_chains.num_rows == 0
    assert result.tables.non_candidate_chains.num_rows == 0
    assert result.tables.dirty_residues.num_rows == 0
    assert result.tables.processing_errors.num_rows == 1

    error = result.tables.processing_errors.to_pylist()[0]
    assert error["processing_stage"] == "source_download_verify"
    assert error["model_id"] is None
    assert error["label_chain_id"] is None


def test_process_manifest_batch_chain_error_keeps_source_successful() -> None:
    from pdbclean.quality_runner import (
        SourceQualityResult,
        process_manifest_batch,
    )

    def source_processor(manifest_row, **kwargs):
        return SourceQualityResult(
            pdb_id="test",
            parsed_silver_chain_count=2,
            selected_silver_chain_count=1,
            candidate_entry_count=1,
            candidate_chain_count=1,
            processing_errors=(
                {
                    "snapshot": "20260101",
                    "pdb_id": "test",
                    "model_id": 1,
                    "label_chain_id": "A",
                    "processing_stage": "quality_cleaning",
                    "error_type": "RuntimeError",
                    "error_message": "synthetic chain failure",
                    "source_mmcif_key": "test.cif.gz",
                    "source_etag": "etag-test",
                    "pipeline_git_commit": "deadbeef",
                },
            ),
            source_failed=False,
        )

    result = process_manifest_batch(
        [{"pdb_id": "TEST"}],
        bucket_url="https://example.invalid",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
        source_processor=source_processor,
    )

    assert result.input_source_object_count == 1
    assert result.successful_source_object_count == 1
    assert result.failed_source_object_count == 0

    assert result.parsed_silver_chain_count == 2
    assert result.selected_silver_chain_count == 1
    assert result.candidate_entry_count == 1
    assert result.candidate_chain_count == 1

    assert result.tables.accepted_chains.num_rows == 0
    assert result.tables.rejected_chains.num_rows == 0
    assert result.tables.non_candidate_chains.num_rows == 0
    assert result.tables.processing_errors.num_rows == 1

    error = result.tables.processing_errors.to_pylist()[0]
    assert error["model_id"] == 1
    assert error["label_chain_id"] == "A"
    assert error["processing_stage"] == "quality_cleaning"


def test_process_manifest_batch_empty_partition() -> None:
    from pdbclean.quality_runner import process_manifest_batch
    from pdbclean.schemas import (
        GOLD_ACCEPTED_CHAIN_SCHEMA,
        GOLD_DIRTY_RESIDUE_SCHEMA,
        GOLD_NON_CANDIDATE_CHAIN_SCHEMA,
        GOLD_PROCESSING_ERROR_SCHEMA,
        GOLD_REJECTED_CHAIN_SCHEMA,
    )

    result = process_manifest_batch(
        [],
        bucket_url="https://example.invalid",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
    )

    assert result.input_source_object_count == 0
    assert result.successful_source_object_count == 0
    assert result.failed_source_object_count == 0
    assert result.parsed_silver_chain_count == 0
    assert result.selected_silver_chain_count == 0
    assert result.candidate_entry_count == 0
    assert result.candidate_chain_count == 0

    assert result.tables.accepted_chains.schema == GOLD_ACCEPTED_CHAIN_SCHEMA
    assert result.tables.rejected_chains.schema == GOLD_REJECTED_CHAIN_SCHEMA
    assert (
        result.tables.non_candidate_chains.schema
        == GOLD_NON_CANDIDATE_CHAIN_SCHEMA
    )
    assert result.tables.dirty_residues.schema == GOLD_DIRTY_RESIDUE_SCHEMA
    assert (
        result.tables.processing_errors.schema
        == GOLD_PROCESSING_ERROR_SCHEMA
    )

    assert result.tables.accepted_chains.num_rows == 0
    assert result.tables.rejected_chains.num_rows == 0
    assert result.tables.non_candidate_chains.num_rows == 0
    assert result.tables.dirty_residues.num_rows == 0
    assert result.tables.processing_errors.num_rows == 0


def _valid_empty_batch():
    from pdbclean.gold import gold_records_to_tables
    from pdbclean.quality_runner import QualityBatchResult

    return QualityBatchResult(
        input_source_object_count=0,
        successful_source_object_count=0,
        failed_source_object_count=0,
        parsed_silver_chain_count=0,
        selected_silver_chain_count=0,
        candidate_entry_count=0,
        candidate_chain_count=0,
        tables=gold_records_to_tables([]),
    )


def test_publish_quality_batch_writes_shards_before_summary(tmp_path) -> None:
    from pdbclean.quality_runner import publish_quality_batch

    calls = []

    def shard_writer(tables, output_root, task_id):
        calls.append(("shards", str(task_id)))
        return {
            "accepted": tmp_path / "accepted.parquet",
            "rejected": tmp_path / "rejected.parquet",
            "non_candidates": tmp_path / "non_candidates.parquet",
            "dirty_residues": tmp_path / "dirty_residues.parquet",
            "errors": tmp_path / "errors.parquet",
        }

    def summary_writer(summary, output_root):
        calls.append(("summary", summary["task_id"]))
        return tmp_path / "summary.json"

    publication = publish_quality_batch(
        _valid_empty_batch(),
        output_root=tmp_path,
        task_id="17",
        snapshot="20260101",
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
        started_at_utc="2026-08-09T20:00:00Z",
        started_perf_counter=10.0,
        utc_now=lambda: "2026-08-09T20:00:01Z",
        perf_counter=lambda: 11.0,
        peak_memory_reader=lambda: 1024,
        shard_writer=shard_writer,
        summary_writer=summary_writer,
    )

    assert calls == [
        ("shards", "17"),
        ("summary", "17"),
    ]
    assert publication.summary["source_object_accounting_valid"] is True
    assert publication.summary["selected_chain_accounting_valid"] is True
    assert publication.summary_path == tmp_path / "summary.json"


def test_publish_quality_batch_rejects_bad_accounting_before_writes(
    tmp_path,
) -> None:
    from pdbclean.gold import gold_records_to_tables
    from pdbclean.quality_runner import (
        QualityBatchResult,
        QualityRunnerError,
        publish_quality_batch,
    )

    batch = QualityBatchResult(
        input_source_object_count=1,
        successful_source_object_count=0,
        failed_source_object_count=0,
        parsed_silver_chain_count=0,
        selected_silver_chain_count=0,
        candidate_entry_count=0,
        candidate_chain_count=0,
        tables=gold_records_to_tables([]),
    )

    calls = []

    def shard_writer(*args, **kwargs):
        calls.append("shards")
        raise AssertionError("Shard writer must not run")

    def summary_writer(*args, **kwargs):
        calls.append("summary")
        raise AssertionError("Summary writer must not run")

    with pytest.raises(
        QualityRunnerError,
        match="source-object accounting failed",
    ):
        publish_quality_batch(
            batch,
            output_root=tmp_path,
            task_id="18",
            snapshot="20260101",
            cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
            pipeline_git_commit="deadbeef",
            started_at_utc="2026-08-09T20:00:00Z",
            started_perf_counter=10.0,
            utc_now=lambda: "2026-08-09T20:00:01Z",
            perf_counter=lambda: 11.0,
            peak_memory_reader=lambda: 1024,
            shard_writer=shard_writer,
            summary_writer=summary_writer,
        )

    assert calls == []


def test_publish_quality_batch_summary_failure_leaves_no_completion_marker(
    tmp_path,
) -> None:
    from pdbclean.quality_runner import publish_quality_batch

    calls = []

    def shard_writer(tables, output_root, task_id):
        calls.append("shards")
        return {
            "accepted": tmp_path / "accepted" / "task_19.parquet",
            "rejected": tmp_path / "rejected" / "task_19.parquet",
            "non_candidates": (
                tmp_path / "non_candidates" / "task_19.parquet"
            ),
            "dirty_residues": (
                tmp_path / "dirty_residues" / "task_19.parquet"
            ),
            "errors": tmp_path / "errors" / "task_19.parquet",
        }

    def summary_writer(summary, output_root):
        calls.append("summary")
        raise OSError("synthetic summary publication failure")

    with pytest.raises(
        OSError,
        match="synthetic summary publication failure",
    ):
        publish_quality_batch(
            _valid_empty_batch(),
            output_root=tmp_path,
            task_id="19",
            snapshot="20260101",
            cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
            pipeline_git_commit="deadbeef",
            started_at_utc="2026-08-09T20:00:00Z",
            started_perf_counter=10.0,
            utc_now=lambda: "2026-08-09T20:00:01Z",
            perf_counter=lambda: 11.0,
            peak_memory_reader=lambda: 1024,
            shard_writer=shard_writer,
            summary_writer=summary_writer,
        )

    assert calls == ["shards", "summary"]
    assert not (tmp_path / "summaries" / "task_19.json").exists()


def test_publish_quality_batch_writes_real_task_outputs(tmp_path) -> None:
    from pdbclean.quality_runner import publish_quality_batch

    publication = publish_quality_batch(
        _valid_empty_batch(),
        output_root=tmp_path,
        task_id="20",
        snapshot="20260101",
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
        started_at_utc="2026-08-09T20:00:00Z",
        started_perf_counter=10.0,
        utc_now=lambda: "2026-08-09T20:00:01Z",
        perf_counter=lambda: 11.0,
        peak_memory_reader=lambda: 1024,
    )

    assert set(publication.shard_paths) == {
        "accepted",
        "rejected",
        "non_candidates",
        "dirty_residues",
        "errors",
    }

    for path in publication.shard_paths.values():
        assert path.exists()

    assert publication.summary_path == (
        tmp_path / "summaries" / "task_20.json"
    )
    assert publication.summary_path.exists()
    assert not list(tmp_path.rglob("*.tmp"))


def test_utc_now_text_is_explicit_utc_timestamp() -> None:
    from datetime import datetime

    from pdbclean.quality_runner import _utc_now_text

    value = _utc_now_text()

    assert value.endswith("Z")

    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    assert parsed.utcoffset().total_seconds() == 0


def test_slurm_environment_reads_ids_without_defaults() -> None:
    from pdbclean.quality_runner import _slurm_environment

    assert _slurm_environment({}) == (None, None)

    assert _slurm_environment(
        {
            "SLURM_JOB_ID": "12345",
            "SLURM_ARRAY_TASK_ID": "7",
        }
    ) == ("12345", "7")


def test_linux_process_peak_memory_bytes_is_non_negative() -> None:
    from pdbclean.quality_runner import _linux_process_peak_memory_bytes

    value = _linux_process_peak_memory_bytes()

    assert value is None or (
        isinstance(value, int)
        and value >= 0
    )



def test_publish_quality_batch_captures_completion_after_shards(
    tmp_path,
) -> None:
    from pdbclean.quality_runner import publish_quality_batch

    calls = []

    def shard_writer(tables, output_root, task_id):
        calls.append("shards")
        return {
            "accepted": tmp_path / "accepted.parquet",
            "rejected": tmp_path / "rejected.parquet",
            "non_candidates": tmp_path / "non_candidates.parquet",
            "dirty_residues": tmp_path / "dirty_residues.parquet",
            "errors": tmp_path / "errors.parquet",
        }

    def utc_now():
        calls.append("completed_at")
        return "2026-08-09T20:00:05Z"

    def perf_counter():
        calls.append("runtime")
        return 15.0

    def peak_memory_reader():
        calls.append("peak_memory")
        return 4096

    def summary_writer(summary, output_root):
        calls.append("summary")
        return tmp_path / "summary.json"

    publication = publish_quality_batch(
        _valid_empty_batch(),
        output_root=tmp_path,
        task_id="timing",
        snapshot="20260101",
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
        started_at_utc="2026-08-09T20:00:00Z",
        started_perf_counter=10.0,
        shard_writer=shard_writer,
        summary_writer=summary_writer,
        utc_now=utc_now,
        perf_counter=perf_counter,
        peak_memory_reader=peak_memory_reader,
    )

    assert calls == [
        "shards",
        "completed_at",
        "runtime",
        "peak_memory",
        "summary",
    ]

    assert publication.summary["started_at_utc"] == (
        "2026-08-09T20:00:00Z"
    )
    assert publication.summary["completed_at_utc"] == (
        "2026-08-09T20:00:05Z"
    )
    assert publication.summary["runtime_seconds"] == 5.0
    assert publication.summary["peak_memory_bytes"] == 4096


def test_execute_quality_task_orchestrates_batch_then_publication(
    tmp_path,
) -> None:
    from pdbclean.quality_runner import (
        QualityTaskPublication,
        execute_quality_task,
    )

    calls = []
    batch = _valid_empty_batch()

    def utc_now():
        calls.append("utc_start")
        return "2026-08-10T00:40:00Z"

    def perf_counter():
        calls.append("perf_start")
        return 100.0

    def batch_processor(manifest_rows, **kwargs):
        calls.append("batch")
        assert list(manifest_rows) == [{"pdb_id": "TEST"}]
        assert kwargs == {
            "bucket_url": "https://example.invalid",
            "selection_config": {
                "models": {
                    "policy": "first_model",
                    "model_id": 1,
                }
            },
            "cleaning_protocol": "Protocol_3.2_BRI_v1.2.2",
            "pipeline_git_commit": "deadbeef",
            "minimum_backbone_distance_angstrom": 0.01,
            "timeout_seconds": 37,
            "max_retries": 4,
            "download_concurrency": 3,
        }
        return batch

    expected = QualityTaskPublication(
        shard_paths={},
        summary_path=tmp_path / "summary.json",
        summary={"task_id": "7"},
    )

    def publisher(received_batch, **kwargs):
        calls.append("publisher")

        assert received_batch is batch
        assert kwargs["output_root"] == tmp_path
        assert kwargs["task_id"] == "7"
        assert kwargs["snapshot"] == "20260101"
        assert kwargs["cleaning_protocol"] == "Protocol_3.2_BRI_v1.2.2"
        assert kwargs["pipeline_git_commit"] == "deadbeef"

        assert kwargs["started_at_utc"] == "2026-08-10T00:40:00Z"
        assert kwargs["started_perf_counter"] == 100.0

        assert kwargs["slurm_job_id"] == "12345"
        assert kwargs["slurm_array_task_id"] == "7"

        # The same clocks are forwarded so publication can capture
        # completion only after Parquet shards have been written.
        assert kwargs["utc_now"] is utc_now
        assert kwargs["perf_counter"] is perf_counter

        return expected

    result = execute_quality_task(
        [{"pdb_id": "TEST"}],
        output_root=tmp_path,
        task_id="7",
        snapshot="20260101",
        bucket_url="https://example.invalid",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
        timeout_seconds=37,
        max_retries=4,
        download_concurrency=3,
        environ={
            "SLURM_JOB_ID": "12345",
            "SLURM_ARRAY_TASK_ID": "7",
        },
        batch_processor=batch_processor,
        publisher=publisher,
        utc_now=utc_now,
        perf_counter=perf_counter,
    )

    assert result is expected
    assert calls == [
        "utc_start",
        "perf_start",
        "batch",
        "publisher",
    ]


def test_execute_quality_task_does_not_publish_when_batch_fails(
    tmp_path,
) -> None:
    from pdbclean.quality_runner import execute_quality_task

    calls = []

    def batch_processor(manifest_rows, **kwargs):
        calls.append("batch")
        raise RuntimeError("synthetic batch failure")

    def publisher(*args, **kwargs):
        calls.append("publisher")
        raise AssertionError("Publisher must not run")

    with pytest.raises(RuntimeError, match="synthetic batch failure"):
        execute_quality_task(
            [],
            output_root=tmp_path,
            task_id="8",
            snapshot="20260101",
            bucket_url="https://example.invalid",
            selection_config={
                "models": {
                    "policy": "first_model",
                    "model_id": 1,
                }
            },
            cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
            pipeline_git_commit="deadbeef",
            environ={},
            batch_processor=batch_processor,
            publisher=publisher,
            utc_now=lambda: "2026-08-10T00:40:00Z",
            perf_counter=lambda: 100.0,
        )

    assert calls == ["batch"]


def test_quality_stage_output_root_is_dynamic() -> None:
    from pathlib import Path

    from pdbclean.quality_runner import quality_stage_output_root

    observed = quality_stage_output_root(
        "outputs/pdbclean",
        snapshot="20310415",
        protocol_version="protocol-next-v2",
    )

    assert observed == Path(
        "outputs/pdbclean/20310415/protocol-next-v2/quality"
    )


@pytest.mark.parametrize(
    "snapshot",
    ["", ".", "..", "../escape", "a/b", "a\\b"],
)
def test_quality_stage_output_root_rejects_unsafe_snapshot(
    snapshot: str,
) -> None:
    from pdbclean.quality_runner import (
        QualityRunnerError,
        quality_stage_output_root,
    )

    with pytest.raises(
        QualityRunnerError,
        match="Unsafe quality-stage snapshot",
    ):
        quality_stage_output_root(
            "outputs/pdbclean",
            snapshot=snapshot,
            protocol_version="protocol-v1",
        )


@pytest.mark.parametrize(
    "protocol_version",
    ["", ".", "..", "../escape", "a/b", "a\\b"],
)
def test_quality_stage_output_root_rejects_unsafe_protocol(
    protocol_version: str,
) -> None:
    from pdbclean.quality_runner import (
        QualityRunnerError,
        quality_stage_output_root,
    )

    with pytest.raises(
        QualityRunnerError,
        match="Unsafe quality-stage protocol version",
    ):
        quality_stage_output_root(
            "outputs/pdbclean",
            snapshot="20310415",
            protocol_version=protocol_version,
        )


def test_process_manifest_source_retries_transport_failure_then_succeeds() -> None:
    from pdbclean.quality_runner import process_manifest_source
    from pdbclean.snapshot import SnapshotTransportError

    compressed = _multimodel_cif_bytes()
    attempts = 0

    def downloader(**kwargs):
        nonlocal attempts
        attempts += 1

        if attempts < 3:
            raise SnapshotTransportError(
                f"synthetic transient failure {attempts}"
            )

        return compressed

    result = process_manifest_source(
        _manifest_row(),
        bucket_url="https://example.invalid",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
        max_retries=3,
        downloader=downloader,
    )

    assert attempts == 3
    assert result.source_failed is False
    assert len(result.gold_records) == 1


def test_process_manifest_source_exhausts_transport_retries() -> None:
    from pdbclean.quality_runner import process_manifest_source
    from pdbclean.snapshot import SnapshotTransportError

    attempts = 0

    def downloader(**kwargs):
        nonlocal attempts
        attempts += 1
        raise SnapshotTransportError("synthetic network outage")

    result = process_manifest_source(
        _manifest_row(),
        bucket_url="https://example.invalid",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
        max_retries=2,
        downloader=downloader,
    )

    # One initial attempt plus two retries.
    assert attempts == 3
    assert result.source_failed is True
    assert len(result.processing_errors) == 1
    assert (
        result.processing_errors[0]["error_type"]
        == "SnapshotTransportError"
    )


def test_process_manifest_source_does_not_retry_verification_failure() -> None:
    from pdbclean.quality_runner import process_manifest_source
    from pdbclean.snapshot import SnapshotVerificationError

    attempts = 0

    def downloader(**kwargs):
        nonlocal attempts
        attempts += 1
        raise SnapshotVerificationError(
            "synthetic ETag mismatch"
        )

    result = process_manifest_source(
        _manifest_row(),
        bucket_url="https://example.invalid",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
        max_retries=3,
        downloader=downloader,
    )

    assert attempts == 1
    assert result.source_failed is True
    assert (
        result.processing_errors[0]["error_type"]
        == "SnapshotVerificationError"
    )


@pytest.mark.parametrize("max_retries", [-1, True, 1.5, "3"])
def test_process_manifest_source_rejects_invalid_max_retries(
    max_retries,
) -> None:
    from pdbclean.quality_runner import (
        QualityRunnerError,
        process_manifest_source,
    )

    with pytest.raises(
        QualityRunnerError,
        match="max_retries must be a non-negative integer",
    ):
        process_manifest_source(
            _manifest_row(),
            bucket_url="https://example.invalid",
            selection_config={
                "models": {
                    "policy": "first_model",
                    "model_id": 1,
                }
            },
            cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
            pipeline_git_commit="deadbeef",
            max_retries=max_retries,
        )


def test_process_manifest_batch_forwards_max_retries() -> None:
    from pdbclean.quality_runner import (
        SourceQualityResult,
        process_manifest_batch,
    )

    observed = []

    def source_processor(manifest_row, **kwargs):
        observed.append(kwargs["max_retries"])

        return SourceQualityResult(
            pdb_id=manifest_row["pdb_id"].lower(),
            parsed_silver_chain_count=0,
            selected_silver_chain_count=0,
            candidate_entry_count=0,
            candidate_chain_count=0,
            source_failed=False,
        )

    result = process_manifest_batch(
        [{"pdb_id": "TEST"}],
        bucket_url="https://example.invalid",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
        max_retries=7,
        download_concurrency=1,
        source_processor=source_processor,
    )

    assert observed == [7]
    assert result.input_source_object_count == 1
    assert result.successful_source_object_count == 1


def test_process_manifest_batch_concurrency_preserves_manifest_order() -> None:
    import threading

    from pdbclean.quality_runner import (
        SourceQualityResult,
        process_manifest_batch,
    )

    second_started = threading.Event()

    rows = [
        {"pdb_id": "FIRST"},
        {"pdb_id": "SECOND"},
    ]

    def source_processor(manifest_row, **kwargs):
        pdb_id = manifest_row["pdb_id"].lower()

        assert kwargs["max_retries"] == 2

        if pdb_id == "first":
            if not second_started.wait(timeout=2.0):
                raise AssertionError(
                    "Second source never started concurrently"
                )
        else:
            second_started.set()

        return SourceQualityResult(
            pdb_id=pdb_id,
            parsed_silver_chain_count=0,
            selected_silver_chain_count=0,
            candidate_entry_count=0,
            candidate_chain_count=0,
            processing_errors=(
                {
                    "snapshot": "20260101",
                    "pdb_id": pdb_id,
                    "model_id": None,
                    "label_chain_id": None,
                    "processing_stage": "source_download_verify",
                    "error_type": "SyntheticError",
                    "error_message": pdb_id,
                    "source_mmcif_key": f"{pdb_id}.cif.gz",
                    "source_etag": f"etag-{pdb_id}",
                    "pipeline_git_commit": "deadbeef",
                },
            ),
            source_failed=True,
        )

    result = process_manifest_batch(
        rows,
        bucket_url="https://example.invalid",
        selection_config={
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
        cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
        pipeline_git_commit="deadbeef",
        max_retries=2,
        download_concurrency=2,
        source_processor=source_processor,
    )

    # SECOND is allowed to finish first internally, but publication order
    # must still follow the immutable manifest.
    assert [
        row["pdb_id"]
        for row in result.tables.processing_errors.to_pylist()
    ] == ["first", "second"]

    assert result.input_source_object_count == 2
    assert result.successful_source_object_count == 0
    assert result.failed_source_object_count == 2


@pytest.mark.parametrize(
    "download_concurrency",
    [0, -1, True, 1.5, "4"],
)
def test_process_manifest_batch_rejects_invalid_download_concurrency(
    download_concurrency,
) -> None:
    from pdbclean.quality_runner import (
        QualityRunnerError,
        process_manifest_batch,
    )

    with pytest.raises(
        QualityRunnerError,
        match="download_concurrency must be a positive integer",
    ):
        process_manifest_batch(
            [],
            bucket_url="https://example.invalid",
            selection_config={
                "models": {
                    "policy": "first_model",
                    "model_id": 1,
                }
            },
            cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
            pipeline_git_commit="deadbeef",
            download_concurrency=download_concurrency,
        )


@pytest.mark.parametrize(
    "max_retries",
    [-1, True, 1.5, "3"],
)
def test_process_manifest_batch_rejects_invalid_max_retries(
    max_retries,
) -> None:
    from pdbclean.quality_runner import (
        QualityRunnerError,
        process_manifest_batch,
    )

    with pytest.raises(
        QualityRunnerError,
        match="max_retries must be a non-negative integer",
    ):
        process_manifest_batch(
            [],
            bucket_url="https://example.invalid",
            selection_config={
                "models": {
                    "policy": "first_model",
                    "model_id": 1,
                }
            },
            cleaning_protocol="Protocol_3.2_BRI_v1.2.2",
            pipeline_git_commit="deadbeef",
            max_retries=max_retries,
        )
