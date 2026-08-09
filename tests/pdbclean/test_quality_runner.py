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


def test_select_configured_model_chains_rejects_unsupported_policy() -> None:
    with pytest.raises(
        QualityRunnerError,
        match="Unsupported model-selection policy",
    ):
        select_configured_model_chains(
            [_chain(1, "A")],
            {
                "models": {
                    "policy": "all_models",
                    "model_id": 1,
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

    def fail_cleaning(chain):
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
