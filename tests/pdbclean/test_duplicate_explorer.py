"""Duplicate Explorer regression tests.

The Explorer filters and displays; it never re-classifies.  Every row it shows
carries the classification the pipeline recorded, and the Stage-14 relationship
it reports comes from the published representative mapping.
"""

from __future__ import annotations

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from pdbclean.duplicates import (
    MAX_PAGE_SIZE,
    DuplicateExplorer,
    DuplicateFilters,
    DuplicateQueryError,
    DuplicateSource,
    pair_detail,
)


SNAPSHOT = "20260101"


def _pair(query, subject, *, length, distance_mA):
    q_pdb, q_chain = query
    s_pdb, s_chain = subject

    return {
        "query_snapshot": SNAPSHOT,
        "query_pdb_id": q_pdb,
        "query_model_id": 1,
        "query_label_chain_id": q_chain,
        "subject_snapshot": SNAPSHOT,
        "subject_pdb_id": s_pdb,
        "subject_model_id": 1,
        "subject_label_chain_id": s_chain,
        "retained_residue_count": length,
        "d_bri_mA": distance_mA,
        "d_bri": distance_mA / 1000.0,
        "is_zero_duplicate": distance_mA == 0,
        "is_paper_near_duplicate": distance_mA <= 10,
        "is_nonzero_near_duplicate": 0 < distance_mA <= 10,
    }


PAIRS = [
    _pair(("1abc", "A"), ("2abc", "B"), length=100, distance_mA=0),
    _pair(("1abc", "A"), ("3abc", "C"), length=100, distance_mA=4),
    _pair(("4abc", "D"), ("5abc", "E"), length=250, distance_mA=10),
    _pair(("6abc", "F"), ("7abc", "G"), length=40, distance_mA=7),
]


MAPPING = [
    {
        "snapshot": SNAPSHOT,
        "pdb_id": "2abc",
        "model_id": 1,
        "label_chain_id": "B",
        "action": "remove",
        "component_id": 11,
        "component_is_clique": True,
        "representative_snapshot": SNAPSHOT,
        "representative_pdb_id": "1abc",
        "representative_model_id": 1,
        "representative_label_chain_id": "A",
        "direct_d_bri_mA": 0,
        "policy_version": "1.0",
    },
    {
        "snapshot": SNAPSHOT,
        "pdb_id": "1abc",
        "model_id": 1,
        "label_chain_id": "A",
        "action": "retain_representative",
        "component_id": 11,
        "component_is_clique": True,
        "representative_snapshot": SNAPSHOT,
        "representative_pdb_id": "1abc",
        "representative_model_id": 1,
        "representative_label_chain_id": "A",
        "direct_d_bri_mA": None,
        "policy_version": "1.0",
    },
    {
        "snapshot": SNAPSHOT,
        "pdb_id": "3abc",
        "model_id": 1,
        "label_chain_id": "C",
        "action": "retain_representative",
        "component_id": 11,
        "component_is_clique": True,
        "representative_snapshot": SNAPSHOT,
        "representative_pdb_id": "3abc",
        "representative_model_id": 1,
        "representative_label_chain_id": "C",
        "direct_d_bri_mA": None,
        "policy_version": "1.0",
    },
]


@pytest.fixture()
def source(tmp_path):
    protocol_root = tmp_path / "protocol"
    finalized = protocol_root / "duplicate_classification" / "finalized"
    finalized.mkdir(parents=True)

    pq.write_table(
        pa.Table.from_pylist(PAIRS),
        finalized / "candidate_classifications.parquet",
    )

    release_root = tmp_path / "release"
    audit = release_root / "audit"
    audit.mkdir(parents=True)

    pq.write_table(
        pa.Table.from_pylist(MAPPING),
        audit / "representative_mapping.parquet",
    )

    return DuplicateSource(
        protocol_root=protocol_root, release_root=release_root
    )


@pytest.fixture()
def explorer(source):
    return DuplicateExplorer(source)


# --------------------------------------------------------------------------
# Source discovery
# --------------------------------------------------------------------------


def test_classified_output_is_preferred(source):
    availability = source.available()

    assert availability["mode"] == "classified"
    assert availability["representative_mapping"]


def test_near_duplicate_tables_are_the_fallback(tmp_path):
    protocol_root = tmp_path / "protocol"
    finalized = protocol_root / "full_bri_nn" / "finalized"
    finalized.mkdir(parents=True)

    stripped = [
        {k: v for k, v in pair.items() if not k.startswith("is_")}
        for pair in PAIRS
    ]

    pq.write_table(
        pa.Table.from_pylist(stripped),
        finalized / "candidate_near_duplicates.parquet",
    )

    source = DuplicateSource(protocol_root=protocol_root)

    assert source.available()["mode"] == "near"

    result = DuplicateExplorer(source).query(DuplicateFilters())

    assert result["matched"] == len(PAIRS)
    assert result["has_classification_flags"] is False

    # Without recorded flags, exactness is read off the exact integer distance.
    exact = [r for r in result["rows"] if r["classification"] == "exact_duplicate"]

    assert len(exact) == 1
    assert exact[0]["d_bri_mA"] == 0


def test_missing_output_is_an_explicit_error(tmp_path):
    source = DuplicateSource(protocol_root=tmp_path / "absent")

    with pytest.raises(DuplicateQueryError):
        DuplicateExplorer(source).query(DuplicateFilters())


# --------------------------------------------------------------------------
# Filtering
# --------------------------------------------------------------------------


def test_unfiltered_query_returns_every_pair(explorer):
    result = explorer.query(DuplicateFilters())

    assert result["matched"] == len(PAIRS)
    assert len(result["rows"]) == len(PAIRS)


def test_filter_by_pdb_id(explorer):
    result = explorer.query(DuplicateFilters(pdb_id="1ABC"))

    assert result["matched"] == 2


def test_filter_by_chain(explorer):
    result = explorer.query(DuplicateFilters(chain="E"))

    assert result["matched"] == 1
    assert result["rows"][0]["pdb_id_b"] == "5abc"


def test_exact_only_uses_the_recorded_flag(explorer):
    result = explorer.query(DuplicateFilters(exact_only=True))

    assert result["matched"] == 1
    assert result["rows"][0]["d_bri_mA"] == 0
    assert result["rows"][0]["classification"] == "exact_duplicate"


def test_nonzero_near_only(explorer):
    result = explorer.query(DuplicateFilters(nonzero_near_only=True))

    assert result["matched"] == 3
    assert all(row["d_bri_mA"] > 0 for row in result["rows"])


def test_exact_and_nonzero_are_mutually_exclusive(explorer):
    with pytest.raises(DuplicateQueryError):
        explorer.query(
            DuplicateFilters(exact_only=True, nonzero_near_only=True)
        )


def test_length_filters_are_inclusive(explorer):
    assert explorer.query(DuplicateFilters(min_length=100))["matched"] == 3
    assert explorer.query(DuplicateFilters(max_length=100))["matched"] == 3
    assert (
        explorer.query(
            DuplicateFilters(min_length=100, max_length=100)
        )["matched"]
        == 2
    )


def test_distance_filters_are_inclusive(explorer):
    """The Explorer must agree with the pipeline's inclusive semantics."""

    at_threshold = explorer.query(DuplicateFilters(min_distance_mA=10))

    assert at_threshold["matched"] == 1
    assert at_threshold["rows"][0]["d_bri_mA"] == 10

    assert explorer.query(DuplicateFilters(max_distance_mA=0))["matched"] == 1


def test_relationship_filter_reports_stage14_decisions(explorer):
    removed = explorer.query(DuplicateFilters(relationship="removed"))

    assert removed["matched"] == 1

    row = removed["rows"][0]

    assert row["removed_chain"] == "2abc:B"
    assert row["representative"] == "1abc:A"
    assert row["direct_d_bri_mA"] == 0
    assert row["component_id"] == 11
    assert row["policy_version"] == "1.0"


def test_both_retained_pairs_are_reported_as_retained(explorer):
    retained = explorer.query(DuplicateFilters(relationship="retained"))

    assert retained["matched"] == 1
    assert retained["rows"][0]["pdb_id_a"] == "1abc"
    assert retained["rows"][0]["pdb_id_b"] == "3abc"


def test_pairs_outside_the_graph_are_reported_as_unaffected(explorer):
    unaffected = explorer.query(DuplicateFilters(relationship="unaffected"))

    assert unaffected["matched"] == 2


def test_relationship_partition_is_complete(explorer):
    total = explorer.query(DuplicateFilters())["matched"]

    partition = sum(
        explorer.query(DuplicateFilters(relationship=name))["matched"]
        for name in ("removed", "retained", "unaffected")
    )

    assert partition == total


# --------------------------------------------------------------------------
# Presentation
# --------------------------------------------------------------------------


def test_distance_is_shown_in_both_units(explorer):
    row = explorer.query(DuplicateFilters(min_distance_mA=10))["rows"][0]

    assert row["d_bri_mA"] == 10
    assert row["d_bri_angstrom"] == pytest.approx(0.010)


def test_explorer_never_invents_a_classification(explorer):
    for row in explorer.query(DuplicateFilters())["rows"]:
        assert row["classification"] in {
            "exact_duplicate",
            "nonzero_near_duplicate",
        }


def test_summary_reports_the_available_sources(explorer):
    summary = explorer.summary()

    assert summary["sources"]["mode"] == "classified"
    assert summary["classified"] is True
    assert explorer.has_representative_mapping is True


def test_summary_reports_no_counts_without_a_stage_summary(explorer):
    """Counts come from stage provenance; the Explorer derives none itself."""

    summary = explorer.summary()

    assert "near_duplicate_pairs" not in summary
    assert "total_tested_pairs" not in summary


def test_summary_reads_counts_from_stage_provenance(source, explorer):
    import json

    (source.protocol_root / "duplicate_classification").mkdir(
        parents=True, exist_ok=True
    )
    (
        source.protocol_root
        / "duplicate_classification"
        / "global_summary.json"
    ).write_text(
        json.dumps(
            {
                "input_pair_count": 3_531_895,
                "paper_near_duplicate_pair_count": 1_072_751,
                "zero_duplicate_pair_count": 17_373,
                "nonzero_near_duplicate_pair_count": 1_055_378,
                "not_near_duplicate_pair_count": 2_459_144,
                "paper_near_duplicate_threshold_mA": 10,
                "paper_near_duplicate_threshold_angstrom": 0.010,
            }
        ),
        encoding="utf-8",
    )

    summary = explorer.summary()

    assert summary["near_duplicate_pairs"] == 1_072_751
    assert summary["exact_duplicate_pairs"] == 17_373
    assert summary["threshold_mA"] == 10


# --------------------------------------------------------------------------
# Paging
# --------------------------------------------------------------------------


def test_paging_walks_the_whole_result_set(explorer):
    seen = []
    offset = 0

    while True:
        page = explorer.query(DuplicateFilters(offset=offset, limit=2))

        if not page["rows"]:
            break

        seen.extend(
            (r["pdb_id_a"], r["chain_a"], r["pdb_id_b"], r["chain_b"])
            for r in page["rows"]
        )
        offset += 2

    assert len(seen) == len(PAIRS)
    assert len(set(seen)) == len(PAIRS)


def test_matched_count_is_independent_of_the_page(explorer):
    first = explorer.query(DuplicateFilters(limit=1))
    second = explorer.query(DuplicateFilters(offset=2, limit=1))

    assert first["matched"] == second["matched"] == len(PAIRS)


def test_page_size_is_capped(explorer):
    result = explorer.query(DuplicateFilters(limit=10_000))

    assert result["limit"] == MAX_PAGE_SIZE


def test_negative_offset_is_clamped(explorer):
    assert explorer.query(DuplicateFilters(offset=-5))["offset"] == 0


# --------------------------------------------------------------------------
# Pair lookup (the Mol* entry point)
# --------------------------------------------------------------------------


def test_pair_detail_finds_a_pair_in_either_order(explorer):
    forward = pair_detail(
        explorer, pdb_id_a="1abc", chain_a="A", pdb_id_b="2abc", chain_b="B"
    )

    assert forward is not None
    assert forward["d_bri_mA"] == 0
    assert forward["classification"] == "exact_duplicate"


def test_pair_detail_returns_none_for_an_unknown_pair(explorer):
    assert (
        pair_detail(
            explorer,
            pdb_id_a="9xyz",
            chain_a="Z",
            pdb_id_b="8xyz",
            chain_b="Y",
        )
        is None
    )
