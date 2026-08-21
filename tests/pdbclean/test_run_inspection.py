"""Historical-run inspection: canonical identity, ordering and immutability.

A historical run must read like a complete audit record:

* stages appear in canonical scientific order, never alphabetical;
* prerequisites are lettered, never numbered as scientific stages;
* canonical Stage 1-14 identities are preserved, including the shared-producer
  pairs (3/4 and 8/9) and the Stage-14 subdivisions;
* Stages 11-13 are shown with an honest status, never silently omitted and
  never mistakable for the Stage-14 deletion relation;
* inspecting a run mutates absolutely nothing.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from pdbclean.artefacts import (
    MAX_PREVIEW_ROWS,
    ArtefactError,
    describe,
    list_directory,
    preview,
    preview_kind,
)
from pdbclean.pipeline import plan_pipeline, record_plan_in_provenance
from pdbclean.run_inspection import (
    NOT_RECORDED,
    STATUS_OPTIONAL,
    duplicate_navigation,
    run_timeline,
    stage_detail,
)
from pdbclean.run_provenance import RunProvenance
from pdbclean.runconfig import resolve_run_config
from pdbclean.stage_registry import canonical_timeline


REPO = Path(__file__).resolve().parents[2]
PROFILE = REPO / "config" / "pdbclean" / "profiles" / "comp702_frozen_20260101.yaml"


@pytest.fixture()
def historical_run(tmp_path):
    """A completed-looking run whose plan reflects the frozen outputs."""

    resolved = resolve_run_config(
        config_path=PROFILE,
        overrides={"storage.run_root": str(tmp_path / "runs")},
    )

    run = RunProvenance.create(
        resolved=resolved,
        run_root=tmp_path / "runs",
        repo_root=REPO,
        snapshot={
            "snapshot_id": "20260101",
            "display": "2026-01-01",
            "selection_mode": "latest_complete",
        },
    )

    record_plan_in_provenance(
        plan_pipeline(resolved, repo_root=REPO), run
    )
    run.flush()

    return run


# --------------------------------------------------------------------------
# 1-2. Canonical order; prerequisites before Stage 1
# --------------------------------------------------------------------------


def test_timeline_is_in_canonical_order_not_alphabetical(historical_run):
    rows = run_timeline(historical_run.record)

    positions = [row["position"] for row in rows]

    assert positions == sorted(positions)
    assert positions == list(range(1, len(rows) + 1))

    labels = [row["label"] for row in rows]

    assert labels != sorted(labels), "timeline must not be alphabetical"


def test_expected_canonical_sequence(historical_run):
    labels = [row["label"] for row in run_timeline(historical_run.record)]

    assert labels == [
        "Prerequisite A",
        "Prerequisite B",
        "Prerequisite C",
        "Stage 1",
        "Stage 2",
        "Stage 3",
        "Stage 4",
        "Stage 5",
        "Stage 6",
        "Stage 7",
        "Stage 8",
        "Stage 9",
        "Stage 10",
        "Stage 11",
        "Stage 12",
        "Stage 13",
        "Stage 14 input",
        "Stage 14a",
        "Stage 14b",
        "Stage 14c",
    ]


def test_prerequisites_precede_stage_one(historical_run):
    rows = run_timeline(historical_run.record)

    prerequisites = [r for r in rows if r["role"] == "prerequisite"]
    stage_one = next(r for r in rows if r["label"] == "Stage 1")

    assert prerequisites
    assert all(r["position"] < stage_one["position"] for r in prerequisites)


# --------------------------------------------------------------------------
# 3. Prerequisites are not numbered as scientific stages
# --------------------------------------------------------------------------


def test_prerequisites_are_lettered_never_numbered(historical_run):
    rows = run_timeline(historical_run.record)

    for row in rows:
        if row["role"] != "prerequisite":
            continue

        assert row["label"].startswith("Prerequisite ")
        assert not row["label"].startswith("Stage ")


def test_no_scientific_stage_number_is_claimed_by_a_prerequisite():
    scientific = {
        entry.label
        for entry in canonical_timeline()
        if entry.role != "prerequisite"
    }
    prerequisites = {
        entry.label
        for entry in canonical_timeline()
        if entry.role == "prerequisite"
    }

    assert not scientific & prerequisites


# --------------------------------------------------------------------------
# 4-7. Canonical identities, shared producers, Stage 14 subdivisions
# --------------------------------------------------------------------------


def test_stage_one_through_fourteen_are_all_present(historical_run):
    labels = {row["label"] for row in run_timeline(historical_run.record)}

    for number in range(1, 14):
        assert f"Stage {number}" in labels, number

    assert {"Stage 14a", "Stage 14b", "Stage 14c"} <= labels


def test_stage_three_and_four_both_appear_and_share_a_producer(historical_run):
    rows = {row["label"]: row for row in run_timeline(historical_run.record)}

    three, four = rows["Stage 3"], rows["Stage 4"]

    assert three["producer"] == four["producer"] == "complete_bri"
    assert three["shared_producer"] and four["shared_producer"]
    assert three["title"] != four["title"]
    assert four["position"] == three["position"] + 1


def test_stage_eight_and_nine_both_appear_and_share_a_producer(historical_run):
    rows = {row["label"]: row for row in run_timeline(historical_run.record)}

    eight, nine = rows["Stage 8"], rows["Stage 9"]

    assert eight["producer"] == nine["producer"] == "complete_bri_nn"
    assert eight["shared_producer"] and nine["shared_producer"]
    assert nine["position"] == eight["position"] + 1


def test_stage_fourteen_subdivisions_map_to_stage_fourteen(historical_run):
    rows = {row["label"]: row for row in run_timeline(historical_run.record)}

    for label in ("Stage 14a", "Stage 14b", "Stage 14c"):
        assert rows[label]["parent"] == "Stage 14"

    assert "Stage 14" not in rows, "Stage 14 is realised by its subdivisions"


# --------------------------------------------------------------------------
# 8. Stages 11-13 are shown, and are never the deletion relation
# --------------------------------------------------------------------------


def test_investigation_stages_are_shown_not_omitted(historical_run):
    rows = {row["label"]: row for row in run_timeline(historical_run.record)}

    for label in ("Stage 11", "Stage 12", "Stage 13"):
        assert label in rows
        assert rows[label]["status"] == STATUS_OPTIONAL
        assert rows[label]["producer"] is None


def test_timeline_never_jumps_from_ten_to_fourteen(historical_run):
    labels = [row["label"] for row in run_timeline(historical_run.record)]

    ten = labels.index("Stage 10")

    assert labels[ten + 1 : ten + 4] == ["Stage 11", "Stage 12", "Stage 13"]


def test_investigation_stages_are_not_the_deletion_relation(historical_run):
    detail = stage_detail(historical_run.record, "stage_13", repo_root=REPO)

    assert "NOT the global Stage-14 deletion relation" in (
        detail["identity"]["note"]
    )
    assert detail["identity"]["producer"] == NOT_RECORDED

    graph = stage_detail(historical_run.record, "stage_14a", repo_root=REPO)

    assert graph["identity"]["producer"] == "redundancy_graph"
    assert graph["identity"]["parent_stage"] == "Stage 14"


# --------------------------------------------------------------------------
# 9. Inspection mutates nothing
# --------------------------------------------------------------------------


def test_inspection_does_not_mutate_provenance(historical_run):
    record_path = historical_run.run_dir / "run.json"
    events_path = historical_run.run_dir / "events.jsonl"

    before_record = hashlib.sha256(record_path.read_bytes()).hexdigest()
    before_events = hashlib.sha256(events_path.read_bytes()).hexdigest()

    record = json.loads(record_path.read_text(encoding="utf-8"))

    run_timeline(record)

    for entry in canonical_timeline():
        stage_detail(record, entry.key, repo_root=REPO)

    assert hashlib.sha256(record_path.read_bytes()).hexdigest() == before_record
    assert hashlib.sha256(events_path.read_bytes()).hexdigest() == before_events


def test_inspection_does_not_reresolve_the_snapshot(historical_run):
    record = historical_run.record

    before = dict(record["snapshot"])

    run_timeline(record)
    stage_detail(record, "stage_1", repo_root=REPO)

    assert record["snapshot"] == before
    assert record["snapshot"]["snapshot_id"] == "20260101"


def test_inspection_does_not_recompute_hashes(historical_run):
    record = historical_run.record

    before = record["scientific_config_sha256"]

    detail = stage_detail(record, "stage_10", repo_root=REPO)

    assert detail["execution"]["scientific_config_sha256"] == before
    assert record["scientific_config_sha256"] == before


# --------------------------------------------------------------------------
# 10-11. Detail comes from that run; missing fields are safe
# --------------------------------------------------------------------------


def test_stage_detail_reports_that_runs_own_values(historical_run):
    detail = stage_detail(historical_run.record, "stage_10", repo_root=REPO)

    assert detail["identity"]["canonical_label"] == "Stage 10"
    assert detail["configuration"]["snapshot.snapshot_id"] == "20260101"
    assert detail["configuration"][
        "duplicate_search.near_duplicate_threshold_angstrom"
    ] == 0.010
    assert detail["configuration"][
        "bri.representation_precision_angstrom"
    ] == 0.001


def test_stage_detail_exposes_every_documented_section(historical_run):
    detail = stage_detail(historical_run.record, "stage_14b", repo_root=REPO)

    for section in (
        "identity",
        "status",
        "configuration",
        "inputs",
        "outputs",
        "validation",
        "execution",
        "reuse",
        "artefacts",
    ):
        assert section in detail, section


def test_missing_fields_are_reported_not_invented(historical_run):
    detail = stage_detail(historical_run.record, "stage_1", repo_root=REPO)

    # This run never executed, so timings were never recorded.
    assert detail["execution"]["started_at"] == NOT_RECORDED
    assert detail["execution"]["runtime_seconds"] == NOT_RECORDED


def test_unknown_canonical_stage_is_an_error(historical_run):
    with pytest.raises(KeyError):
        stage_detail(historical_run.record, "stage_99", repo_root=REPO)


def test_stage_detail_names_its_sibling_identity(historical_run):
    detail = stage_detail(historical_run.record, "stage_3", repo_root=REPO)

    assert detail["identity"]["sibling_identities"] == [
        "Stage 4 — BRI numerical representation"
    ]
    assert detail["identity"]["shared_producer"] is True


# --------------------------------------------------------------------------
# 12. Previews are bounded
# --------------------------------------------------------------------------


def test_preview_row_limit_is_capped(tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq

    target = tmp_path / "rows.parquet"
    pq.write_table(
        pa.table({"n": list(range(5000))}), target
    )

    result = preview(target, limit=10_000)["preview"]

    assert result["row_count"] == 5000
    assert result["row_preview_count"] <= MAX_PREVIEW_ROWS
    assert result["truncated"] is True


def test_csv_preview_is_bounded(tmp_path):
    target = tmp_path / "rows.csv"
    target.write_text(
        "a,b\n" + "".join(f"{i},{i}\n" for i in range(500)), encoding="utf-8"
    )

    result = preview(target, limit=25)["preview"]

    assert result["columns"] == ["a", "b"]
    assert result["row_preview_count"] == 25
    assert result["truncated"] is True


def test_large_text_preview_is_truncated(tmp_path):
    from pdbclean.artefacts import MAX_TEXT_BYTES

    target = tmp_path / "big.log"
    target.write_text("x" * (MAX_TEXT_BYTES + 2048), encoding="utf-8")

    result = preview(target)["preview"]

    assert result["truncated"] is True
    assert len(result["content"]) <= MAX_TEXT_BYTES


def test_binary_artefacts_report_metadata_only(tmp_path):
    target = tmp_path / "weights.bin"
    target.write_bytes(b"\x00\x01\x02" * 100)

    assert preview_kind(target) == "metadata"

    result = preview(target)["preview"]

    assert result["kind"] == "metadata"
    assert "content" not in result


def test_json_preview_round_trips(tmp_path):
    target = tmp_path / "summary.json"
    target.write_text(json.dumps({"retained": 499770}), encoding="utf-8")

    result = preview(target)["preview"]

    assert result["kind"] == "json"
    assert result["content"]["retained"] == 499770
    assert result["truncated"] is False


def test_preview_refuses_a_missing_artefact(tmp_path):
    with pytest.raises(ArtefactError):
        preview(tmp_path / "absent.json")


def test_describe_skips_hashing_enormous_files(tmp_path, monkeypatch):
    import pdbclean.artefacts as module

    target = tmp_path / "huge.parquet"
    target.write_bytes(b"0" * 1024)

    monkeypatch.setattr(module, "MAX_HASH_BYTES", 10)

    assert module.describe(target)["sha256"] is None


def test_listing_is_bounded(tmp_path):
    for index in range(50):
        (tmp_path / f"file{index}.json").write_text("{}", encoding="utf-8")

    assert len(list_directory(tmp_path, limit=10)) == 10


def test_listing_a_missing_directory_is_empty(tmp_path):
    assert list_directory(tmp_path / "absent") == []


# --------------------------------------------------------------------------
# 13. UI / CLI / provenance agree on canonical identity
# --------------------------------------------------------------------------


def test_provenance_and_timeline_agree_on_canonical_identity(historical_run):
    from pdbclean.stage_registry import producer_canonical_label

    recorded = {
        stage["stage_id"]: stage["canonical_stage"]
        for stage in historical_run.record["stages"]
    }

    for row in run_timeline(historical_run.record):
        if row["producer"] is None:
            continue

        assert row["producer"] in recorded
        assert recorded[row["producer"]] == producer_canonical_label(
            row["producer"]
        )


def test_registry_label_matches_the_canonical_timeline():
    """Provenance records exactly the identity the timeline displays."""

    from pdbclean.stage_registry import STAGES, producer_canonical_label

    for stage in STAGES:
        assert stage.canonical_stage == producer_canonical_label(
            stage.stage_id
        ), stage.stage_id


def test_every_canonical_identity_has_a_display_label():
    for entry in canonical_timeline():
        assert entry.display.startswith(entry.label)
        assert entry.title in entry.display
        assert "—" in entry.display


# --------------------------------------------------------------------------
# 14-15. Duplicate navigation; Mol* stays inspection-only
# --------------------------------------------------------------------------


def test_duplicate_stages_offer_explorer_navigation():
    for key in ("stage_8", "stage_9", "stage_10", "stage_14a", "stage_14b"):
        navigation = duplicate_navigation(key)

        assert navigation is not None
        assert navigation["available"] is True
        assert "Duplicate Explorer" in navigation["label"]


def test_non_duplicate_stages_offer_no_navigation():
    for key in ("prerequisite_a", "stage_1", "stage_5", "stage_14c"):
        assert duplicate_navigation(key) is None


def test_navigation_states_that_molstar_never_classifies():
    navigation = duplicate_navigation("stage_10")

    assert "never determines" in navigation["caveat"]
    assert "inspection only" in navigation["caveat"]


def test_removal_stages_filter_to_removed_chains():
    assert duplicate_navigation("stage_14b")["filters"] == {
        "relationship": "removed"
    }


# --------------------------------------------------------------------------
# Detailed scientific descriptions (section 31)
# --------------------------------------------------------------------------


def _by_key():
    return {entry.key: entry for entry in canonical_timeline()}


def test_every_canonical_stage_has_a_detailed_description():
    for entry in canonical_timeline():
        assert entry.rationale, entry.key
        assert entry.scientific_method, entry.key
        assert entry.stage_input, entry.key
        assert entry.stage_output, entry.key
        assert entry.downstream_role, entry.key
        assert entry.implementation_note, entry.key


def test_descriptions_are_substantial_enough_to_explain_the_method():
    for entry in canonical_timeline():
        if entry.role == "input":
            continue

        assert len(entry.scientific_method) > 150, entry.key


def test_every_prerequisite_is_described():
    entries = _by_key()

    for key in ("prerequisite_a", "prerequisite_b", "prerequisite_c"):
        assert entries[key].scientific_method
        assert entries[key].label.startswith("Prerequisite ")


def test_stage_three_and_four_have_distinct_descriptions():
    entries = _by_key()

    three, four = entries["stage_3"], entries["stage_4"]

    assert three.scientific_method != four.scientific_method
    assert "m x 9" in three.scientific_method
    assert "BRI_units = round(BRI / p)" in four.scientific_method


def test_stage_eight_and_nine_have_distinct_descriptions():
    entries = _by_key()

    eight, nine = entries["stage_8"], entries["stage_9"]

    assert eight.scientific_method != nine.scientific_method
    assert "cover tree" in eight.scientific_method
    assert "representation units" in nine.scientific_method


def test_ckdtree_is_described_only_for_brain_filtering():
    mentions = [
        entry.key
        for entry in canonical_timeline()
        if "cKDTree" in (entry.scientific_method + entry.implementation_note)
    ]

    assert mentions == ["stage_7"]

    stage_7 = _by_key()["stage_7"]

    assert "never the final complete-BRI search engine" in (
        stage_7.implementation_note
    )


def test_complete_bri_is_described_as_the_final_classification_basis():
    entries = _by_key()

    assert "complete BRI" in entries["stage_8"].rationale
    assert "complete-BRI L-infinity" in entries["stage_10"].scientific_method
    assert "Brain distance plays no part" in (
        entries["stage_10"].scientific_method
    )


def test_brain_is_described_as_filtering_only():
    stage_5 = _by_key()["stage_5"]

    assert "filtering and indexing layer only" in stage_5.scientific_method
    assert "never classifies duplicates" in stage_5.scientific_method


def test_components_are_described_as_not_equivalence_classes():
    stage_14a = _by_key()["stage_14a"]

    assert "NOT a duplicate equivalence class" in (
        stage_14a.scientific_method
    )
    assert "not transitive" in stage_14a.scientific_method


def test_direct_edge_requirement_is_described():
    stage_14b = _by_key()["stage_14b"]

    assert "DIRECT-EDGE" in stage_14b.scientific_method
    assert "NO transitive removal" in stage_14b.scientific_method
    assert "m = 1 are retained" in stage_14b.scientific_method


def test_investigation_stages_are_described_as_not_deletion_relations():
    entries = _by_key()

    for key in ("stage_11", "stage_12", "stage_13"):
        text = entries[key].downstream_role + entries[key].note

        assert "NOT" in text or "not a deletion relation" in text.lower()


def test_inclusive_threshold_is_described():
    stage_10 = _by_key()["stage_10"]

    assert "inclusive" in stage_10.scientific_method
    assert "d <= tau" in stage_10.scientific_method


def test_pair_counts_are_distinguished_from_chain_counts():
    stage_10 = _by_key()["stage_10"]

    assert "not the number of chains" in stage_10.scientific_method


def test_precision_is_distinguished_from_the_threshold():
    stage_4 = _by_key()["stage_4"]

    assert "NOT a duplicate threshold" in stage_4.scientific_method


def test_comp702_choices_are_not_attributed_to_papers():
    """A computational choice must be labelled as one."""

    entries = _by_key()

    assert "not a quantisation step prescribed by the paper" in (
        entries["stage_4"].implementation_note
    )
    assert "COMP702 engineering implementation" in (
        entries["stage_7"].implementation_note
    )
    assert "COMP702" in entries["stage_14b"].implementation_note


def test_references_resolve_to_real_records():
    from pdbclean.stage_registry import METHOD_REFERENCES, method_references

    for entry in canonical_timeline():
        for key in entry.references:
            assert key in METHOD_REFERENCES, (entry.key, key)

        for reference in method_references(entry.references):
            assert reference["doi"]
            assert reference["authors"]
            assert reference["year"]


def test_bri_stages_cite_the_match_paper():
    entries = _by_key()

    for key in ("stage_3", "stage_5"):
        assert "anosova_match_2025" in entries[key].references


def test_duplicate_stages_cite_the_acta_paper():
    entries = _by_key()

    for key in ("stage_8", "stage_10", "stage_13"):
        assert "wlodawer_acta_2025" in entries[key].references


def test_stage_detail_exposes_the_description(historical_run):
    detail = stage_detail(historical_run.record, "stage_7", repo_root=REPO)

    description = detail["description"]

    assert "cKDTree" in description["implementation_note"]
    assert description["rationale"]
    assert description["references"][0]["doi"] == "10.46793/match.94-1.097A"
