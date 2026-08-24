"""The configuration a worker executes must be the one the operator froze.

These tests exist because of a specific, silent failure: every generated stage
command used to point at ``config/pdbclean/protocol_3_2_comp702_v1.yaml``.
That file has no ``brain_filter`` section, so ``brain_prefilter_production``
fell back to its module constant and computed with 0.010 A no matter what the
operator entered.  It also pins ``duplicate_search`` to 0.010 and
``snapshot.mode`` to ``latest_complete``, so an entered tau and an entered
snapshot were discarded too.  Nothing raised.

Everything below is written to make that class of failure loud.
"""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import pytest
import yaml

from pdbclean.config import load_config
from pdbclean.runconfig import resolve_run_config
from pdbclean.stage_config import (
    STAGE_CONFIG_BASENAME,
    StageConfigError,
    cached_stage_config,
    compare_values,
    executed_values,
    frozen_values,
    stage_config_document,
    stage_config_path,
    stage_config_text,
    verify_stage_config,
    write_stage_config,
)


REPO_ROOT = Path(__file__).resolve().parents[2]

#: The sections the production stages read.  A projection that omits any of
#: them lets a stage fall back to a source-code default.
REQUIRED_SECTIONS = (
    "release",
    "snapshot",
    "selection",
    "quality_rules",
    "post_cleaning_geometric_validation",
    "bri",
    "brain",
    "brain_filter",
    "duplicate_search",
    "graph",
    "representative_selection",
    "execution",
    "storage",
)


def _resolved(**overrides):
    items = {
        "snapshot.mode": "fixed",
        "snapshot.snapshot_id": "20260101",
        **overrides,
    }

    return resolve_run_config(
        overrides=[f"{key}={value}" for key, value in items.items()]
    )


def _tau(value):
    """A scientifically valid threshold pair.

    Brain is a lossless prefilter for the complete-BRI search, so its radius
    may never be smaller than tau. Lowering tau therefore means lowering both
    -- which is exactly what the tau-sensitivity experiment does, and what
    `validate_resolved_config` enforces.
    """

    return _resolved(
        **{
            "brain_filter.threshold_angstrom": value,
            "duplicate_search.near_duplicate_threshold_angstrom": value,
        }
    )


# --------------------------------------------------------------------------
# Completeness
# --------------------------------------------------------------------------


@pytest.mark.parametrize("section", REQUIRED_SECTIONS)
def test_projection_carries_every_section_a_stage_reads(section):
    document = stage_config_document(_resolved())

    assert section in document, (
        f"{section} is missing from the executable configuration; a stage "
        "would fall back to its own default"
    )


def test_projection_carries_the_brain_threshold():
    """The exact omission that silently replaced an entered 0.005 with 0.010."""

    document = stage_config_document(_tau(0.005))

    assert document["brain_filter"]["threshold_angstrom"] == 0.005
    assert document["duplicate_search"][
        "near_duplicate_threshold_angstrom"
    ] == 0.005


def test_projection_carries_the_representation_precision():
    document = stage_config_document(_resolved())

    assert document["bri"]["representation_precision_angstrom"] == 0.001


def test_projection_pins_the_snapshot_as_fixed():
    """The base YAML says latest_complete; a run must never inherit that."""

    document = stage_config_document(_resolved())

    assert document["snapshot"]["mode"] == "fixed"
    assert document["snapshot"]["snapshot_id"] == "20260101"


def test_projection_embeds_the_canonical_identity():
    resolved = _resolved()
    document = stage_config_document(resolved)

    assert document["resolved_run"]["resolved_config_sha256"] == resolved.sha256
    assert document["resolved_run"]["scientific_config_sha256"] == (
        resolved.scientific_sha256
    )


# --------------------------------------------------------------------------
# The legacy compatibility section
# --------------------------------------------------------------------------


def test_geometric_search_is_present_but_marked_superseded():
    document = stage_config_document(_resolved())

    assert document["geometric_search"]["status"] == "superseded_not_executed"
    assert "1.0" not in yaml.safe_dump(document["geometric_search"])


def test_no_production_module_reads_geometric_search():
    """The compatibility block is inert, and this is what proves it.

    ``pdbclean.config.load_config`` lists ``geometric_search`` as required, so
    the projection has to carry it.  It describes the superseded PDB707K
    design.  If any stage ever started reading it, carrying an inert marker
    would become a silent scientific change -- so that is a test failure.
    """

    offenders = []

    for path in sorted(
        list((REPO_ROOT / "src" / "pdbclean").rglob("*.py"))
        + list((REPO_ROOT / "scripts").rglob("*.py"))
    ):
        source = path.read_text(encoding="utf-8")

        if "geometric_search" not in source:
            continue

        tree = ast.parse(source)

        for node in ast.walk(tree):
            reads = (
                isinstance(node, ast.Subscript)
                and isinstance(node.slice, ast.Constant)
                and node.slice.value == "geometric_search"
            ) or (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "get"
                and node.args
                and isinstance(node.args[0], ast.Constant)
                and node.args[0].value == "geometric_search"
            )

            if reads:
                offenders.append(f"{path.relative_to(REPO_ROOT)}:{node.lineno}")

    assert not offenders, (
        "geometric_search is now read by production code: " + ", ".join(offenders)
    )


# --------------------------------------------------------------------------
# Determinism
# --------------------------------------------------------------------------


def test_projection_is_deterministic():
    first = stage_config_text(_resolved())
    second = stage_config_text(_resolved())

    assert first == second


def test_projection_does_not_depend_on_the_writing_host(monkeypatch):
    monkeypatch.setenv("TMPDIR", "/scratch/login-node")
    first = stage_config_text(_resolved())

    monkeypatch.setenv("TMPDIR", "/local/compute-node/7788")
    second = stage_config_text(_resolved())

    assert first == second


def test_different_thresholds_give_different_documents():
    assert stage_config_text(_tau(0.005)) != stage_config_text(_resolved())


# --------------------------------------------------------------------------
# The production loader
# --------------------------------------------------------------------------


def test_projection_loads_through_the_production_loader(tmp_path):
    written = write_stage_config(_resolved(), tmp_path)

    loaded = load_config(written["stage_config_path"])

    assert loaded.data["release"]["protocol_version"] == (
        "protocol3.2-comp702-v1"
    )


def test_worker_loads_the_entered_thresholds(tmp_path):
    """The whole point, measured at the worker boundary."""

    resolved = _resolved(
        **{
            "bri.representation_precision_angstrom": 0.001,
            "brain_filter.threshold_angstrom": 0.005,
            "duplicate_search.near_duplicate_threshold_angstrom": 0.005,
        }
    )

    written = write_stage_config(resolved, tmp_path)

    values = executed_values(load_config(written["stage_config_path"]).data)

    assert values["representation_precision_angstrom"] == 0.001
    assert values["brain_filter_threshold_angstrom"] == 0.005
    assert values["brain_filter_threshold_units"] == 5
    assert values["complete_bri_near_duplicate_threshold_angstrom"] == 0.005
    assert values["complete_bri_near_duplicate_threshold_units"] == 5


def test_the_frozen_base_yaml_would_have_executed_the_wrong_values():
    """Document the defect this architecture removes.

    Driving a stage from the byte-frozen base configuration -- which is what
    every generated command used to do -- makes the Brain worker fall back to
    its module constant. This test pins that behaviour so the regression is
    visible rather than theoretical.
    """

    base = load_config(
        REPO_ROOT / "config/pdbclean/protocol_3_2_comp702_v1.yaml"
    ).data

    assert "brain_filter" not in base

    values = executed_values(base)

    # The operator's 0.005 never reaches this; 0.010 comes from the module.
    assert values["brain_filter_threshold_angstrom"] == 0.010
    assert values["complete_bri_near_duplicate_threshold_angstrom"] == 0.010
    assert values["snapshot"] is None


# --------------------------------------------------------------------------
# Writing, verifying, tampering
# --------------------------------------------------------------------------


def test_write_records_the_digest_of_what_it_wrote(tmp_path):
    resolved = _resolved()
    written = write_stage_config(resolved, tmp_path)

    on_disk = Path(written["stage_config_path"]).read_bytes()

    assert written["stage_config_sha256"] == hashlib.sha256(on_disk).hexdigest()
    assert Path(written["stage_config_path"]).name == STAGE_CONFIG_BASENAME


def test_writing_the_same_configuration_twice_is_idempotent(tmp_path):
    resolved = _resolved()

    first = write_stage_config(resolved, tmp_path)
    second = write_stage_config(resolved, tmp_path)

    assert first == second


def test_refuses_to_overwrite_a_frozen_configuration(tmp_path):
    write_stage_config(_resolved(), tmp_path)

    with pytest.raises(StageConfigError, match="different content"):
        write_stage_config(_tau(0.005), tmp_path)


def test_verification_detects_an_edited_configuration(tmp_path):
    resolved = _resolved()
    write_stage_config(resolved, tmp_path)

    path = stage_config_path(tmp_path)
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    document["brain_filter"]["threshold_angstrom"] = 0.010
    path.write_text(yaml.safe_dump(document), encoding="utf-8")

    with pytest.raises(StageConfigError, match="not the projection"):
        verify_stage_config(tmp_path, resolved)


def test_verification_reports_a_missing_configuration(tmp_path):
    with pytest.raises(StageConfigError, match="no executable stage"):
        verify_stage_config(tmp_path, _resolved())


def test_cached_projection_is_content_addressed(tmp_path):
    resolved = _resolved()

    first = cached_stage_config(resolved, tmp_path)
    second = cached_stage_config(resolved, tmp_path)

    assert first == second
    assert first.stem == resolved.sha256

    other = cached_stage_config(_tau(0.005), tmp_path)

    assert other != first


# --------------------------------------------------------------------------
# Comparing executed values with frozen values
# --------------------------------------------------------------------------


def test_matching_values_produce_no_mismatches(tmp_path):
    resolved = _tau(0.005)

    written = write_stage_config(resolved, tmp_path)

    executed = executed_values(load_config(written["stage_config_path"]).data)

    assert compare_values(executed, frozen_values(resolved)) == []


def test_a_lost_brain_override_is_detected(tmp_path):
    """A projection that stopped carrying brain_filter must fail, not compute."""

    resolved = _tau(0.005)

    written = write_stage_config(resolved, tmp_path)

    path = Path(written["stage_config_path"])
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    del document["brain_filter"]
    path.write_text(yaml.safe_dump(document), encoding="utf-8")

    executed = executed_values(load_config(path).data)
    mismatches = compare_values(executed, frozen_values(resolved))

    keys = {item["key"] for item in mismatches}

    assert "brain_filter_threshold_angstrom" in keys

    lost = next(
        item
        for item in mismatches
        if item["key"] == "brain_filter_threshold_angstrom"
    )

    assert lost["executed"] == 0.010
    assert lost["frozen"] == 0.005
