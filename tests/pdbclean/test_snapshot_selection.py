"""Snapshot-selection regression tests.

The pipeline must default to the latest complete snapshot while keeping every
historical snapshot -- in particular the frozen 2026-01-01 production snapshot
-- reproducible by explicit selection.  A run is always pinned to a concrete
snapshot identity before any processing begins.
"""

from __future__ import annotations

import pytest

from pdbclean.runconfig import resolve_run_config
from pdbclean.snapshot import ResolvedSnapshot
from pdbclean.snapshot_selection import (
    SnapshotSelectionError,
    format_snapshot_id,
    interpret_menu_response,
    interpret_selection,
    list_available_snapshots,
    normalise_snapshot_id,
    render_snapshot_menu,
    resolve_snapshot_for_run,
    snapshot_provenance,
)


DISCOVERED = ["20250101", "20260101", "20260401", "not-a-date", "20260701"]


def _discover(**_kwargs):
    return list(DISCOVERED)


def _choices():
    return list_available_snapshots(
        bucket_url="https://example.invalid",
        discover=_discover,
    )


def _fake_resolver(snapshot_config, **_kwargs):
    snapshot_id = snapshot_config.get("snapshot_id") or "20260701"

    return ResolvedSnapshot(
        snapshot_id=snapshot_id,
        layout="pub/pdb/data/structures/divided/mmCIF",
        source_prefix=f"{snapshot_id}/pub/pdb/data/structures/divided/mmCIF",
        sample_mmcif_key=f"{snapshot_id}/.../1abc.cif.gz",
        selection_mode=snapshot_config.get("mode", "fixed"),
    )


# --------------------------------------------------------------------------
# Identity normalisation
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value", ["2026-01-01", "20260101", " 2026-01-01 ", "20260101\n"]
)
def test_dates_normalise_to_the_archive_identity(value):
    assert normalise_snapshot_id(value) == "20260101"


@pytest.mark.parametrize("value", ["2026-1-1", "26-01-01", "January", "202601"])
def test_malformed_dates_are_rejected(value):
    with pytest.raises(SnapshotSelectionError):
        normalise_snapshot_id(value)


def test_display_form_is_the_iso_date():
    assert format_snapshot_id("20260101") == "2026-01-01"


# --------------------------------------------------------------------------
# Discovery and the picker
# --------------------------------------------------------------------------


def test_discovery_is_newest_first_and_drops_non_dates():
    choices = _choices()

    assert [c.snapshot_id for c in choices] == [
        "20260701",
        "20260401",
        "20260101",
        "20250101",
    ]
    assert choices[0].is_latest is True
    assert all(c.is_latest is False for c in choices[1:])


def test_discovery_limit_keeps_the_newest():
    choices = list_available_snapshots(
        bucket_url="https://example.invalid",
        limit=2,
        discover=_discover,
    )

    assert [c.snapshot_id for c in choices] == ["20260701", "20260401"]


def test_menu_offers_latest_complete_as_position_one():
    menu = render_snapshot_menu(_choices())

    assert "1. latest complete snapshot [default]" in menu
    assert "2. 2026-07-01" in menu
    assert "4. 2026-01-01" in menu


@pytest.mark.parametrize("response", [None, "", "   ", "latest", "default"])
def test_empty_response_means_latest_complete(response):
    assert interpret_menu_response(response, _choices()) is None


def test_menu_position_one_means_latest_complete():
    assert interpret_menu_response("1", _choices()) is None


def test_menu_position_selects_a_historical_snapshot():
    assert interpret_menu_response("4", _choices()) == "20260101"


def test_menu_accepts_an_explicit_date():
    assert interpret_menu_response("2026-01-01", _choices()) == "20260101"


def test_menu_rejects_an_out_of_range_position():
    with pytest.raises(SnapshotSelectionError):
        interpret_menu_response("99", _choices())


def test_list_position_selection_is_one_based_on_the_bare_list():
    assert interpret_selection("3", _choices()) == "20260101"
    assert interpret_selection(None, _choices()) is None


# --------------------------------------------------------------------------
# Resolution for a run
# --------------------------------------------------------------------------


def test_default_run_resolves_to_the_latest_complete_snapshot():
    resolved, snapshot = resolve_snapshot_for_run(
        resolve_run_config(),
        resolver=_fake_resolver,
    )

    assert snapshot.snapshot_id == "20260701"
    assert resolved.get("snapshot.snapshot_id") == "20260701"
    assert resolved.get("snapshot.resolved_selection_mode") == "latest_complete"


def test_frozen_snapshot_is_reproducible_by_explicit_selection():
    resolved, snapshot = resolve_snapshot_for_run(
        resolve_run_config(
            overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
        ),
        resolver=_fake_resolver,
    )

    assert snapshot.snapshot_id == "20260101"
    assert resolved.get("snapshot.snapshot_id") == "20260101"


def test_explicit_request_beats_the_configuration_file():
    resolved, snapshot = resolve_snapshot_for_run(
        resolve_run_config(
            overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
        ),
        requested="2026-04-01",
        resolver=_fake_resolver,
    )

    assert snapshot.snapshot_id == "20260401"
    assert resolved.get("snapshot.snapshot_id") == "20260401"


def test_fixed_mode_without_an_identity_is_refused():
    """Never silently fall back to "latest" when a fixed run was requested."""

    unpinned = resolve_run_config(overrides=["snapshot.mode=fixed"], validate=False)

    assert unpinned.get("snapshot.snapshot_id") is None

    with pytest.raises(SnapshotSelectionError):
        resolve_snapshot_for_run(unpinned, resolver=_fake_resolver)


def test_run_is_pinned_before_any_processing():
    """A resolved run must be projectable onto the stage schema at once."""

    resolved, _ = resolve_snapshot_for_run(
        resolve_run_config(),
        resolver=_fake_resolver,
    )

    projected = resolved.to_protocol_config()

    assert projected["snapshot"]["mode"] == "fixed"
    assert projected["snapshot"]["snapshot_id"] == "20260701"


def test_provenance_records_the_full_snapshot_identity():
    _, snapshot = resolve_snapshot_for_run(
        resolve_run_config(overrides=["snapshot.snapshot_id=20260101"]),
        resolver=_fake_resolver,
    )

    payload = snapshot_provenance(snapshot)

    assert payload["snapshot_id"] == "20260101"
    assert payload["display"] == "2026-01-01"
    assert payload["layout"]
    assert payload["source_prefix"]
    assert payload["selection_mode"]
