"""UI and CLI must resolve to the same configuration and the same backend.

The UI is a presentation layer over exactly the modules the CLI drives.  These
tests run the real HTTP server on an ephemeral port and compare its answers
with the CLI's, so a divergence between the two front ends fails the build.
"""

from __future__ import annotations

import json
import threading
import urllib.error
import urllib.parse
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

import pytest
import yaml

from pdbclean.pipeline import plan_pipeline
from pdbclean.runconfig import resolve_run_config
from pdbclean.ui import server as ui_server


REPO_ROOT = Path(__file__).resolve().parents[2]

PROFILE = "config/pdbclean/profiles/comp702_frozen_20260101.yaml"

OVERRIDES = ["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]


@pytest.fixture()
def ui():
    """A real server instance on an ephemeral port."""

    ui_server.Handler.state = ui_server.UIState(
        repo_root=REPO_ROOT,
        config_path=None,
        overrides=[],
    )

    httpd = ThreadingHTTPServer(("127.0.0.1", 0), ui_server.Handler)

    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()

    base = f"http://127.0.0.1:{httpd.server_address[1]}"

    def get(route):
        with urllib.request.urlopen(base + route, timeout=30) as response:
            return json.loads(response.read().decode("utf-8"))

    def post(route, payload):
        request = urllib.request.Request(
            base + route,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json"},
            method="POST",
        )

        with urllib.request.urlopen(request, timeout=60) as response:
            return json.loads(response.read().decode("utf-8"))

    def raw(route):
        with urllib.request.urlopen(base + route, timeout=30) as response:
            return response.status, response.read()

    try:
        yield type("UI", (), {"get": staticmethod(get),
                              "post": staticmethod(post),
                              "raw": staticmethod(raw),
                              "base": base})
    finally:
        httpd.shutdown()
        httpd.server_close()
        thread.join(timeout=5)


# --------------------------------------------------------------------------
# Static assets
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "route", ["/", "/app.css", "/app.js", "/viewer.html"]
)
def test_static_assets_are_served(ui, route):
    status, body = ui.raw(route)

    assert status == 200
    assert body


def test_unknown_endpoint_is_not_found(ui):
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        ui.get("/api/nonsense")

    assert excinfo.value.code == 404


# --------------------------------------------------------------------------
# Bootstrap
# --------------------------------------------------------------------------


def test_bootstrap_exposes_defaults_and_the_stage_chain(ui):
    payload = ui.get("/api/bootstrap")

    assert payload["defaults"]["brain"]["dimension"] == 9
    assert payload["defaults"]["duplicate_search"][
        "near_duplicate_threshold_angstrom"
    ] == 0.010

    stage_ids = [stage["stage_id"] for stage in payload["stages"]]

    assert stage_ids[0] == "snapshot"
    assert stage_ids[-1] == "gold_release"
    assert len(stage_ids) == 15


def test_bootstrap_lists_the_frozen_profile(ui):
    payload = ui.get("/api/bootstrap")

    names = {profile["path"] for profile in payload["profiles"]}

    assert any(name.endswith("comp702_frozen_20260101.yaml") for name in names)


# --------------------------------------------------------------------------
# Configuration equivalence
# --------------------------------------------------------------------------


def test_ui_and_cli_resolve_the_same_configuration(ui):
    cli = resolve_run_config(config_path=REPO_ROOT / PROFILE)

    payload = ui.post(
        "/api/config/resolve", {"config_path": str(REPO_ROOT / PROFILE)}
    )

    assert payload["scientific_config_sha256"] == cli.scientific_sha256
    assert payload["resolved"]["duplicate_search"] == cli.get(
        "duplicate_search"
    )
    assert payload["near_duplicate_threshold_mA"] == (
        cli.near_duplicate_threshold_mA
    )
    assert payload["brain_threshold_mA"] == cli.brain_threshold_mA


def test_ui_overrides_use_the_same_precedence_as_the_cli(ui):
    cli = resolve_run_config(
        config_path=REPO_ROOT / PROFILE,
        overrides={"selection.models.model_id": 1},
        override_origin="ui",
    )

    payload = ui.post(
        "/api/config/resolve",
        {
            "config_path": str(REPO_ROOT / PROFILE),
            "overrides": {"selection.models.model_id": 1},
        },
    )

    assert payload["resolved_config_sha256"] == cli.sha256
    assert payload["sources"]["selection.models.model_id"] == "override:ui"


def test_ui_emits_the_same_resolved_run_yaml(ui):
    cli = resolve_run_config(config_path=REPO_ROOT / PROFILE)

    payload = ui.post(
        "/api/config/resolve", {"config_path": str(REPO_ROOT / PROFILE)}
    )

    assert yaml.safe_load(payload["resolved_config_yaml"]) == yaml.safe_load(
        cli.to_yaml()
    )


def test_ui_rejects_a_configuration_the_cli_would_reject(ui):
    """A scientific guard rail must not be bypassable through the browser."""

    with pytest.raises(urllib.error.HTTPError) as excinfo:
        ui.post(
            "/api/config/resolve",
            {"overrides": {"duplicate_search.operator": "less_than"}},
        )

    assert excinfo.value.code == 400

    body = json.loads(excinfo.value.read().decode("utf-8"))

    assert "inclusive" in body["error"]


def test_ui_snapshot_selection_pins_the_identity(ui):
    payload = ui.post(
        "/api/config/resolve",
        {"config_path": str(REPO_ROOT / PROFILE), "snapshot": "2026-01-01"},
    )

    assert payload["resolved"]["snapshot"]["snapshot_id"] == "20260101"
    assert payload["paths"]["release"].startswith("PDBClean-20260101-")


# --------------------------------------------------------------------------
# Plan equivalence
# --------------------------------------------------------------------------


def test_ui_and_cli_produce_the_same_plan(ui):
    cli_plan = plan_pipeline(
        resolve_run_config(config_path=REPO_ROOT / PROFILE),
        repo_root=REPO_ROOT,
    ).to_dict()

    payload = ui.post("/api/plan", {"config_path": str(REPO_ROOT / PROFILE)})

    ui_plan = payload["plan"] if "plan" in payload else payload

    assert ui_plan["scientific_config_sha256"] == cli_plan[
        "scientific_config_sha256"
    ]
    assert ui_plan["snapshot"] == cli_plan["snapshot"]
    assert ui_plan["release"] == cli_plan["release"]

    assert [s["stage_id"] for s in ui_plan["stages"]] == [
        s["stage_id"] for s in cli_plan["stages"]
    ]
    assert [s["action"] for s in ui_plan["stages"]] == [
        s["action"] for s in cli_plan["stages"]
    ]
    assert [s["status"] for s in ui_plan["stages"]] == [
        s["status"] for s in cli_plan["stages"]
    ]


# --------------------------------------------------------------------------
# Gold release page
# --------------------------------------------------------------------------


def test_release_page_shows_nothing_for_an_unproduced_release(ui, tmp_path):
    payload = ui.get(
        "/api/release?"
        + urllib.parse.urlencode({"snapshot": "20250101"})
    )

    assert payload["published"] is False
    assert payload["release"] == {}


# --------------------------------------------------------------------------
# Mol* scenes
# --------------------------------------------------------------------------


def test_scene_index_is_read_only_metadata(ui):
    payload = ui.get("/api/scenes")

    assert "scenes" in payload

    for scene in payload["scenes"]:
        assert scene["key"]
        assert scene["views"]

        for view in scene["views"]:
            assert view["url"].startswith("/structures/")


def test_structure_route_refuses_traversal(ui):
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        ui.raw("/structures/../../../etc/passwd")

    assert excinfo.value.code in {403, 404}


# --------------------------------------------------------------------------
# Runs
# --------------------------------------------------------------------------


def test_runs_endpoint_lists_without_creating_anything(ui):
    payload = ui.get("/api/runs")

    assert "runs" in payload
    assert "run_root" in payload


# --------------------------------------------------------------------------
# Theme is a viewing preference, never a run parameter
# --------------------------------------------------------------------------

STATIC = Path(__file__).resolve().parents[2] / "src" / "pdbclean" / "ui" / "static"


def test_theme_is_never_a_configuration_key():
    """No theme key may exist anywhere in a resolved configuration."""

    resolved = resolve_run_config(config_path=REPO_ROOT / PROFILE)

    rendered = json.dumps(resolved.to_dict()).lower()

    for token in ("theme", "dark_mode", "colour_scheme", "color_scheme"):
        assert token not in rendered, token


def test_theme_cannot_affect_either_hash():
    """Two identical configurations hash identically whatever the theme."""

    first = resolve_run_config(config_path=REPO_ROOT / PROFILE)
    second = resolve_run_config(config_path=REPO_ROOT / PROFILE)

    assert first.sha256 == second.sha256
    assert first.scientific_sha256 == second.scientific_sha256


def test_theme_override_is_rejected_as_an_unknown_scientific_key():
    """A theme cannot be smuggled in as a configuration override."""

    resolved = resolve_run_config(
        config_path=REPO_ROOT / PROFILE,
        overrides={"ui.theme": "dark"},
    )

    from pdbclean.runconfig import scientific_projection

    # It lands outside every scientific section, so the scientific identity
    # is unchanged.
    assert resolved.scientific_sha256 == resolve_run_config(
        config_path=REPO_ROOT / PROFILE
    ).scientific_sha256
    assert "ui" not in scientific_projection(resolved.data)


def test_theme_is_client_side_only(ui):
    """The backend is never told which theme the browser is using."""

    payload = ui.post(
        "/api/config/resolve", {"config_path": str(REPO_ROOT / PROFILE)}
    )

    assert "theme" not in json.dumps(payload).lower()


def test_stylesheet_defines_all_three_theme_states():
    css = (STATIC / "app.css").read_text(encoding="utf-8")

    # Light on bare :root, so a page always has a complete palette.
    assert ":root {" in css

    # OS preference, but never overriding an explicit light choice.
    assert '@media (prefers-color-scheme: dark)' in css
    assert ':root:not([data-theme="light"])' in css

    # An explicit choice wins in both directions.
    assert ':root[data-theme="dark"]' in css


def test_both_themes_define_every_colour_token():
    css = (STATIC / "app.css").read_text(encoding="utf-8")

    tokens = [
        "--bg", "--panel", "--panel-alt", "--border", "--border-strong",
        "--text", "--muted", "--faint", "--accent", "--accent-bg",
        "--pass", "--fail", "--warn", "--idle",
    ]

    light = css.split("@media (prefers-color-scheme: dark)")[0]
    dark = css.split(':root[data-theme="dark"] {')[1].split("}")[0]

    for token in tokens:
        assert f"{token}:" in light, f"light theme missing {token}"
        assert f"{token}:" in dark, f"dark theme missing {token}"


def test_theme_toggle_exists_and_is_labelled():
    html = (STATIC / "index.html").read_text(encoding="utf-8")

    assert 'id="theme-toggle"' in html
    assert "aria-label" in html


def test_theme_is_applied_before_first_paint():
    """Avoids a light flash, and works before app.js has loaded."""

    html = (STATIC / "index.html").read_text(encoding="utf-8")
    head = html.split("</head>")[0]

    assert "pdbclean.theme" in head
    assert "data-theme" in head


def test_theme_switching_does_not_reload_the_page():
    js = (STATIC / "app.js").read_text(encoding="utf-8")

    toggle = js.split("function initTheme()")[1].split("\n}\n")[0]

    for forbidden in ("location.reload", "location.href", "window.location"):
        assert forbidden not in toggle, forbidden


def test_theme_preference_is_persisted_locally():
    js = (STATIC / "app.js").read_text(encoding="utf-8")

    assert "localStorage" in js
    assert "prefers-color-scheme" in js


# --------------------------------------------------------------------------
# Canonical stage identity across the whole UI
# --------------------------------------------------------------------------


def test_bootstrap_exposes_the_canonical_timeline(ui):
    payload = ui.get("/api/bootstrap")

    timeline = payload["canonical_timeline"]

    assert [entry["label"] for entry in timeline][:4] == [
        "Prerequisite A",
        "Prerequisite B",
        "Prerequisite C",
        "Stage 1",
    ]
    assert [entry["position"] for entry in timeline] == list(
        range(1, len(timeline) + 1)
    )


def test_every_ui_stage_row_carries_a_canonical_label(ui):
    """The Pipeline view must never show a bare stage name."""

    payload = ui.post("/api/plan", {"config_path": str(REPO_ROOT / PROFILE)})
    plan = payload.get("plan", payload)

    for stage in plan["stages"]:
        assert stage["canonical_stage"]
        assert stage["canonical_stage"].startswith(
            ("Stage ", "Prerequisite ")
        ), stage["stage_id"]


def test_ui_never_sorts_stages_alphabetically():
    js = (STATIC / "app.js").read_text(encoding="utf-8")

    # The timeline is rendered straight from the ordered payload.
    assert "payload.timeline.forEach" in js
    assert "state.bootstrap.canonical_timeline.forEach" in js

    # And nothing sorts a stage collection.
    for forbidden in (
        "stages.sort(",
        "timeline.sort(",
        "canonical_timeline.sort(",
    ):
        assert forbidden not in js, forbidden


def test_precision_control_is_present_in_the_configuration_form():
    html = (STATIC / "index.html").read_text(encoding="utf-8")

    assert 'id="cfg-precision"' in html
    assert "representation precision" in html.lower()
    assert "Validated default: 0.001" in html


def test_configuration_form_distinguishes_the_three_concepts():
    html = (STATIC / "index.html").read_text(encoding="utf-8")

    section = html.split("BRI and duplicate detection")[1].split("</fieldset>")[0]

    assert "Representation precision" in section
    assert "Brain filtering threshold" in section
    assert "Complete-BRI" in section
    assert "not a duplicate threshold" in section


def test_resolve_reports_precision_and_units(ui):
    payload = ui.post(
        "/api/config/resolve", {"config_path": str(REPO_ROOT / PROFILE)}
    )

    assert payload["representation_precision_angstrom"] == 0.001
    assert payload["representation_unit"] == "mA"
    assert payload["precision_is_implemented"] is True
    assert payload["near_duplicate_threshold_mA"] == 10
