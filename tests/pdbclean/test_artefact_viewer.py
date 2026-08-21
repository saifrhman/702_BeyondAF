"""Artefact Viewer: bounded inspection, faithful download, strict paths.

The viewer is a provenance tool, not a filesystem browser. These tests pin
that it shows the real records, never materialises a whole table, hands back
the original bytes unmodified, and refuses everything outside its allowlist.
"""

from __future__ import annotations

import hashlib
import json
import threading
import urllib.error
import urllib.parse
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from pdbclean.artefacts import (
    DEFAULT_PAGE_SIZE,
    MAX_PAGE_SIZE,
    ArtefactError,
    csv_page,
    parquet_page,
    parquet_schema,
    rows_to_csv,
    table_page,
)
from pdbclean.ui import server as ui_server


REPO = Path(__file__).resolve().parents[2]
PROFILE = "config/pdbclean/profiles/comp702_frozen_20260101.yaml"

RELEASE = (
    REPO
    / "outputs"
    / "releases"
    / "PDBClean-20260101-protocol3.2-comp702-v1-dedup-v1"
)


def _frozen(path: Path) -> Path:
    if not path.is_file():
        pytest.skip(f"frozen artefact not present: {path}")

    return path


@pytest.fixture()
def table(tmp_path):
    """A deterministic Parquet table with known contents."""

    target = tmp_path / "rows.parquet"

    pq.write_table(
        pa.table(
            {
                "n": list(range(1000)),
                "name": [f"chain-{index:04d}" for index in range(1000)],
                "score": [float(1000 - index) for index in range(1000)],
            }
        ),
        target,
    )

    return target


@pytest.fixture()
def ui():
    ui_server.Handler.state = ui_server.UIState(
        repo_root=REPO, config_path=PROFILE, overrides=[]
    )

    httpd = ThreadingHTTPServer(("127.0.0.1", 0), ui_server.Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()

    base = f"http://127.0.0.1:{httpd.server_address[1]}"

    def get(route, **params):
        url = base + route

        if params:
            url += "?" + urllib.parse.urlencode(params)

        with urllib.request.urlopen(url, timeout=120) as response:
            return json.loads(response.read().decode("utf-8"))

    def raw(route, **params):
        url = base + route

        if params:
            url += "?" + urllib.parse.urlencode(params)

        with urllib.request.urlopen(url, timeout=120) as response:
            return response.status, response.read(), dict(response.headers)

    try:
        yield type(
            "UI",
            (),
            {"get": staticmethod(get), "raw": staticmethod(raw)},
        )
    finally:
        httpd.shutdown()
        httpd.server_close()
        thread.join(timeout=5)


# --------------------------------------------------------------------------
# Schema and row count
# --------------------------------------------------------------------------


def test_schema_is_returned_without_reading_rows(table):
    info = parquet_schema(table)

    assert info["row_count"] == 1000
    assert info["column_count"] == 3
    assert [column["name"] for column in info["columns"]] == [
        "n",
        "name",
        "score",
    ]
    assert "int64" in info["columns"][0]["type"]


def test_frozen_artefact_schema_matches_the_release():
    target = _frozen(RELEASE / "audit" / "removed_chain_audit.parquet")

    info = parquet_schema(target)

    assert info["row_count"] == 78_754
    assert any(c["name"] == "pdb_id" for c in info["columns"])
    assert any(c["name"] == "label_chain_id" for c in info["columns"])


# --------------------------------------------------------------------------
# Pagination
# --------------------------------------------------------------------------


def test_first_page_returns_the_first_rows(table):
    page = parquet_page(table, page=1, page_size=5)

    assert page["returned"] == 5
    assert [row["n"] for row in page["rows"]] == [0, 1, 2, 3, 4]
    assert page["row_count"] == 1000
    assert page["page_count"] == 200


def test_pages_are_disjoint_and_ordered(table):
    first = parquet_page(table, page=1, page_size=10)
    second = parquet_page(table, page=2, page_size=10)

    assert [row["n"] for row in first["rows"]] == list(range(10))
    assert [row["n"] for row in second["rows"]] == list(range(10, 20))


def test_rows_correspond_to_the_actual_file(table):
    page = parquet_page(table, page=3, page_size=4)

    for row in page["rows"]:
        assert row["name"] == f"chain-{row['n']:04d}"
        assert row["score"] == float(1000 - row["n"])


def test_page_size_is_capped(table):
    page = parquet_page(table, page=1, page_size=10_000)

    assert page["page_size"] == MAX_PAGE_SIZE
    assert page["returned"] <= MAX_PAGE_SIZE


def test_page_size_and_page_are_floored(table):
    page = parquet_page(table, page=0, page_size=0)

    assert page["page"] == 1
    assert page["page_size"] >= 1


def test_default_page_size_is_modest():
    assert DEFAULT_PAGE_SIZE <= 50
    assert MAX_PAGE_SIZE <= 200


def test_a_large_table_is_never_fully_materialised():
    """Page 1 of a multi-million-row table must not read the whole file."""

    target = _frozen(
        REPO
        / "outputs"
        / "pdbclean"
        / "20260101"
        / "protocol3.2-comp702-v1"
        / "duplicate_classification"
        / "finalized"
        / "candidate_classifications.parquet"
    )

    page = parquet_page(target, page=1, page_size=10)

    assert page["row_count"] > 1_000_000
    assert page["returned"] == 10
    assert page["scanned_rows"] < 100_000


# --------------------------------------------------------------------------
# Search and sort
# --------------------------------------------------------------------------


def test_search_filters_rows(table):
    page = parquet_page(table, page=1, page_size=10, search="chain-0042")

    assert page["matched_rows"] == 1
    assert page["rows"][0]["n"] == 42


def test_sort_orders_the_scanned_window(table):
    page = parquet_page(
        table, page=1, page_size=5, sort_by="score", descending=True
    )

    scores = [row["score"] for row in page["rows"]]

    assert scores == sorted(scores, reverse=True)
    assert page["sort_scope"] == "scanned window"
    assert page["note"]


def test_scan_is_bounded_for_a_filtered_query(table):
    from pdbclean.artefacts import MAX_FILTER_SCAN_ROWS

    page = parquet_page(table, page=1, page_size=5, search="nothing-matches")

    assert page["matched_rows"] == 0
    assert page["scanned_rows"] <= MAX_FILTER_SCAN_ROWS


# --------------------------------------------------------------------------
# Other formats
# --------------------------------------------------------------------------


def test_csv_pagination(tmp_path):
    target = tmp_path / "rows.csv"
    target.write_text(
        "a,b\n" + "".join(f"{i},{i * 2}\n" for i in range(200)),
        encoding="utf-8",
    )

    page = csv_page(target, page=2, page_size=10)

    assert page["columns"] == ["a", "b"]
    assert page["returned"] == 10
    assert page["rows"][0] == ["10", "20"]


def test_tsv_is_detected(tmp_path):
    target = tmp_path / "rows.tsv"
    target.write_text("a\tb\n1\t2\n", encoding="utf-8")

    page = csv_page(target, page=1, page_size=10)

    assert page["columns"] == ["a", "b"]
    assert page["rows"] == [["1", "2"]]


def test_table_page_dispatches_by_type(table, tmp_path):
    assert table_page(table, page=1, page_size=2)["kind"] == "parquet"

    csv_target = tmp_path / "x.csv"
    csv_target.write_text("a\n1\n", encoding="utf-8")

    assert table_page(csv_target, page=1, page_size=2)["kind"] == "table"


def test_unsupported_file_fails_safely(tmp_path):
    target = tmp_path / "weights.bin"
    target.write_bytes(b"\x00\x01\x02")

    with pytest.raises(ArtefactError):
        table_page(target)


def test_csv_export_is_a_view_not_the_artefact():
    rendered = rows_to_csv(
        ["a", "b"], [{"a": 1, "b": 2}, {"a": 3, "b": 4}]
    )

    assert rendered.splitlines()[0] == "a,b"
    assert rendered.splitlines()[1] == "1,2"


# --------------------------------------------------------------------------
# HTTP surface
# --------------------------------------------------------------------------


def test_artefact_endpoint_returns_schema_and_provenance(ui):
    target = _frozen(RELEASE / "audit" / "representative_mapping.parquet")

    payload = ui.get("/api/artefact", path=str(target))

    assert payload["name"] == "representative_mapping.parquet"
    assert payload["preview_kind"] == "parquet"
    assert payload["schema"]["row_count"] > 0
    assert payload["download_url"]
    assert payload["provenance"]["snapshot_id"] == "20260101"


def test_table_endpoint_paginates(ui):
    target = _frozen(RELEASE / "audit" / "removed_chain_audit.parquet")

    first = ui.get(
        "/api/artefact/table", path=str(target), page=1, page_size=3
    )
    second = ui.get(
        "/api/artefact/table", path=str(target), page=2, page_size=3
    )

    assert first["returned"] == 3
    assert second["returned"] == 3
    assert first["rows"] != second["rows"]
    assert first["row_count"] == 78_754


def test_download_returns_the_original_bytes(ui):
    target = _frozen(RELEASE / "audit" / "removed_chain_audit.parquet")

    status, body, headers = ui.raw(
        "/api/artefact/download", path=str(target)
    )

    assert status == 200
    assert len(body) == target.stat().st_size

    expected = hashlib.sha256(target.read_bytes()).hexdigest()

    assert hashlib.sha256(body).hexdigest() == expected
    assert 'filename="removed_chain_audit.parquet"' in headers[
        "Content-Disposition"
    ]


def test_downloaded_parquet_is_still_readable(ui, tmp_path):
    """The download is the artefact, not a re-encoded export."""

    target = _frozen(RELEASE / "audit" / "representatives.parquet")

    _, body, _ = ui.raw("/api/artefact/download", path=str(target))

    copied = tmp_path / "copy.parquet"
    copied.write_bytes(body)

    assert (
        pq.ParquetFile(str(copied)).metadata.num_rows
        == pq.ParquetFile(str(target)).metadata.num_rows
    )


def test_json_artefact_is_rendered(ui):
    target = _frozen(RELEASE / "release_manifest.json")

    payload = ui.get("/api/artefact", path=str(target))

    assert payload["preview_kind"] == "json"
    assert payload["preview"]["content"]["retained_chain_count"] == 499_770


def test_text_artefact_is_rendered(ui):
    target = _frozen(RELEASE / "_SUCCESS")

    payload = ui.get("/api/artefact", path=str(target))

    assert payload["preview_kind"] == "text"
    assert payload["preview"]["content"]


# --------------------------------------------------------------------------
# Security
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "path",
    [
        "/etc/passwd",
        "/etc/shadow",
        "../../../etc/passwd",
        "/root/.ssh/id_rsa",
    ],
)
def test_absolute_and_traversal_paths_are_refused(ui, path):
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        ui.get("/api/artefact", path=path)

    assert excinfo.value.code in {400, 403, 404}


@pytest.mark.parametrize(
    "relative",
    [
        "README.md",
        "src/pdbclean/cli.py",
        "pyproject.toml",
        "config/pdbclean/protocol_3_2_comp702_v1.yaml",
    ],
)
def test_repository_source_is_not_browsable(ui, relative):
    """This is a provenance viewer, not a filesystem browser."""

    with pytest.raises(urllib.error.HTTPError) as excinfo:
        ui.get("/api/artefact", path=str(REPO / relative))

    assert excinfo.value.code == 403


def test_traversal_out_of_an_allowed_root_is_refused(ui):
    escape = str(RELEASE / ".." / ".." / ".." / "README.md")

    with pytest.raises(urllib.error.HTTPError) as excinfo:
        ui.get("/api/artefact", path=escape)

    assert excinfo.value.code == 403


def test_symlink_escape_is_refused(ui, tmp_path):
    """resolve() follows the link before the containment check."""

    secret = tmp_path / "secret.json"
    secret.write_text("{}", encoding="utf-8")

    link = RELEASE / "escape_link.json"

    if not RELEASE.is_dir():
        pytest.skip("frozen release not present")

    try:
        link.symlink_to(secret)
    except OSError:  # pragma: no cover - symlinks unavailable
        pytest.skip("cannot create symlinks here")

    try:
        with pytest.raises(urllib.error.HTTPError) as excinfo:
            ui.get("/api/artefact", path=str(link))

        assert excinfo.value.code == 403
    finally:
        link.unlink()


def test_download_is_subject_to_the_same_allowlist(ui):
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        ui.raw("/api/artefact/download", path="/etc/passwd")

    assert excinfo.value.code in {400, 403, 404}


def test_missing_artefact_is_not_found(ui):
    with pytest.raises(urllib.error.HTTPError) as excinfo:
        ui.get("/api/artefact", path=str(RELEASE / "absent.parquet"))

    assert excinfo.value.code == 404


# --------------------------------------------------------------------------
# Read-only guarantees
# --------------------------------------------------------------------------


def test_viewing_does_not_modify_the_artefact(ui):
    target = _frozen(RELEASE / "audit" / "removed_chain_audit.parquet")

    before = hashlib.sha256(target.read_bytes()).hexdigest()
    before_mtime = target.stat().st_mtime

    ui.get("/api/artefact", path=str(target))
    ui.get("/api/artefact/table", path=str(target), page=1, page_size=5)
    ui.raw("/api/artefact/download", path=str(target))

    assert hashlib.sha256(target.read_bytes()).hexdigest() == before
    assert target.stat().st_mtime == before_mtime


def test_viewing_does_not_modify_run_provenance(ui, tmp_path):
    from pdbclean.pipeline import plan_pipeline, record_plan_in_provenance
    from pdbclean.run_provenance import RunProvenance
    from pdbclean.runconfig import resolve_run_config

    resolved = resolve_run_config(
        config_path=REPO / PROFILE,
        overrides={"storage.run_root": str(tmp_path / "runs")},
    )

    run = RunProvenance.create(
        resolved=resolved, run_root=tmp_path / "runs", repo_root=REPO
    )
    record_plan_in_provenance(plan_pipeline(resolved, repo_root=REPO), run)
    run.flush()

    record = run.run_dir / "run.json"
    events = run.run_dir / "events.jsonl"

    before_record = hashlib.sha256(record.read_bytes()).hexdigest()
    before_events = hashlib.sha256(events.read_bytes()).hexdigest()

    target = _frozen(RELEASE / "audit" / "representatives.parquet")

    ui.get("/api/artefact", path=str(target))
    ui.get("/api/artefact/table", path=str(target), page=1, page_size=5)

    assert hashlib.sha256(record.read_bytes()).hexdigest() == before_record
    assert hashlib.sha256(events.read_bytes()).hexdigest() == before_events
