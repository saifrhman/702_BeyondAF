"""Tests for the quality-merge command-line entrypoint."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from pdbclean.manifest import ManifestError


SCRIPT_PATH = Path(
    "scripts/pdbclean/merge_quality_outputs.py"
)


def _load_cli_module():
    spec = importlib.util.spec_from_file_location(
        "pdbclean_merge_quality_outputs_test",
        SCRIPT_PATH,
    )
    assert spec is not None
    assert spec.loader is not None

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_config(
    tmp_path: Path,
    *,
    batch_size: int,
    snapshot_mode: str = "latest_complete",
    snapshot_id: str | None = None,
    expected_count: int | None = None,
    expected_total_bytes: int | None = None,
) -> Path:
    source = Path(
        "config/pdbclean/protocol_3_2_comp702_v1.yaml"
    )
    data = yaml.safe_load(source.read_text())

    data["execution"]["batch_size"] = batch_size
    data["storage"]["output_root"] = str(
        tmp_path / "pipeline-output"
    )

    snapshot = data["snapshot"]
    snapshot["mode"] = snapshot_mode

    if snapshot_mode == "latest_complete":
        snapshot.pop("snapshot_id", None)
        snapshot.pop("expected_mmcif_count", None)
        snapshot.pop("expected_total_bytes", None)
    else:
        assert snapshot_id is not None
        snapshot["snapshot_id"] = snapshot_id

        if expected_count is None:
            snapshot.pop("expected_mmcif_count", None)
        else:
            snapshot["expected_mmcif_count"] = expected_count

        if expected_total_bytes is None:
            snapshot.pop("expected_total_bytes", None)
        else:
            snapshot["expected_total_bytes"] = (
                expected_total_bytes
            )

    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    return config_path


def _write_manifest(
    tmp_path: Path,
    *,
    snapshot: str,
    row_count: int,
) -> tuple[Path, int]:
    rows = []

    for index in range(row_count):
        pdb_id = f"{index:04x}"[-4:]
        size_bytes = 100 + index

        rows.append(
            {
                "snapshot": snapshot,
                "source_layout": (
                    "recursive_coordinate_files"
                ),
                "pdb_id": pdb_id,
                "s3_key": (
                    f"{snapshot}/coordinates/"
                    f"{pdb_id}.cif.gz"
                ),
                "size_bytes": size_bytes,
                "etag": f"etag-{index}",
            }
        )

    table = pa.Table.from_pylist(rows)

    path = tmp_path / "source_manifest.parquet"
    pq.write_table(table, path)

    return path, sum(
        row["size_bytes"] for row in rows
    )


def _publication(
    quality_root: Path,
    *,
    task_count: int,
) -> SimpleNamespace:
    return SimpleNamespace(
        global_summary={
            "task_count": task_count,
            "accepted_chain_count": 10,
            "rejected_chain_count": 2,
            "non_candidate_chain_count": 3,
            "processing_error_count": 1,
        },
        global_summary_path=(
            quality_root / "global_summary.json"
        ),
        success_path=quality_root / "_SUCCESS",
    )


def test_cli_forwards_runtime_manifest_and_batch_size(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    cli = _load_cli_module()

    config_path = _write_config(
        tmp_path,
        batch_size=2,
    )
    manifest_path, _ = _write_manifest(
        tmp_path,
        snapshot="20310415",
        row_count=5,
    )

    monkeypatch.setattr(
        cli,
        "parse_args",
        lambda: argparse.Namespace(
            config=str(config_path),
            manifest=str(manifest_path),
        ),
    )

    monkeypatch.setattr(
        cli,
        "resolve_clean_git_commit",
        lambda repository_root: "a" * 40,
    )

    observed = {}

    def merge_quality_stage(**kwargs):
        observed.update(kwargs)

        return _publication(
            Path(kwargs["quality_root"]),
            task_count=3,
        )

    monkeypatch.setattr(
        cli,
        "merge_quality_stage",
        merge_quality_stage,
    )

    cli.main()

    assert observed["manifest"].num_rows == 5
    assert observed["manifest_row_count"] == 5
    assert observed["batch_size"] == 2
    assert observed["snapshot"] == "20310415"
    assert observed["cleaning_protocol"] == (
        "protocol3.2-comp702-v1"
    )
    assert observed["pipeline_git_commit"] == "a" * 40

    assert observed["quality_root"] == (
        tmp_path
        / "pipeline-output"
        / "20310415"
        / "protocol3.2-comp702-v1"
        / "quality"
    )

    output = capsys.readouterr().out

    assert "Manifest rows: 5" in output
    assert "Completed tasks: 3" in output
    assert "Accepted chains: 10" in output
    assert "Rejected chains: 2" in output
    assert "Non-candidate chains: 3" in output
    assert "Processing errors: 1" in output
    assert "Stage completion:" in output


@pytest.mark.parametrize(
    ("row_count", "batch_size"),
    [
        (1, 500),
        (4, 2),
        (5, 2),
        (11, 4),
        (37, 7),
    ],
)
def test_cli_passes_dynamic_manifest_size_to_merger(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    row_count: int,
    batch_size: int,
) -> None:
    cli = _load_cli_module()

    config_path = _write_config(
        tmp_path,
        batch_size=batch_size,
    )
    manifest_path, _ = _write_manifest(
        tmp_path,
        snapshot="20400101",
        row_count=row_count,
    )

    monkeypatch.setattr(
        cli,
        "parse_args",
        lambda: argparse.Namespace(
            config=str(config_path),
            manifest=str(manifest_path),
        ),
    )

    monkeypatch.setattr(
        cli,
        "resolve_clean_git_commit",
        lambda repository_root: "b" * 40,
    )

    observed = {}

    def merge_quality_stage(**kwargs):
        observed["manifest_row_count"] = kwargs[
            "manifest_row_count"
        ]
        observed["batch_size"] = kwargs["batch_size"]

        task_count = (
            row_count + batch_size - 1
        ) // batch_size

        return _publication(
            Path(kwargs["quality_root"]),
            task_count=task_count,
        )

    monkeypatch.setattr(
        cli,
        "merge_quality_stage",
        merge_quality_stage,
    )

    cli.main()

    assert observed["manifest_row_count"] == row_count
    assert observed["batch_size"] == batch_size


def test_cli_fixed_snapshot_validates_expected_manifest_totals(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cli = _load_cli_module()

    manifest_path, total_bytes = _write_manifest(
        tmp_path,
        snapshot="20310415",
        row_count=3,
    )

    config_path = _write_config(
        tmp_path,
        batch_size=2,
        snapshot_mode="fixed",
        snapshot_id="20310415",
        expected_count=4,
        expected_total_bytes=total_bytes,
    )

    monkeypatch.setattr(
        cli,
        "parse_args",
        lambda: argparse.Namespace(
            config=str(config_path),
            manifest=str(manifest_path),
        ),
    )

    called = False

    def merge_quality_stage(**kwargs):
        nonlocal called
        called = True
        raise AssertionError("merge must not run")

    monkeypatch.setattr(
        cli,
        "merge_quality_stage",
        merge_quality_stage,
    )

    with pytest.raises(
        ManifestError,
        match=r"Expected 4 rows, found 3",
    ):
        cli.main()

    assert called is False


def test_merge_cli_contains_no_hardcoded_production_task_count() -> None:
    source = SCRIPT_PATH.read_text()

    assert "494" not in source
