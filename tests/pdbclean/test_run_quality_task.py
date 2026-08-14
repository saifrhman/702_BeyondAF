"""Tests for the quality-task command-line entrypoint."""

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


SCRIPT_PATH = Path("scripts/pdbclean/run_quality_task.py")


def _load_cli_module():
    spec = importlib.util.spec_from_file_location(
        "pdbclean_run_quality_task_test",
        SCRIPT_PATH,
    )
    assert spec is not None
    assert spec.loader is not None

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_latest_complete_config(
    tmp_path: Path,
    *,
    batch_size: int,
    minimum_backbone_distance_angstrom: float | None = None,
) -> Path:
    source = Path(
        "config/pdbclean/protocol_3_2_comp702_v1.yaml"
    )
    data = yaml.safe_load(source.read_text())

    data["snapshot"]["mode"] = "latest_complete"
    data["snapshot"].pop("snapshot_id", None)

    data["execution"]["batch_size"] = batch_size

    if minimum_backbone_distance_angstrom is not None:
        data["quality_rules"]["backbone_distance"][
            "minimum_distance_angstrom"
        ] = minimum_backbone_distance_angstrom

    data["storage"]["output_root"] = str(
        tmp_path / "pipeline-output"
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
) -> Path:
    rows = []

    for index in range(row_count):
        pdb_id = f"{index:04x}"[-4:]

        rows.append(
            {
                "snapshot": snapshot,
                "source_layout": "recursive_coordinate_files",
                "pdb_id": pdb_id,
                "s3_key": (
                    f"{snapshot}/coordinates/{pdb_id}.cif.gz"
                ),
                "size_bytes": 100 + index,
                "etag": f"etag-{index}",
            }
        )

    table = pa.Table.from_pylist(rows)

    path = tmp_path / "source_manifest.parquet"
    pq.write_table(table, path)

    return path


def test_cli_derives_partition_count_from_manifest_runtime_size(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    cli = _load_cli_module()

    config_path = _write_latest_complete_config(
        tmp_path,
        batch_size=2,
        minimum_backbone_distance_angstrom=0.02,
    )
    manifest_path = _write_manifest(
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
            task_id=2,
        ),
    )

    monkeypatch.setattr(
        cli,
        "resolve_clean_git_commit",
        lambda repository_root: "a" * 40,
    )

    observed = {}

    def execute_quality_task(rows, **kwargs):
        observed["rows"] = rows
        observed["kwargs"] = kwargs

        return SimpleNamespace(
            summary_path=(
                Path(kwargs["output_root"])
                / "summaries"
                / "task_2.json"
            )
        )

    monkeypatch.setattr(
        cli,
        "execute_quality_task",
        execute_quality_task,
    )

    cli.main()

    # 5 rows with batch_size 2 => tasks 0, 1, 2.
    # Task 2 receives only the final source row.
    assert len(observed["rows"]) == 1
    assert observed["rows"][0]["pdb_id"] == "0004"

    kwargs = observed["kwargs"]

    assert kwargs["task_id"] == 2
    assert kwargs["snapshot"] == "20310415"
    assert kwargs["cleaning_protocol"] == (
        "protocol3.2-comp702-v1"
    )
    assert kwargs["pipeline_git_commit"] == "a" * 40
    assert kwargs["max_retries"] == 3
    assert kwargs["download_concurrency"] == 4
    assert kwargs["timeout_seconds"] == 60
    assert kwargs["minimum_backbone_distance_angstrom"] == 0.02

    assert kwargs["output_root"] == (
        tmp_path
        / "pipeline-output"
        / "20310415"
        / "protocol3.2-comp702-v1"
        / "quality"
    )

    output = capsys.readouterr().out

    assert "Manifest rows: 5" in output
    assert "Partition: 2 / 2" in output
    assert "Partition rows: 1" in output


@pytest.mark.parametrize(
    ("row_count", "batch_size", "task_id", "expected_rows"),
    [
        (1, 500, 0, 1),
        (4, 2, 0, 2),
        (4, 2, 1, 2),
        (5, 2, 2, 1),
        (11, 4, 2, 3),
    ],
)
def test_cli_partition_size_changes_with_manifest_size(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    row_count: int,
    batch_size: int,
    task_id: int,
    expected_rows: int,
) -> None:
    cli = _load_cli_module()

    config_path = _write_latest_complete_config(
        tmp_path,
        batch_size=batch_size,
    )
    manifest_path = _write_manifest(
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
            task_id=task_id,
        ),
    )

    monkeypatch.setattr(
        cli,
        "resolve_clean_git_commit",
        lambda repository_root: "b" * 40,
    )

    observed = {}

    def execute_quality_task(rows, **kwargs):
        observed["row_count"] = len(rows)

        return SimpleNamespace(
            summary_path=tmp_path / "summary.json"
        )

    monkeypatch.setattr(
        cli,
        "execute_quality_task",
        execute_quality_task,
    )

    cli.main()

    assert observed["row_count"] == expected_rows


def test_cli_rejects_task_outside_dynamic_manifest_range(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cli = _load_cli_module()

    config_path = _write_latest_complete_config(
        tmp_path,
        batch_size=2,
    )
    manifest_path = _write_manifest(
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
            task_id=3,
        ),
    )

    with pytest.raises(
        ManifestError,
        match=r"task_id 3 is out of range",
    ):
        cli.main()
