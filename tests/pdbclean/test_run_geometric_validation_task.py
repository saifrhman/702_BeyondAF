"""Tests for the Step-2 geometric-validation CLI."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from pdbclean.geometric_validation import (
    GeometricValidationConfig,
)
from pdbclean.geometric_validation_runner import (
    GeometricValidationRunnerError,
)


SCRIPT_PATH = Path(
    "scripts/pdbclean/run_geometric_validation_task.py"
)


def _load_cli_module():
    spec = importlib.util.spec_from_file_location(
        "pdbclean_run_geometric_validation_task_test",
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
    batch_size: int = 2,
    distance: float = 0.02,
    angle: float = 5.0,
    enabled: bool = True,
) -> Path:
    source = Path(
        "config/pdbclean/protocol_3_2_comp702_v1.yaml"
    )

    data = yaml.safe_load(source.read_text())

    data["snapshot"]["mode"] = "latest_complete"
    data["snapshot"].pop("snapshot_id", None)

    data["execution"]["batch_size"] = batch_size
    data["storage"]["output_root"] = str(
        tmp_path / "pipeline-output"
    )

    data["quality_rules"]["backbone_distance"][
        "minimum_distance_angstrom"
    ] = distance

    data["post_cleaning_geometric_validation"][
        "enabled"
    ] = enabled

    data["post_cleaning_geometric_validation"][
        "minimum_triangle_angle_degrees"
    ] = angle

    path = tmp_path / "config.yaml"

    path.write_text(
        yaml.safe_dump(
            data,
            sort_keys=False,
        )
    )

    return path


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
                "source_layout": (
                    "recursive_coordinate_files"
                ),
                "pdb_id": pdb_id,
                "s3_key": (
                    f"{snapshot}/coordinates/{pdb_id}.cif.gz"
                ),
                "size_bytes": 100 + index,
                "etag": f"etag-{index}",
            }
        )

    path = tmp_path / "source_manifest.parquet"

    pq.write_table(
        pa.Table.from_pylist(rows),
        path,
    )

    return path


def _write_quality_task(
    tmp_path: Path,
    *,
    snapshot: str,
    task_id: int,
    partition_source_count: int,
    accepted_rows: list[dict],
) -> None:
    protocol = "protocol3.2-comp702-v1"

    quality_root = (
        tmp_path
        / "pipeline-output"
        / snapshot
        / protocol
        / "quality"
    )

    accepted_path = (
        quality_root
        / "accepted"
        / f"task_{task_id}.parquet"
    )
    summary_path = (
        quality_root
        / "summaries"
        / f"task_{task_id}.json"
    )

    accepted_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    summary_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    table = pa.Table.from_pylist(
        accepted_rows,
        schema=pa.schema(
            [
                pa.field(
                    "snapshot",
                    pa.string(),
                    nullable=False,
                ),
                pa.field(
                    "pdb_id",
                    pa.string(),
                    nullable=False,
                ),
                pa.field(
                    "model_id",
                    pa.int32(),
                    nullable=False,
                ),
                pa.field(
                    "label_chain_id",
                    pa.string(),
                    nullable=False,
                ),
                pa.field(
                    "cleaning_protocol",
                    pa.string(),
                    nullable=False,
                ),
                pa.field(
                    "pipeline_git_commit",
                    pa.string(),
                    nullable=False,
                ),
            ]
        ),
    )

    pq.write_table(
        table,
        accepted_path,
    )

    summary = {
        "summary_schema_name": (
            "pdbclean_quality_task_summary"
        ),
        "summary_schema_version": "1.0",
        "task_id": str(task_id),
        "snapshot": snapshot,
        "cleaning_protocol": protocol,
        "pipeline_git_commit": "q" * 40,
        "input_source_object_count": (
            partition_source_count
        ),
        "accepted_chain_count": len(
            accepted_rows
        ),
        "source_object_accounting_valid": True,
        "selected_chain_accounting_valid": True,
    }

    summary_path.write_text(
        json.dumps(summary),
        encoding="utf-8",
    )


def test_cli_resolves_dynamic_snapshot_and_upstream_quality_task(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    cli = _load_cli_module()

    config_path = _write_config(
        tmp_path,
        batch_size=2,
        distance=0.02,
        angle=5.0,
    )

    manifest_path = _write_manifest(
        tmp_path,
        snapshot="20310415",
        row_count=5,
    )

    # Task 2 is the final one-row partition (pdb_id 0004).
    accepted_row = {
        "snapshot": "20310415",
        "pdb_id": "0004",
        "model_id": 1,
        "label_chain_id": "A",
        "cleaning_protocol": (
            "protocol3.2-comp702-v1"
        ),
        "pipeline_git_commit": "q" * 40,
    }

    _write_quality_task(
        tmp_path,
        snapshot="20310415",
        task_id=2,
        partition_source_count=1,
        accepted_rows=[accepted_row],
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
        lambda repository_root: "g" * 40,
    )

    observed = {}

    def execute_geometric_validation_task(
        manifest_rows,
        accepted_rows,
        **kwargs,
    ):
        observed["manifest_rows"] = list(
            manifest_rows
        )
        observed["accepted_rows"] = list(
            accepted_rows
        )
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
        "execute_geometric_validation_task",
        execute_geometric_validation_task,
    )

    cli.main()

    assert len(observed["manifest_rows"]) == 1
    assert (
        observed["manifest_rows"][0]["pdb_id"]
        == "0004"
    )

    assert observed["accepted_rows"] == [
        accepted_row
    ]

    kwargs = observed["kwargs"]

    assert kwargs["task_id"] == 2
    assert kwargs["snapshot"] == "20310415"

    assert kwargs["config"] == (
        GeometricValidationConfig(
            minimum_backbone_distance_angstrom=0.02,
            minimum_triangle_angle_degrees=5.0,
        )
    )

    assert kwargs["cleaning_protocol"] == (
        "protocol3.2-comp702-v1"
    )
    assert (
        kwargs["quality_pipeline_git_commit"]
        == "q" * 40
    )
    assert (
        kwargs[
            "geometric_validation_pipeline_git_commit"
        ]
        == "g" * 40
    )

    assert kwargs["timeout_seconds"] == 60
    assert kwargs["max_retries"] == 3
    assert kwargs["download_concurrency"] == 4

    assert kwargs["output_root"] == (
        tmp_path
        / "pipeline-output"
        / "20310415"
        / "protocol3.2-comp702-v1"
        / "geometric_validation"
    )

    output = capsys.readouterr().out

    assert "Snapshot: 20310415" in output
    assert "Partition: 2 / 2" in output
    assert "Partition rows: 1" in output
    assert "Accepted input chains: 1" in output


def test_cli_requires_upstream_quality_completion_summary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cli = _load_cli_module()

    config_path = _write_config(
        tmp_path,
    )

    manifest_path = _write_manifest(
        tmp_path,
        snapshot="20310415",
        row_count=1,
    )

    monkeypatch.setattr(
        cli,
        "parse_args",
        lambda: argparse.Namespace(
            config=str(config_path),
            manifest=str(manifest_path),
            task_id=0,
        ),
    )

    with pytest.raises(
        FileNotFoundError,
        match="completion summary",
    ):
        cli.main()


def test_cli_refuses_when_geometric_validation_is_disabled(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cli = _load_cli_module()

    config_path = _write_config(
        tmp_path,
        enabled=False,
    )

    manifest_path = _write_manifest(
        tmp_path,
        snapshot="20310415",
        row_count=1,
    )

    monkeypatch.setattr(
        cli,
        "parse_args",
        lambda: argparse.Namespace(
            config=str(config_path),
            manifest=str(manifest_path),
            task_id=0,
        ),
    )

    with pytest.raises(
        GeometricValidationRunnerError,
        match="disabled",
    ):
        cli.main()
