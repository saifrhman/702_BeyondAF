"""Tests for the Step-2 finalization command-line entrypoint."""

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


SCRIPT_PATH = Path(
    "scripts/pdbclean/finalize_geometric_validation.py"
)


def _load_cli_module():
    spec = importlib.util.spec_from_file_location(
        "pdbclean_finalize_geometry_test",
        SCRIPT_PATH,
    )
    assert spec is not None
    assert spec.loader is not None

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


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

    path = tmp_path / "source_manifest.parquet"

    pq.write_table(
        pa.Table.from_pylist(rows),
        path,
    )

    return path


def _write_config(
    tmp_path: Path,
    *,
    snapshot: str,
    batch_size: int,
) -> Path:
    source = Path(
        "config/pdbclean/protocol_3_2_comp702_v1.yaml"
    )

    data = yaml.safe_load(
        source.read_text(encoding="utf-8")
    )

    data["snapshot"]["mode"] = "fixed"
    data["snapshot"]["snapshot_id"] = snapshot
    data["snapshot"].pop(
        "expected_mmcif_count",
        None,
    )
    data["snapshot"].pop(
        "expected_total_bytes",
        None,
    )

    data["execution"]["batch_size"] = batch_size
    data["storage"]["output_root"] = str(
        tmp_path / "pipeline-output"
    )

    path = tmp_path / "config.yaml"

    path.write_text(
        yaml.safe_dump(data, sort_keys=False),
        encoding="utf-8",
    )

    return path


def test_cli_forwards_dynamic_manifest_and_configuration(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cli = _load_cli_module()

    snapshot = "20310415"

    config_path = _write_config(
        tmp_path,
        snapshot=snapshot,
        batch_size=2,
    )

    manifest_path = _write_manifest(
        tmp_path,
        snapshot=snapshot,
        row_count=5,
    )

    protocol = "protocol3.2-comp702-v1"

    quality_root = (
        tmp_path
        / "pipeline-output"
        / snapshot
        / protocol
        / "quality"
    )

    quality_root.mkdir(parents=True)
    (quality_root / "_SUCCESS").write_text(
        "{}\n",
        encoding="utf-8",
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
        "resolve_stage2_producer_commits",
        lambda root: ("a" * 40, "b" * 40),
    )

    monkeypatch.setattr(
        cli,
        "resolve_clean_git_commit",
        lambda repository_root: "c" * 40,
    )

    observed = {}

    def finalize_geometric_validation_stage(**kwargs):
        observed.update(kwargs)

        return SimpleNamespace(
            eligible_path=Path("eligible.parquet"),
            quarantined_path=Path(
                "quarantined.parquet"
            ),
            global_summary_path=Path(
                "global_summary.json"
            ),
            success_path=Path("_SUCCESS"),
            global_summary={
                "task_count": 3,
                "input_accepted_chain_count": 10,
                "eligible_chain_count": 9,
                "quarantined_chain_count": 1,
                "processing_error_count": 0,
                "violation_event_count": 1,
            },
        )

    monkeypatch.setattr(
        cli,
        "finalize_geometric_validation_stage",
        finalize_geometric_validation_stage,
    )

    cli.main()

    assert observed["manifest_row_count"] == 5
    assert observed["batch_size"] == 2
    assert observed["snapshot"] == snapshot
    assert observed["quality_root"] == quality_root
    assert observed[
        "geometric_validation_root"
    ] == quality_root.parent / "geometric_validation"

    assert observed[
        "quality_pipeline_git_commit"
    ] == "a" * 40

    assert observed[
        "geometric_validation_pipeline_git_commit"
    ] == "b" * 40

    assert observed[
        "finalizer_pipeline_git_commit"
    ] == "c" * 40

    assert observed[
        "minimum_backbone_distance_angstrom"
    ] == 0.01

    assert observed[
        "minimum_triangle_angle_degrees"
    ] == 3.0


def test_resolve_stage2_commits_rejects_mixed_geometry_commits(
    tmp_path: Path,
) -> None:
    cli = _load_cli_module()

    root = tmp_path / "geometry"
    summaries = root / "summaries"
    summaries.mkdir(parents=True)

    for task_id, geometry_commit in (
        (0, "b" * 40),
        (1, "c" * 40),
    ):
        (summaries / f"task_{task_id}.json").write_text(
            json.dumps(
                {
                    "quality_pipeline_git_commit": (
                        "a" * 40
                    ),
                    "geometric_validation_pipeline_git_commit": (
                        geometry_commit
                    ),
                }
            ),
            encoding="utf-8",
        )

    with pytest.raises(
        cli.GeometricValidationFinalizeError,
        match="multiple geometry producer",
    ):
        cli.resolve_stage2_producer_commits(root)


def test_cli_requires_completed_quality_stage(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cli = _load_cli_module()

    snapshot = "20310415"

    config_path = _write_config(
        tmp_path,
        snapshot=snapshot,
        batch_size=2,
    )

    manifest_path = _write_manifest(
        tmp_path,
        snapshot=snapshot,
        row_count=1,
    )

    monkeypatch.setattr(
        cli,
        "parse_args",
        lambda: argparse.Namespace(
            config=str(config_path),
            manifest=str(manifest_path),
        ),
    )

    with pytest.raises(
        cli.GeometricValidationFinalizeError,
        match="no _SUCCESS marker",
    ):
        cli.main()


def test_cli_contains_no_production_snapshot_or_task_count() -> None:
    source = SCRIPT_PATH.read_text(
        encoding="utf-8"
    )

    assert "20260101" not in source
    assert "494" not in source
