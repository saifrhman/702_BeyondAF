"""Tests for the Stage-3 BRI finalization command-line entrypoint."""

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
    "scripts/pdbclean/finalize_bri.py"
)


def _load_cli_module():
    spec = importlib.util.spec_from_file_location(
        "pdbclean_finalize_bri_test",
        SCRIPT_PATH,
    )

    assert spec is not None
    assert spec.loader is not None

    module = importlib.util.module_from_spec(
        spec
    )

    spec.loader.exec_module(
        module
    )

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
                "source_layout":
                    "recursive_coordinate_files",
                "pdb_id": pdb_id,
                "s3_key": (
                    f"{snapshot}/coordinates/"
                    f"{pdb_id}.cif.gz"
                ),
                "size_bytes": 100 + index,
                "etag": f"etag-{index}",
            }
        )

    path = (
        tmp_path
        / "source_manifest.parquet"
    )

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
        "config/pdbclean/"
        "protocol_3_2_comp702_v1.yaml"
    )

    data = yaml.safe_load(
        source.read_text(
            encoding="utf-8"
        )
    )

    data["snapshot"]["mode"] = "fixed"
    data["snapshot"]["snapshot_id"] = (
        snapshot
    )
    data["snapshot"].pop(
        "expected_mmcif_count",
        None,
    )
    data["snapshot"].pop(
        "expected_total_bytes",
        None,
    )

    data["execution"]["batch_size"] = (
        batch_size
    )

    data["storage"]["output_root"] = str(
        tmp_path
        / "pipeline-output"
    )

    path = (
        tmp_path
        / "config.yaml"
    )

    path.write_text(
        yaml.safe_dump(
            data,
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    return path


def test_resolve_bri_producer_commit_accepts_unique_commit(
    tmp_path: Path,
) -> None:
    cli = _load_cli_module()

    summary_dir = (
        tmp_path
        / "bri"
        / "summaries"
    )

    summary_dir.mkdir(
        parents=True
    )

    for task_id in (0, 1):
        (
            summary_dir
            / f"task_{task_id}.json"
        ).write_text(
            json.dumps(
                {
                    "bri_pipeline_git_commit":
                        "d" * 40,
                }
            )
            + "\n",
            encoding="utf-8",
        )

    assert cli.resolve_bri_producer_commit(
        tmp_path / "bri"
    ) == "d" * 40


def test_resolve_bri_producer_commit_rejects_mixed_commits(
    tmp_path: Path,
) -> None:
    cli = _load_cli_module()

    summary_dir = (
        tmp_path
        / "bri"
        / "summaries"
    )

    summary_dir.mkdir(
        parents=True
    )

    (
        summary_dir
        / "task_0.json"
    ).write_text(
        json.dumps(
            {
                "bri_pipeline_git_commit":
                    "d" * 40,
            }
        )
        + "\n",
        encoding="utf-8",
    )

    (
        summary_dir
        / "task_1.json"
    ).write_text(
        json.dumps(
            {
                "bri_pipeline_git_commit":
                    "e" * 40,
            }
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(
        RuntimeError,
        match="multiple producer Git commits",
    ):
        cli.resolve_bri_producer_commit(
            tmp_path / "bri"
        )


def test_cli_forwards_dynamic_manifest_and_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cli = _load_cli_module()

    snapshot = "20310415"
    protocol = "protocol3.2-comp702-v1"

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

    quality_root = (
        tmp_path
        / "pipeline-output"
        / snapshot
        / protocol
        / "quality"
    )

    geometric_validation_root = (
        quality_root.parent
        / "geometric_validation"
    )

    eligible_path = (
        geometric_validation_root
        / "finalized"
        / "eligible.parquet"
    )

    bri_root = (
        quality_root.parent
        / "bri"
    )

    monkeypatch.setattr(
        cli,
        "parse_args",
        lambda: argparse.Namespace(
            config=str(config_path),
            manifest=str(manifest_path),
        ),
    )

    upstream = SimpleNamespace(
        eligible_path=eligible_path,
        quality_pipeline_git_commit=(
            "a" * 40
        ),
        geometric_validation_pipeline_git_commit=(
            "b" * 40
        ),
        geometric_validation_finalizer_git_commit=(
            "c" * 40
        ),
    )

    observed_upstream = {}

    def validate_upstream(
        root,
        **kwargs,
    ):
        observed_upstream["root"] = root
        observed_upstream.update(
            kwargs
        )

        return upstream

    monkeypatch.setattr(
        cli,
        "validate_upstream_geometric_validation_stage",
        validate_upstream,
    )

    monkeypatch.setattr(
        cli,
        "resolve_bri_producer_commit",
        lambda root: "d" * 40,
    )

    monkeypatch.setattr(
        cli,
        "resolve_clean_git_commit",
        lambda repository_root: "e" * 40,
    )

    observed_finalizer = {}

    def finalize_bri_stage(**kwargs):
        observed_finalizer.update(
            kwargs
        )

        return SimpleNamespace(
            bri_path=Path(
                "finalized/bri.parquet"
            ),
            global_summary_path=Path(
                "global_summary.json"
            ),
            success_path=Path(
                "_SUCCESS"
            ),
            global_summary={
                "task_count": 3,
                "input_eligible_chain_count": 10,
                "bri_chain_count": 10,
                "processing_error_count": 0,
                "unique_bri_identity_count": 10,
                "minimum_retained_residue_count": 1,
                "maximum_retained_residue_count": 20,
            },
        )

    monkeypatch.setattr(
        cli,
        "finalize_bri_stage",
        finalize_bri_stage,
    )

    cli.main()

    assert observed_upstream[
        "root"
    ] == geometric_validation_root

    assert observed_upstream[
        "expected_snapshot"
    ] == snapshot

    assert observed_upstream[
        "expected_cleaning_protocol"
    ] == protocol

    # Five manifest objects / batch size two = three tasks.
    assert observed_upstream[
        "expected_task_count"
    ] == 3

    assert observed_finalizer[
        "bri_root"
    ] == bri_root

    assert observed_finalizer[
        "eligible_path"
    ] == eligible_path

    assert observed_finalizer[
        "manifest_row_count"
    ] == 5

    assert observed_finalizer[
        "batch_size"
    ] == 2

    assert observed_finalizer[
        "snapshot"
    ] == snapshot

    assert observed_finalizer[
        "cleaning_protocol"
    ] == protocol

    assert observed_finalizer[
        "quality_pipeline_git_commit"
    ] == "a" * 40

    assert observed_finalizer[
        "geometric_validation_pipeline_git_commit"
    ] == "b" * 40

    assert observed_finalizer[
        "geometric_validation_finalizer_git_commit"
    ] == "c" * 40

    assert observed_finalizer[
        "bri_pipeline_git_commit"
    ] == "d" * 40

    assert observed_finalizer[
        "finalizer_pipeline_git_commit"
    ] == "e" * 40
