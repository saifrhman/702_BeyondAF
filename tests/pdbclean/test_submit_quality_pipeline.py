"""Tests for dynamic Slurm submission of the PDBClean quality stage."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml


SUBMIT_SCRIPT = Path(
    "scripts/pdbclean/submit_quality_pipeline.sh"
)
QUALITY_SCRIPT = Path(
    "scripts/pdbclean/run_quality_array.sbatch"
)
MERGE_SCRIPT = Path(
    "scripts/pdbclean/merge_quality_outputs.sbatch"
)


def _write_config(
    tmp_path: Path,
    *,
    batch_size: int,
) -> Path:
    source = Path(
        "config/pdbclean/protocol_3_2_comp702_v1.yaml"
    )
    data = yaml.safe_load(source.read_text())

    data["snapshot"]["mode"] = "latest_complete"
    data["snapshot"].pop("snapshot_id", None)
    data["snapshot"].pop("expected_mmcif_count", None)
    data["snapshot"].pop("expected_total_bytes", None)

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
                    f"{snapshot}/coordinates/"
                    f"{pdb_id}.cif.gz"
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


def _prepare_clean_repo(tmp_path: Path) -> tuple[Path, Path]:
    repo = tmp_path / "submission-repo"
    scripts = repo / "scripts" / "pdbclean"
    scripts.mkdir(parents=True)

    for source in (
        SUBMIT_SCRIPT,
        QUALITY_SCRIPT,
        MERGE_SCRIPT,
    ):
        destination = scripts / source.name
        shutil.copy2(source, destination)
        destination.chmod(0o755)

    subprocess.run(
        ["git", "init", "-q"],
        cwd=repo,
        check=True,
    )
    subprocess.run(
        ["git", "config", "user.email", "test@example.com"],
        cwd=repo,
        check=True,
    )
    subprocess.run(
        ["git", "config", "user.name", "PDBClean Test"],
        cwd=repo,
        check=True,
    )
    subprocess.run(
        ["git", "add", "scripts"],
        cwd=repo,
        check=True,
    )
    subprocess.run(
        ["git", "commit", "-q", "-m", "test fixture"],
        cwd=repo,
        check=True,
    )

    return repo, scripts / SUBMIT_SCRIPT.name


def _write_fake_sbatch(tmp_path: Path) -> tuple[Path, Path, Path]:
    fake_bin = tmp_path / "fake-bin"
    fake_bin.mkdir()

    log_path = tmp_path / "sbatch.log"
    counter_path = tmp_path / "sbatch.counter"

    fake = fake_bin / "sbatch"
    fake.write_text(
        """#!/bin/bash
set -euo pipefail

printf '%s\n' "$*" >> "$FAKE_SBATCH_LOG"

count=0
if [[ -f "$FAKE_SBATCH_COUNTER" ]]; then
    count="$(cat "$FAKE_SBATCH_COUNTER")"
fi

count=$((count + 1))
printf '%s\n' "$count" > "$FAKE_SBATCH_COUNTER"

if [[ "$count" -eq 1 ]]; then
    echo "12345;fakecluster"
else
    echo "12346;fakecluster"
fi
""",
        encoding="utf-8",
    )
    fake.chmod(0o755)

    return fake_bin, log_path, counter_path


@pytest.mark.parametrize(
    ("row_count", "batch_size", "expected_array"),
    [
        (1, 500, "0-0%24"),
        (5, 2, "0-2%24"),
        (11, 4, "0-2%24"),
        (37, 7, "0-5%24"),
    ],
)
def test_submission_derives_array_and_afterok_dependency(
    tmp_path: Path,
    row_count: int,
    batch_size: int,
    expected_array: str,
) -> None:
    config = _write_config(
        tmp_path,
        batch_size=batch_size,
    )
    manifest = _write_manifest(
        tmp_path,
        snapshot="20400101",
        row_count=row_count,
    )

    repo, wrapper = _prepare_clean_repo(tmp_path)
    fake_bin, log_path, counter_path = (
        _write_fake_sbatch(tmp_path)
    )

    home = tmp_path / "home"
    (home / "fastscratch").mkdir(parents=True)

    env = os.environ.copy()
    env["HOME"] = str(home)
    env["PATH"] = (
        str(fake_bin)
        + os.pathsep
        + env["PATH"]
    )
    env["FAKE_SBATCH_LOG"] = str(log_path)
    env["FAKE_SBATCH_COUNTER"] = str(counter_path)

    result = subprocess.run(
        [
            "bash",
            str(wrapper),
            str(config),
            str(manifest),
        ],
        cwd=repo,
        env=env,
        text=True,
        capture_output=True,
        check=True,
    )

    lines = log_path.read_text(
        encoding="utf-8"
    ).splitlines()

    assert len(lines) == 2

    array_args = shlex.split(lines[0])
    merge_args = shlex.split(lines[1])

    assert "--parsable" in array_args
    assert f"--array={expected_array}" in array_args
    assert any(
        value.endswith("run_quality_array.sbatch")
        for value in array_args
    )

    assert "--parsable" in merge_args
    assert "--dependency=afterok:12345" in merge_args
    assert any(
        value.endswith("merge_quality_outputs.sbatch")
        for value in merge_args
    )

    expected_task_count = (
        row_count + batch_size - 1
    ) // batch_size

    assert (
        f"Task count:     {expected_task_count}"
        in result.stdout
    )
    assert "Concurrency:    24" in result.stdout
    assert (
        f"Array range:    {expected_array}"
        in result.stdout
    )
    assert (
        "Quality array submitted: 12345"
        in result.stdout
    )
    assert (
        "Quality merge submitted: 12346"
        in result.stdout
    )
    assert (
        "Dependency: afterok:12345"
        in result.stdout
    )


def test_submission_scripts_do_not_hardcode_494() -> None:
    for path in (
        SUBMIT_SCRIPT,
        QUALITY_SCRIPT,
        MERGE_SCRIPT,
    ):
        assert "494" not in path.read_text(
            encoding="utf-8"
        )
