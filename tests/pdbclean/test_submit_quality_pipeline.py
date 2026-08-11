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



def test_quality_worker_requests_16g_memory() -> None:
    text = QUALITY_SCRIPT.read_text(encoding="utf-8")

    assert "#SBATCH --mem=16G" in text



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
    (
        "row_count",
        "batch_size",
        "expected_logical_tasks",
        "expected_workers",
        "expected_array",
    ),
    [
        (1, 500, 1, 1, "0-0%1"),
        (5, 2, 3, 3, "0-2%3"),
        (11, 4, 3, 3, "0-2%3"),
        (37, 7, 6, 6, "0-5%4"),
        (130, 1, 130, 64, "0-63%4"),
    ],
)
def test_submission_derives_array_and_afterok_dependency(
    tmp_path: Path,
    row_count: int,
    batch_size: int,
    expected_logical_tasks: int,
    expected_workers: int,
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

    worker_script_index = next(
        index
        for index, value in enumerate(array_args)
        if value.endswith("run_quality_array.sbatch")
    )

    assert array_args[worker_script_index + 3] == str(
        expected_logical_tasks
    )
    assert array_args[worker_script_index + 4] == str(
        expected_workers
    )

    repository_root = Path(
        array_args[worker_script_index + 5]
    )
    assert repository_root.is_absolute()
    assert repository_root == repo.resolve()

    assert "--parsable" in merge_args
    assert "--dependency=afterok:12345" in merge_args

    merge_script_index = next(
        index
        for index, value in enumerate(merge_args)
        if value.endswith("merge_quality_outputs.sbatch")
    )

    merge_repository_root = Path(
        merge_args[merge_script_index + 3]
    )
    assert merge_repository_root.is_absolute()
    assert merge_repository_root == repo.resolve()

    assert (
        f"Logical tasks:  {expected_logical_tasks}"
        in result.stdout
    )
    assert (
        f"Workers:        {expected_workers}"
        in result.stdout
    )
    expected_concurrency = min(4, expected_workers)
    assert (
        f"Concurrency:    {expected_concurrency}"
        in result.stdout
    )
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


def test_physical_worker_processes_expected_logical_stride(
    tmp_path: Path,
) -> None:
    fake_bin = tmp_path / "fake-bin"
    fake_bin.mkdir()

    python_log = tmp_path / "python.log"

    fake_module = fake_bin / "module"
    fake_module.write_text(
        "#!/bin/bash\nexit 0\n",
        encoding="utf-8",
    )
    fake_module.chmod(0o755)

    fake_conda = fake_bin / "conda"
    fake_conda.write_text(
        """#!/bin/bash
set -euo pipefail
if [[ "${1:-}" == "shell.bash" && "${2:-}" == "hook" ]]; then
    echo ":"
    exit 0
fi
exit 0
""",
        encoding="utf-8",
    )
    fake_conda.chmod(0o755)

    fake_python = fake_bin / "python"
    fake_python.write_text(
        """#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "$FAKE_PYTHON_LOG"
""",
        encoding="utf-8",
    )
    fake_python.chmod(0o755)

    config = tmp_path / "config.yaml"
    config.write_text("test\n", encoding="utf-8")

    manifest = tmp_path / "manifest.parquet"
    manifest.write_text("test\n", encoding="utf-8")

    env = os.environ.copy()

    # Barkla defines module/conda through shell initialisation. Force a
    # controlled non-interactive shell environment so the worker cannot
    # alter PATH and accidentally invoke the real Python environment.
    bash_env = tmp_path / "bash_env.sh"
    bash_env.write_text(
        """module() {
    return 0
}

conda() {
    if [[ "${1:-}" == "shell.bash" && "${2:-}" == "hook" ]]; then
        printf ':\\n'
    fi
    return 0
}
""",
        encoding="utf-8",
    )

    for key in tuple(env):
        if (
            key.startswith("BASH_FUNC_conda")
            or key.startswith("BASH_FUNC_module")
        ):
            env.pop(key)

    env["BASH_ENV"] = str(bash_env)
    env["PATH"] = (
        str(fake_bin)
        + os.pathsep
        + env["PATH"]
    )
    env["FAKE_PYTHON_LOG"] = str(python_log)
    env["SLURM_JOB_ID"] = "12345"
    env["SLURM_ARRAY_JOB_ID"] = "12345"
    env["SLURM_ARRAY_TASK_ID"] = "3"

    # Slurm executes a copied batch script from its spool area rather
    # than from the repository. Reproduce that condition explicitly.
    spool_dir = tmp_path / "var" / "spool" / "slurmd"
    spool_dir.mkdir(parents=True)
    spool_worker = spool_dir / "run_quality_array.sbatch"
    spool_worker.write_text(
        QUALITY_SCRIPT.read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    spool_worker.chmod(0o755)

    repository_root = QUALITY_SCRIPT.resolve().parents[2]

    result = subprocess.run(
        [
            "bash",
            str(spool_worker),
            str(config),
            str(manifest),
            "130",
            "64",
            str(repository_root),
        ],
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"worker stdout:\n{result.stdout}\n"
        f"worker stderr:\n{result.stderr}"
    )
    assert f"Repository:           {repository_root}" in result.stdout
    assert str(spool_dir) not in (
        result.stdout.split("Repository:", 1)[1].splitlines()[0]
    )

    calls = python_log.read_text(
        encoding="utf-8"
    ).splitlines()

    task_ids = []

    for call in calls:
        args = shlex.split(call)
        index = args.index("--task-id")
        task_ids.append(int(args[index + 1]))

    assert task_ids == [3, 67]


def test_physical_worker_striding_covers_logical_tasks_once() -> None:
    logical_task_count = 494
    physical_worker_count = 64

    assignments = [
        logical_task_id
        for worker_id in range(physical_worker_count)
        for logical_task_id in range(
            worker_id,
            logical_task_count,
            physical_worker_count,
        )
    ]

    assert len(assignments) == logical_task_count
    assert len(set(assignments)) == logical_task_count
    assert sorted(assignments) == list(
        range(logical_task_count)
    )

    per_worker_counts = [
        len(
            range(
                worker_id,
                logical_task_count,
                physical_worker_count,
            )
        )
        for worker_id in range(physical_worker_count)
    ]

    assert min(per_worker_counts) == 7
    assert max(per_worker_counts) == 8


def test_merge_worker_ignores_slurm_spool_location(
    tmp_path: Path,
) -> None:
    fake_bin = tmp_path / "fake-bin"
    fake_bin.mkdir()

    python_log = tmp_path / "python.log"

    fake_python = fake_bin / "python"
    fake_python.write_text(
        """#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "$FAKE_PYTHON_LOG"
""",
        encoding="utf-8",
    )
    fake_python.chmod(0o755)

    bash_env = tmp_path / "bash_env.sh"
    bash_env.write_text(
        """module() {
    return 0
}

conda() {
    if [[ "${1:-}" == "shell.bash" && "${2:-}" == "hook" ]]; then
        printf ':\\n'
    fi
    return 0
}
""",
        encoding="utf-8",
    )

    config = tmp_path / "config.yaml"
    config.write_text("test\n", encoding="utf-8")

    manifest = tmp_path / "manifest.parquet"
    manifest.write_text("test\n", encoding="utf-8")

    spool_dir = tmp_path / "var" / "spool" / "slurmd"
    spool_dir.mkdir(parents=True)

    spool_merge = spool_dir / "merge_quality_outputs.sbatch"
    spool_merge.write_text(
        MERGE_SCRIPT.read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    spool_merge.chmod(0o755)

    repository_root = MERGE_SCRIPT.resolve().parents[2]

    env = os.environ.copy()

    for key in tuple(env):
        if (
            key.startswith("BASH_FUNC_conda")
            or key.startswith("BASH_FUNC_module")
        ):
            env.pop(key)

    env["BASH_ENV"] = str(bash_env)
    env["PATH"] = (
        str(fake_bin)
        + os.pathsep
        + env["PATH"]
    )
    env["FAKE_PYTHON_LOG"] = str(python_log)
    env["SLURM_JOB_ID"] = "54321"

    result = subprocess.run(
        [
            "bash",
            str(spool_merge),
            str(config),
            str(manifest),
            str(repository_root),
        ],
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"merge stdout:\n{result.stdout}\n"
        f"merge stderr:\n{result.stderr}"
    )

    assert (
        f"Repository: {repository_root}"
        in result.stdout
    )

    calls = python_log.read_text(
        encoding="utf-8"
    ).splitlines()

    assert len(calls) == 1

    args = shlex.split(calls[0])

    assert args[0] == (
        "scripts/pdbclean/merge_quality_outputs.py"
    )
    assert args[1:3] == [
        "--config",
        str(config.resolve()),
    ]
    assert args[3:5] == [
        "--manifest",
        str(manifest.resolve()),
    ]

    assert str(spool_dir) not in calls[0]
