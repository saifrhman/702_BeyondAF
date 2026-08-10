"""Tests for reproducible execution provenance."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from pdbclean.provenance import (
    ProvenanceError,
    resolve_clean_git_commit,
)


def _git(
    repository: Path,
    *args: str,
) -> str:
    result = subprocess.run(
        ["git", "-C", str(repository), *args],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def _clean_repository(tmp_path: Path) -> Path:
    repository = tmp_path / "repo"
    repository.mkdir()

    _git(repository, "init")
    _git(repository, "config", "user.name", "PDBClean Test")
    _git(
        repository,
        "config",
        "user.email",
        "pdbclean-test@example.invalid",
    )

    tracked = repository / "tracked.txt"
    tracked.write_text("initial\n")

    _git(repository, "add", "tracked.txt")
    _git(repository, "commit", "-m", "Initial commit")

    return repository


def test_resolve_clean_git_commit_returns_full_head(
    tmp_path: Path,
) -> None:
    repository = _clean_repository(tmp_path)

    expected = _git(repository, "rev-parse", "HEAD")
    observed = resolve_clean_git_commit(repository)

    assert observed == expected.lower()
    assert len(observed) == 40


def test_resolve_clean_git_commit_rejects_modified_file(
    tmp_path: Path,
) -> None:
    repository = _clean_repository(tmp_path)

    (repository / "tracked.txt").write_text("modified\n")

    with pytest.raises(
        ProvenanceError,
        match="worktree is not clean",
    ):
        resolve_clean_git_commit(repository)


def test_resolve_clean_git_commit_rejects_staged_change(
    tmp_path: Path,
) -> None:
    repository = _clean_repository(tmp_path)

    (repository / "tracked.txt").write_text("staged\n")
    _git(repository, "add", "tracked.txt")

    with pytest.raises(
        ProvenanceError,
        match="worktree is not clean",
    ):
        resolve_clean_git_commit(repository)


def test_resolve_clean_git_commit_rejects_untracked_file(
    tmp_path: Path,
) -> None:
    repository = _clean_repository(tmp_path)

    (repository / "untracked.txt").write_text("untracked\n")

    with pytest.raises(
        ProvenanceError,
        match="worktree is not clean",
    ):
        resolve_clean_git_commit(repository)


def test_resolve_clean_git_commit_rejects_non_repository(
    tmp_path: Path,
) -> None:
    with pytest.raises(
        ProvenanceError,
        match="Could not resolve Git provenance",
    ):
        resolve_clean_git_commit(tmp_path)
