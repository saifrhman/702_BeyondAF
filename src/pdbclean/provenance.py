"""Execution provenance helpers for reproducible PDBClean releases."""

from __future__ import annotations

import subprocess
from pathlib import Path


class ProvenanceError(RuntimeError):
    """Raised when reproducible execution provenance cannot be established."""


def resolve_clean_git_commit(
    repository_root: str | Path,
) -> str:
    """Return full HEAD SHA only when the repository worktree is clean."""

    root = Path(repository_root).resolve()

    try:
        commit = subprocess.run(
            [
                "git",
                "-C",
                str(root),
                "rev-parse",
                "--verify",
                "HEAD",
            ],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()

        status = subprocess.run(
            [
                "git",
                "-C",
                str(root),
                "status",
                "--porcelain=v1",
                "--untracked-files=all",
            ],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ProvenanceError(
            f"Could not resolve Git provenance for {root}"
        ) from exc

    if len(commit) != 40 or any(
        character not in "0123456789abcdefABCDEF"
        for character in commit
    ):
        raise ProvenanceError(
            f"Git HEAD is not a full 40-character SHA: {commit!r}"
        )

    if status.strip():
        raise ProvenanceError(
            "Refusing reproducible production execution because "
            "the Git worktree is not clean"
        )

    return commit.lower()
