"""CLI for Stage-6 exact chain-length bucketing."""

from __future__ import annotations

import argparse
from pathlib import Path

from pdbclean.config import load_config
from pdbclean.length_buckets import (
    finalize_length_buckets,
    validate_stage5_brain_publication,
)


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--length-bucket-pipeline-git-commit",
        required=True,
    )

    args = parser.parse_args()

    repo = Path.cwd()

    config = load_config(
        args.config
    ).data

    protocol = (
        config["release"]["protocol_version"]
    )

    storage_root = Path(
        config["storage"]["output_root"]
    )

    if not storage_root.is_absolute():
        storage_root = repo / storage_root

    completed = sorted(
        path.parent
        for path in storage_root.glob(
            f"*/{protocol}/brain/_SUCCESS"
        )
    )

    if len(completed) != 1:
        raise RuntimeError(
            "Expected exactly one completed Stage-5 Brain "
            f"publication; found {len(completed)}"
        )

    upstream = validate_stage5_brain_publication(
        completed[0]
    )

    publication = finalize_length_buckets(
        upstream=upstream,
        output_root=(
            completed[0].parent
            / "length_buckets"
        ),
        length_bucket_pipeline_git_commit=(
            args.length_bucket_pipeline_git_commit
        ),
    )

    summary = publication.global_summary

    print(
        "Input chains:",
        summary["input_chain_count"],
    )
    print(
        "Distinct m buckets:",
        summary["distinct_length_bucket_count"],
    )
    print(
        "Minimum m:",
        summary["minimum_retained_residue_count"],
    )
    print(
        "Maximum m:",
        summary["maximum_retained_residue_count"],
    )
    print(
        "m=1 chains:",
        summary["m1_chain_count"],
    )
    print()
    print(
        "STAGE-6 EXACT LENGTH BUCKETING: PASS"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
