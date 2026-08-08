#!/usr/bin/env python3
"""Generate an immutable manifest for a fixed PDB snapshot."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.manifest import file_sha256, validate_manifest_table
from pdbclean.schemas import SOURCE_MANIFEST_SCHEMA
from pdbclean.snapshot import iter_snapshot_objects, resolve_snapshot


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate and validate a PDB snapshot mmCIF manifest."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to the versioned PDBClean YAML configuration.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory in which manifest outputs will be written.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    loaded = load_config(args.config)
    config = loaded.data

    snapshot_config = config["snapshot"]
    execution_config = config["execution"]

    resolved = resolve_snapshot(
        snapshot_config,
        timeout_seconds=execution_config[
            "connection_timeout_seconds"
        ],
    )

    snapshot = resolved.snapshot_id
    source_prefix = resolved.source_prefix
    bucket_url = snapshot_config["bucket_url"]

    generated_at = datetime.now(timezone.utc)

    rows: list[dict[str, object]] = []

    for obj in iter_snapshot_objects(
        bucket_url=bucket_url,
        snapshot=snapshot,
        source_prefix=source_prefix,
        page_size=1000,
        timeout_seconds=execution_config[
            "connection_timeout_seconds"
        ],
    ):
        rows.append(
            {
                "snapshot": obj.snapshot,
                "source_layout": resolved.layout,
                "pdb_id": obj.pdb_id,
                "s3_key": obj.s3_key,
                "size_bytes": obj.size_bytes,
                "etag": obj.etag,
                "last_modified_utc": obj.last_modified_utc,
                "manifest_generated_at_utc": generated_at,
            }
        )

    rows.sort(key=lambda row: str(row["pdb_id"]))

    table = pa.Table.from_pylist(
        rows,
        schema=SOURCE_MANIFEST_SCHEMA,
    )

    summary = validate_manifest_table(
        table,
        expected_snapshot=snapshot,
        expected_count=snapshot_config.get(
            "expected_mmcif_count"
        ),
        expected_total_bytes=snapshot_config.get(
            "expected_total_bytes"
        ),
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_path = output_dir / "source_manifest.csv"
    parquet_path = output_dir / "source_manifest.parquet"
    summary_path = output_dir / "source_manifest_summary.json"

    csv_tmp = output_dir / "source_manifest.csv.tmp"
    parquet_tmp = output_dir / "source_manifest.parquet.tmp"
    summary_tmp = output_dir / "source_manifest_summary.json.tmp"

    pacsv.write_csv(table, csv_tmp)

    pq.write_table(
        table,
        parquet_tmp,
        compression="zstd",
        version="2.6",
    )

    csv_tmp.replace(csv_path)
    parquet_tmp.replace(parquet_path)

    result = {
        "snapshot_selection_mode": resolved.selection_mode,
        "resolved_snapshot": snapshot,
        "source_prefix": source_prefix,
        "sample_coordinate_mmcif": resolved.sample_mmcif_key,
        "row_count": summary.row_count,
        "total_bytes": summary.total_bytes,
        "unique_pdb_ids": summary.unique_pdb_ids,
        "unique_s3_keys": summary.unique_s3_keys,
        "config_sha256": loaded.sha256,
        "csv_sha256": file_sha256(csv_path),
        "parquet_sha256": file_sha256(parquet_path),
        "manifest_generated_at_utc": generated_at.isoformat(),
    }

    summary_tmp.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    summary_tmp.replace(summary_path)

    print(f"Snapshot: {snapshot}")
    print(f"mmCIF files: {summary.row_count:,}")
    print(f"Total bytes: {summary.total_bytes:,}")
    print(
        "Total GiB: "
        f"{summary.total_bytes / (1024 ** 3):.2f}"
    )
    print(f"CSV: {csv_path}")
    print(f"Parquet: {parquet_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
