#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import socket
import sys
import time

import pyarrow.compute as pc
import pyarrow.parquet as pq

from pdbclean.snapshot import download_verified_s3_object_bytes


BUCKET_URL = "https://pdbsnapshots.s3.us-west-2.amazonaws.com"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()

    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)

    return h.hexdigest()


def atomic_write_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    tmp = path.with_name(
        f".{path.name}.tmp.{os.getpid()}"
    )

    try:
        with tmp.open("wb") as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())

        os.replace(tmp, path)

    finally:
        if tmp.exists():
            tmp.unlink()


def atomic_write_json(path: Path, value: object) -> None:
    encoded = (
        json.dumps(value, indent=2, sort_keys=True)
        + "\n"
    ).encode("utf-8")

    atomic_write_bytes(path, encoded)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()

    p.add_argument(
        "--worklist",
        required=True,
    )

    p.add_argument(
        "--batch-id",
        required=True,
        type=int,
    )

    p.add_argument(
        "--status-dir",
        required=True,
    )

    p.add_argument(
        "--expected-worklist-sha256",
        required=True,
    )

    return p.parse_args()


def validate_existing_pass(
    *,
    report_path: Path,
    batch_id: int,
    worklist_sha256: str,
) -> bool:

    if not report_path.is_file():
        return False

    try:
        report = json.loads(
            report_path.read_text()
        )
    except Exception:
        return False

    if report.get("status") != "PASS":
        return False

    if report.get("batch_id") != batch_id:
        return False

    if (
        report.get("worklist_sha256")
        != worklist_sha256
    ):
        return False

    outputs = report.get("outputs")

    if not isinstance(outputs, list):
        return False

    for row in outputs:
        path = Path(row["output_cif_path"])

        if not path.is_file():
            return False

        if path.stat().st_size != row["cif_size_bytes"]:
            return False

    return True


def main() -> None:
    args = parse_args()

    worklist = Path(args.worklist).resolve()
    status_dir = Path(args.status_dir).resolve()

    status_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not worklist.is_file():
        raise RuntimeError(
            f"Worklist does not exist: {worklist}"
        )

    actual_worklist_sha = sha256_file(worklist)

    if (
        actual_worklist_sha
        != args.expected_worklist_sha256
    ):
        raise RuntimeError(
            "Worklist SHA256 mismatch: "
            f"expected={args.expected_worklist_sha256} "
            f"actual={actual_worklist_sha}"
        )

    batch_id = args.batch_id

    if batch_id < 0:
        raise RuntimeError(
            f"Invalid batch ID: {batch_id}"
        )

    report_path = (
        status_dir
        / f"batch_{batch_id:04d}.json"
    )

    if validate_existing_pass(
        report_path=report_path,
        batch_id=batch_id,
        worklist_sha256=actual_worklist_sha,
    ):
        print(
            f"Batch {batch_id}: existing PASS "
            "report and output files verified; skipping."
        )
        return

    table = pq.read_table(worklist)

    mask = pc.equal(
        table["batch_id"],
        batch_id,
    )

    batch = table.filter(mask)

    if batch.num_rows == 0:
        raise RuntimeError(
            f"No worklist rows for batch {batch_id}"
        )

    rows = batch.to_pylist()

    print(
        "===== SOURCE MATERIALISATION BATCH ====="
    )
    print("batch_id:", batch_id)
    print("rows:", len(rows))
    print("host:", socket.gethostname())
    print("worklist:", worklist)
    print(
        "worklist sha256:",
        actual_worklist_sha,
    )
    print("bucket:", BUCKET_URL)
    print()

    started = time.time()
    outputs = []

    for ordinal, row in enumerate(rows, start=1):
        pdb_id = row["pdb_id"]
        s3_key = row["s3_key"]
        expected_size = int(
            row["size_bytes"]
        )
        expected_etag = row["etag"]

        output_path = Path(
            row["output_cif_path"]
        )

        print(
            f"[{ordinal:03d}/{len(rows):03d}] "
            f"{pdb_id} "
            f"compressed={expected_size}"
        )

        compressed = (
            download_verified_s3_object_bytes(
                bucket_url=BUCKET_URL,
                s3_key=s3_key,
                expected_size_bytes=expected_size,
                expected_etag=expected_etag,
                timeout_seconds=60,
            )
        )

        if len(compressed) != expected_size:
            raise RuntimeError(
                f"{pdb_id}: compressed size changed "
                "after verified download"
            )

        try:
            cif_bytes = gzip.decompress(
                compressed
            )
        except Exception as exc:
            raise RuntimeError(
                f"{pdb_id}: gzip decompression failed"
            ) from exc

        if not cif_bytes:
            raise RuntimeError(
                f"{pdb_id}: decompressed mmCIF is empty"
            )

        stripped = cif_bytes.lstrip()

        if not stripped.startswith(b"data_"):
            raise RuntimeError(
                f"{pdb_id}: decompressed object does "
                "not begin with an mmCIF data block"
            )

        atomic_write_bytes(
            output_path,
            cif_bytes,
        )

        if (
            output_path.stat().st_size
            != len(cif_bytes)
        ):
            raise RuntimeError(
                f"{pdb_id}: output size mismatch "
                "after atomic write"
            )

        cif_sha = sha256_file(output_path)

        outputs.append(
            {
                "work_index":
                    int(row["work_index"]),
                "pdb_id":
                    pdb_id,
                "s3_key":
                    s3_key,
                "source_etag":
                    expected_etag,
                "compressed_size_bytes":
                    expected_size,
                "output_cif_path":
                    str(output_path),
                "cif_size_bytes":
                    len(cif_bytes),
                "cif_sha256":
                    cif_sha,
            }
        )

    elapsed = time.time() - started

    report = {
        "status": "PASS",
        "snapshot": "20260101",
        "batch_id": batch_id,
        "row_count": len(rows),
        "host": socket.gethostname(),
        "worklist": str(worklist),
        "worklist_sha256":
            actual_worklist_sha,
        "bucket_url": BUCKET_URL,
        "elapsed_seconds": elapsed,
        "outputs": outputs,
    }

    atomic_write_json(
        report_path,
        report,
    )

    print()
    print(
        "===== SOURCE MATERIALISATION "
        "BATCH COMPLETE ====="
    )
    print("batch_id:", batch_id)
    print("materialised:", len(outputs))
    print(
        "decompressed bytes:",
        sum(
            r["cif_size_bytes"]
            for r in outputs
        ),
    )
    print(
        "elapsed seconds:",
        f"{elapsed:.2f}",
    )
    print("report:", report_path)
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
