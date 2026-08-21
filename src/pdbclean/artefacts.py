"""Bounded, read-only previews of pipeline artefacts.

Historical inspection must never load a multi-gigabyte dataset into a browser,
and must never modify anything it looks at.  Every function here opens files
read-only, caps how much it will materialise, and reports what it truncated.

Supported:

``.json``      parsed and returned whole when small, else metadata only
``.csv``/``.tsv``  header plus a bounded number of rows
``.parquet``   schema, row count, key/value metadata, bounded row preview
``.txt``/``.log``/``.md``/``.yaml``  bounded text
anything else  metadata only (path, size, hash) -- never the bytes
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
from typing import Any


#: Never preview more rows than this in one request, whatever the caller asks.
MAX_PREVIEW_ROWS = 200

#: Never read more than this many bytes of a text-like artefact.
MAX_TEXT_BYTES = 256 * 1024

#: Parse a JSON artefact whole only below this size.
MAX_JSON_BYTES = 2 * 1024 * 1024

#: Above this, report size only rather than hashing a huge file inline.
MAX_HASH_BYTES = 256 * 1024 * 1024

TEXT_SUFFIXES = frozenset({".txt", ".log", ".md", ".yaml", ".yml", ".out", ".err"})
TABLE_SUFFIXES = frozenset({".csv", ".tsv"})


class ArtefactError(RuntimeError):
    """Raised when an artefact cannot be previewed."""


def _sha256(path: Path) -> str | None:
    if path.stat().st_size > MAX_HASH_BYTES:
        return None

    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def describe(path: str | Path, *, compute_hash: bool = True) -> dict[str, Any]:
    """Return metadata for one artefact without reading its contents."""

    target = Path(path)

    if not target.exists():
        return {
            "path": str(target),
            "name": target.name,
            "exists": False,
            "kind": "missing",
        }

    if target.is_dir():
        return {
            "path": str(target),
            "name": target.name,
            "exists": True,
            "kind": "directory",
        }

    stat = target.stat()

    return {
        "path": str(target),
        "name": target.name,
        "exists": True,
        "kind": "file",
        "suffix": target.suffix.lower(),
        "bytes": stat.st_size,
        "modified": stat.st_mtime,
        "sha256": _sha256(target) if compute_hash else None,
        "previewable": preview_kind(target) != "metadata",
    }


def preview_kind(path: str | Path) -> str:
    """Return how an artefact would be previewed."""

    suffix = Path(path).suffix.lower()

    if suffix == ".json":
        return "json"

    if suffix == ".parquet":
        return "parquet"

    if suffix in TABLE_SUFFIXES:
        return "table"

    if suffix in TEXT_SUFFIXES or Path(path).name == "_SUCCESS":
        return "text"

    return "metadata"


def _preview_json(path: Path) -> dict[str, Any]:
    size = path.stat().st_size

    if size > MAX_JSON_BYTES:
        return {
            "kind": "json",
            "truncated": True,
            "reason": (
                f"JSON artefact is {size} bytes, above the "
                f"{MAX_JSON_BYTES}-byte inline limit"
            ),
        }

    try:
        content = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ArtefactError(f"Could not parse JSON: {exc}") from exc

    return {"kind": "json", "truncated": False, "content": content}


def _preview_text(path: Path) -> dict[str, Any]:
    size = path.stat().st_size

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        text = handle.read(MAX_TEXT_BYTES)

    return {
        "kind": "text",
        "truncated": size > MAX_TEXT_BYTES,
        "bytes": size,
        "content": text,
    }


def _preview_table(path: Path, *, limit: int) -> dict[str, Any]:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","

    rows: list[list[str]] = []
    columns: list[str] = []
    truncated = False

    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)

        try:
            columns = next(reader)
        except StopIteration:
            columns = []

        for index, row in enumerate(reader):
            if index >= limit:
                truncated = True
                break

            rows.append(row)

    return {
        "kind": "table",
        "columns": columns,
        "rows": rows,
        "row_preview_count": len(rows),
        "truncated": truncated,
    }


def _preview_parquet(path: Path, *, limit: int) -> dict[str, Any]:
    import pyarrow.parquet as pq

    try:
        parquet_file = pq.ParquetFile(str(path))
    except Exception as exc:  # pragma: no cover - corrupt file
        raise ArtefactError(f"Could not open Parquet file: {exc}") from exc

    schema = parquet_file.schema_arrow
    total_rows = parquet_file.metadata.num_rows

    metadata: dict[str, str] = {}

    if schema.metadata:
        for key, value in schema.metadata.items():
            try:
                metadata[key.decode("utf-8")] = value.decode("utf-8")
            except UnicodeDecodeError:
                metadata[repr(key)] = repr(value)

    rows: list[dict[str, Any]] = []

    if limit > 0 and total_rows > 0:
        for batch in parquet_file.iter_batches(batch_size=min(limit, 1024)):
            for row in batch.to_pylist():
                rows.append(_json_safe(row))

                if len(rows) >= limit:
                    break

            if len(rows) >= limit:
                break

    return {
        "kind": "parquet",
        "columns": [
            {"name": field.name, "type": str(field.type)}
            for field in schema
        ],
        "row_count": total_rows,
        "row_group_count": parquet_file.num_row_groups,
        "schema_metadata": metadata,
        "rows": rows,
        "row_preview_count": len(rows),
        "truncated": len(rows) < total_rows,
    }


def _json_safe(value: Any) -> Any:
    """Make a Parquet row JSON-serialisable without changing its meaning."""

    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}

    if isinstance(value, (list, tuple)):
        # Bound nested arrays -- a BRI cell is (m, 9) and must not be inlined.
        rendered = [_json_safe(v) for v in value[:16]]

        if len(value) > 16:
            rendered.append(f"... {len(value) - 16} more")

        return rendered

    if isinstance(value, (bytes, bytearray)):
        return value.decode("utf-8", errors="replace")

    if isinstance(value, (str, int, float, bool)) or value is None:
        return value

    return str(value)


def preview(
    path: str | Path,
    *,
    limit: int = 50,
) -> dict[str, Any]:
    """Return a bounded, read-only preview of one artefact."""

    target = Path(path)

    if not target.is_file():
        raise ArtefactError(f"Not a readable artefact: {target}")

    bounded = max(0, min(int(limit), MAX_PREVIEW_ROWS))

    payload = describe(target)
    kind = preview_kind(target)

    if kind == "json":
        payload["preview"] = _preview_json(target)
    elif kind == "text":
        payload["preview"] = _preview_text(target)
    elif kind == "table":
        payload["preview"] = _preview_table(target, limit=bounded)
    elif kind == "parquet":
        payload["preview"] = _preview_parquet(target, limit=bounded)
    else:
        payload["preview"] = {
            "kind": "metadata",
            "truncated": True,
            "reason": (
                "No safe inline viewer for this file type; showing metadata "
                "only."
            ),
        }

    return payload


def list_directory(
    directory: str | Path,
    *,
    limit: int = 200,
    compute_hash: bool = False,
) -> list[dict[str, Any]]:
    """List a stage output directory, newest-looking artefacts first.

    Hashing is off by default: listing a stage must stay cheap enough to run
    on a login node.
    """

    root = Path(directory)

    if not root.is_dir():
        return []

    entries: list[dict[str, Any]] = []

    for child in sorted(root.rglob("*")):
        if len(entries) >= limit:
            break

        if child.is_dir():
            continue

        described = describe(child, compute_hash=compute_hash)
        described["relative_path"] = str(child.relative_to(root))
        entries.append(described)

    return entries
