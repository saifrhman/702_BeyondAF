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


# ===========================================================================
# Paginated table access
# ===========================================================================
#
# The Artefact Viewer must stay usable against tables with millions of rows.
# Nothing here ever materialises a whole file: Parquet is read row-group by
# row-group and stops as soon as the requested page is filled, and CSV is
# streamed line by line.

#: Hard cap on one page, whatever the caller asks for.
MAX_PAGE_SIZE = 200

DEFAULT_PAGE_SIZE = 50

#: Rows scanned before a server-side filter gives up, so a filter over a
#: multi-million-row table can never hang the UI.
MAX_FILTER_SCAN_ROWS = 200_000


def _clamp_page(page: int, page_size: int) -> tuple[int, int]:
    size = max(1, min(int(page_size or DEFAULT_PAGE_SIZE), MAX_PAGE_SIZE))
    number = max(1, int(page or 1))

    return number, size


def parquet_schema(path: str | Path) -> dict[str, Any]:
    """Return schema and metadata without reading any row data."""

    import pyarrow.parquet as pq

    target = Path(path)

    try:
        parquet_file = pq.ParquetFile(str(target))
    except Exception as exc:  # pragma: no cover - corrupt file
        raise ArtefactError(f"Could not open Parquet file: {exc}") from exc

    schema = parquet_file.schema_arrow

    metadata: dict[str, str] = {}

    if schema.metadata:
        for key, value in schema.metadata.items():
            try:
                metadata[key.decode("utf-8")] = value.decode("utf-8")
            except UnicodeDecodeError:
                metadata[repr(key)] = repr(value)

    return {
        "columns": [
            {"name": field.name, "type": str(field.type)}
            for field in schema
        ],
        "column_count": len(schema),
        "row_count": parquet_file.metadata.num_rows,
        "row_group_count": parquet_file.num_row_groups,
        "schema_metadata": metadata,
        "created_by": parquet_file.metadata.created_by,
    }


def _matches(row: dict[str, Any], needle: str, columns: list[str]) -> bool:
    lowered = needle.lower()

    for column in columns:
        value = row.get(column)

        if value is None:
            continue

        if lowered in str(value).lower():
            return True

    return False


def parquet_page(
    path: str | Path,
    *,
    page: int = 1,
    page_size: int = DEFAULT_PAGE_SIZE,
    columns: list[str] | None = None,
    search: str | None = None,
    sort_by: str | None = None,
    descending: bool = False,
) -> dict[str, Any]:
    """Return one bounded page of a Parquet table.

    Reads only what the page needs. Sorting is applied to the *scanned window*
    rather than the whole file, and says so, because globally sorting a
    multi-million-row table on demand is not something a UI request should do.
    """

    import pyarrow.parquet as pq

    target = Path(path)
    number, size = _clamp_page(page, page_size)

    info = parquet_schema(target)
    available = [column["name"] for column in info["columns"]]

    selected = [c for c in (columns or available) if c in available]

    if not selected:
        selected = available

    parquet_file = pq.ParquetFile(str(target))

    start = (number - 1) * size
    end = start + size

    rows: list[dict[str, Any]] = []
    scanned = 0
    matched = 0
    truncated_scan = False

    filtering = bool(search) or bool(sort_by)

    # Sorting needs a window larger than one page to be meaningful.
    window_end = (
        min(MAX_FILTER_SCAN_ROWS, max(end, 5 * size)) if sort_by else end
    )

    for batch in parquet_file.iter_batches(
        batch_size=min(4096, max(size, 512)),
        columns=selected,
    ):
        if not filtering and scanned >= end:
            break

        for row in batch.to_pylist():
            scanned += 1

            if search and not _matches(row, search, selected):
                continue

            matched += 1

            if sort_by:
                if len(rows) < window_end:
                    rows.append(_json_safe(row))
            elif start < matched <= end:
                rows.append(_json_safe(row))

            if not sort_by and matched >= end and not search:
                break

        if scanned >= MAX_FILTER_SCAN_ROWS and filtering:
            truncated_scan = True
            break

        if not filtering and scanned >= end:
            break

    if sort_by and sort_by in selected:
        def _key(item: dict[str, Any]) -> tuple[int, Any]:
            value = item.get(sort_by)

            # None sorts last in both directions, and mixed types never raise.
            return (1, "") if value is None else (0, _sortable(value))

        rows.sort(key=_key, reverse=descending)
        rows = rows[start:end]

    total = info["row_count"]

    return {
        "columns": [
            column for column in info["columns"] if column["name"] in selected
        ],
        "rows": rows,
        "page": number,
        "page_size": size,
        "row_count": total,
        "returned": len(rows),
        "matched_rows": matched if search else total,
        "page_count": max(1, (total + size - 1) // size),
        "search": search,
        "sort_by": sort_by,
        "descending": descending,
        "scanned_rows": scanned,
        "scan_truncated": truncated_scan,
        "sort_scope": (
            "scanned window" if sort_by else None
        ),
        "note": (
            "Sorting and search are applied to a bounded scan of the file, "
            "not to the whole table."
            if filtering
            else None
        ),
    }


def _sortable(value: Any) -> Any:
    if isinstance(value, (int, float, bool)):
        return value

    return str(value)


def csv_page(
    path: str | Path,
    *,
    page: int = 1,
    page_size: int = DEFAULT_PAGE_SIZE,
    search: str | None = None,
) -> dict[str, Any]:
    """Return one bounded page of a CSV/TSV file, streamed line by line."""

    target = Path(path)
    number, size = _clamp_page(page, page_size)

    delimiter = "\t" if target.suffix.lower() == ".tsv" else ","

    start = (number - 1) * size
    end = start + size

    rows: list[list[str]] = []
    columns: list[str] = []
    matched = 0
    scanned = 0

    with target.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)

        try:
            columns = next(reader)
        except StopIteration:
            columns = []

        for row in reader:
            scanned += 1

            if search:
                joined = delimiter.join(row).lower()

                if search.lower() not in joined:
                    continue

            matched += 1

            if start < matched <= end:
                rows.append(row)

            if matched >= end and not search:
                break

            if scanned >= MAX_FILTER_SCAN_ROWS:
                break

    return {
        "columns": columns,
        "rows": rows,
        "page": number,
        "page_size": size,
        "returned": len(rows),
        "matched_rows": matched,
        "scanned_rows": scanned,
        "search": search,
    }


def table_page(
    path: str | Path,
    *,
    page: int = 1,
    page_size: int = DEFAULT_PAGE_SIZE,
    search: str | None = None,
    sort_by: str | None = None,
    descending: bool = False,
) -> dict[str, Any]:
    """Dispatch to the right paginated reader for a tabular artefact."""

    kind = preview_kind(path)

    if kind == "parquet":
        payload = parquet_page(
            path,
            page=page,
            page_size=page_size,
            search=search,
            sort_by=sort_by,
            descending=descending,
        )
        payload["kind"] = "parquet"

        return payload

    if kind == "table":
        payload = csv_page(path, page=page, page_size=page_size, search=search)
        payload["kind"] = "table"

        return payload

    raise ArtefactError(f"Not a paginated tabular artefact: {path}")


def rows_to_csv(columns: list[str], rows: list[Any]) -> str:
    """Render the currently displayed rows as CSV.

    This is a **convenience export of the current view**, not the authoritative
    scientific artefact. The original Parquet remains authoritative.
    """

    buffer = io.StringIO()
    writer = csv.writer(buffer)

    writer.writerow(columns)

    for row in rows:
        if isinstance(row, dict):
            writer.writerow([row.get(column, "") for column in columns])
        else:
            writer.writerow(row)

    return buffer.getvalue()
