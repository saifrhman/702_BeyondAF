"""Read-only query layer over the duplicate-detection results.

This module powers the Duplicate Explorer in the CLI and the UI.  It is a
*view*: it reads the Parquet artefacts the pipeline already produced and never
recomputes, re-derives or re-classifies anything.  The complete-BRI result
remains authoritative, and no new duplicate category is introduced here.

Sources, in order of preference:

``duplicate_classification/finalized/candidate_classifications.parquet``
    Stage-10 classifications, carrying the authoritative
    ``is_zero_duplicate`` / ``is_paper_near_duplicate`` /
    ``is_nonzero_near_duplicate`` flags alongside ``d_bri_mA``.

``full_bri_nn/finalized/candidate_near_duplicates.parquet``
    Stage-8 near-duplicate pairs, used when Stage 10 has not been published.

``<release>/audit/representative_mapping.parquet``
    Stage-14 retained/removed relationships, joined in when available.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq


PAIR_COLUMNS = [
    "query_snapshot",
    "query_pdb_id",
    "query_model_id",
    "query_label_chain_id",
    "subject_snapshot",
    "subject_pdb_id",
    "subject_model_id",
    "subject_label_chain_id",
    "retained_residue_count",
    "d_bri_mA",
    "d_bri",
]

CLASSIFICATION_COLUMNS = [
    "is_zero_duplicate",
    "is_paper_near_duplicate",
    "is_nonzero_near_duplicate",
]

#: Hard cap on a single page, so a UI request can never materialise the whole
#: pair table into memory.
MAX_PAGE_SIZE = 500


class DuplicateQueryError(RuntimeError):
    """Raised when the duplicate results cannot be queried."""


@dataclass(frozen=True)
class DuplicateSource:
    """Where the Duplicate Explorer reads from for one run."""

    protocol_root: Path
    release_root: Path | None = None

    @property
    def classification_paths(self) -> list[Path]:
        base = self.protocol_root / "duplicate_classification" / "finalized"

        return [
            base / "candidate_classifications.parquet",
            base / "m1_classifications.parquet",
        ]

    @property
    def near_duplicate_paths(self) -> list[Path]:
        base = self.protocol_root / "full_bri_nn" / "finalized"

        return [
            base / "candidate_near_duplicates.parquet",
            base / "m1_near_duplicates.parquet",
        ]

    @property
    def representative_mapping_path(self) -> Path | None:
        if self.release_root is None:
            return None

        candidate = (
            self.release_root / "audit" / "representative_mapping.parquet"
        )

        if candidate.is_file():
            return candidate

        fallback = (
            self.protocol_root
            / "stage14_representative_selection_v1"
            / "representative_mapping.parquet"
        )

        return fallback if fallback.is_file() else None

    def available(self) -> dict[str, Any]:
        classified = [p for p in self.classification_paths if p.is_file()]
        near = [p for p in self.near_duplicate_paths if p.is_file()]

        return {
            "classified": [str(p) for p in classified],
            "near_duplicates": [str(p) for p in near],
            "representative_mapping": (
                str(self.representative_mapping_path)
                if self.representative_mapping_path
                else None
            ),
            "mode": "classified" if classified else ("near" if near else None),
        }


@dataclass
class DuplicateFilters:
    """Filters exposed by the Duplicate Explorer.

    Every filter narrows the existing result set.  None of them changes how a
    pair was classified.
    """

    pdb_id: str | None = None
    chain: str | None = None
    exact_only: bool = False
    nonzero_near_only: bool = False
    min_length: int | None = None
    max_length: int | None = None
    min_distance_mA: int | None = None
    max_distance_mA: int | None = None
    relationship: str | None = None  # "removed" | "retained" | "unaffected"
    offset: int = 0
    limit: int = 100

    def normalised(self) -> "DuplicateFilters":
        limit = max(1, min(int(self.limit or 100), MAX_PAGE_SIZE))
        offset = max(0, int(self.offset or 0))

        if self.exact_only and self.nonzero_near_only:
            raise DuplicateQueryError(
                "exact_only and nonzero_near_only are mutually exclusive"
            )

        return DuplicateFilters(
            pdb_id=(self.pdb_id or "").strip().lower() or None,
            chain=(self.chain or "").strip() or None,
            exact_only=bool(self.exact_only),
            nonzero_near_only=bool(self.nonzero_near_only),
            min_length=self.min_length,
            max_length=self.max_length,
            min_distance_mA=self.min_distance_mA,
            max_distance_mA=self.max_distance_mA,
            relationship=(self.relationship or None),
            offset=offset,
            limit=limit,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "pdb_id": self.pdb_id,
            "chain": self.chain,
            "exact_only": self.exact_only,
            "nonzero_near_only": self.nonzero_near_only,
            "min_length": self.min_length,
            "max_length": self.max_length,
            "min_distance_mA": self.min_distance_mA,
            "max_distance_mA": self.max_distance_mA,
            "relationship": self.relationship,
            "offset": self.offset,
            "limit": self.limit,
        }


def _chain_key(snapshot: Any, pdb_id: Any, model_id: Any, chain: Any) -> tuple:
    return (
        str(snapshot),
        str(pdb_id).lower(),
        int(model_id),
        str(chain),
    )


@lru_cache(maxsize=8)
def _load_representative_mapping(path: str) -> dict[tuple, dict[str, Any]]:
    """Load the Stage-14 mapping keyed by canonical chain identity."""

    table = pq.read_table(
        path,
        columns=[
            "snapshot",
            "pdb_id",
            "model_id",
            "label_chain_id",
            "action",
            "component_id",
            "component_is_clique",
            "representative_snapshot",
            "representative_pdb_id",
            "representative_model_id",
            "representative_label_chain_id",
            "direct_d_bri_mA",
            "policy_version",
        ],
    )

    mapping: dict[tuple, dict[str, Any]] = {}

    for row in table.to_pylist():
        mapping[
            _chain_key(
                row["snapshot"],
                row["pdb_id"],
                row["model_id"],
                row["label_chain_id"],
            )
        ] = row

    return mapping


class DuplicateExplorer:
    """Paginated, filterable view over the detected duplicate pairs."""

    def __init__(self, source: DuplicateSource) -> None:
        self.source = source
        self._availability = source.available()

        if self._availability["mode"] is None:
            raise DuplicateQueryError(
                "No duplicate results are available for this configuration. "
                f"Looked in {source.protocol_root}"
            )

        self._classified = self._availability["mode"] == "classified"

        self._paths = [
            Path(p)
            for p in (
                self._availability["classified"]
                if self._classified
                else self._availability["near_duplicates"]
            )
        ]

    # -- metadata -------------------------------------------------------

    @property
    def has_classification_flags(self) -> bool:
        return self._classified

    @property
    def has_representative_mapping(self) -> bool:
        return self._availability["representative_mapping"] is not None

    def summary(self) -> dict[str, Any]:
        """Return the pair-population summary from stage provenance.

        Counts come from the stage summaries the pipeline already wrote, so the
        Explorer never presents a number it derived itself.
        """

        payload: dict[str, Any] = {
            "sources": self._availability,
            "classified": self._classified,
        }

        stage10 = (
            self.source.protocol_root
            / "duplicate_classification"
            / "global_summary.json"
        )

        if stage10.is_file():
            summary = json.loads(stage10.read_text(encoding="utf-8"))

            payload["total_tested_pairs"] = summary.get("input_pair_count")
            payload["near_duplicate_pairs"] = summary.get(
                "paper_near_duplicate_pair_count"
            )
            payload["exact_duplicate_pairs"] = summary.get(
                "zero_duplicate_pair_count"
            )
            payload["nonzero_near_duplicate_pairs"] = summary.get(
                "nonzero_near_duplicate_pair_count"
            )
            payload["non_near_duplicate_pairs"] = summary.get(
                "not_near_duplicate_pair_count"
            )
            payload["threshold_mA"] = summary.get(
                "paper_near_duplicate_threshold_mA"
            )
            payload["threshold_angstrom"] = summary.get(
                "paper_near_duplicate_threshold_angstrom"
            )

        return payload

    # -- querying -------------------------------------------------------

    def _expression(self, filters: DuplicateFilters):
        clauses = []

        if self._classified:
            # Only near duplicates are shown; the classification flag is the
            # authoritative complete-BRI decision.
            clauses.append(pc.field("is_paper_near_duplicate"))

            if filters.exact_only:
                clauses.append(pc.field("is_zero_duplicate"))

            if filters.nonzero_near_only:
                clauses.append(pc.field("is_nonzero_near_duplicate"))
        else:
            if filters.exact_only:
                clauses.append(pc.field("d_bri_mA") == 0)

            if filters.nonzero_near_only:
                clauses.append(pc.field("d_bri_mA") > 0)

        if filters.min_length is not None:
            clauses.append(
                pc.field("retained_residue_count") >= int(filters.min_length)
            )

        if filters.max_length is not None:
            clauses.append(
                pc.field("retained_residue_count") <= int(filters.max_length)
            )

        if filters.min_distance_mA is not None:
            clauses.append(pc.field("d_bri_mA") >= int(filters.min_distance_mA))

        if filters.max_distance_mA is not None:
            clauses.append(pc.field("d_bri_mA") <= int(filters.max_distance_mA))

        if filters.pdb_id:
            lowered = filters.pdb_id
            clauses.append(
                (pc.field("query_pdb_id") == lowered)
                | (pc.field("subject_pdb_id") == lowered)
            )

        if filters.chain:
            chain = filters.chain
            clauses.append(
                (pc.field("query_label_chain_id") == chain)
                | (pc.field("subject_label_chain_id") == chain)
            )

        if not clauses:
            return None

        expression = clauses[0]

        for clause in clauses[1:]:
            expression = expression & clause

        return expression

    def _columns(self) -> list[str]:
        if self._classified:
            return PAIR_COLUMNS + CLASSIFICATION_COLUMNS

        return PAIR_COLUMNS

    def query(self, filters: DuplicateFilters) -> dict[str, Any]:
        """Return one page of duplicate pairs plus the matched-row count."""

        active = filters.normalised()
        expression = self._expression(active)
        columns = self._columns()

        available = [p for p in self._paths if p.is_file()]

        if not available:
            raise DuplicateQueryError("No duplicate parquet files present")

        dataset = ds.dataset(
            [str(p) for p in available],
            format="parquet",
        )

        # Relationship filtering needs the Stage-14 mapping, which is not part
        # of the pair tables, so it is applied after the columnar scan.
        needs_post_filter = bool(active.relationship)

        scanner = dataset.scanner(
            columns=[c for c in columns if c in dataset.schema.names],
            filter=expression,
            batch_size=65536,
        )

        mapping: dict[tuple, dict[str, Any]] = {}

        mapping_path = self._availability["representative_mapping"]

        if mapping_path:
            mapping = _load_representative_mapping(mapping_path)

        rows: list[dict[str, Any]] = []

        start = active.offset
        end = active.offset + active.limit

        if needs_post_filter:
            # The relationship filter needs the Stage-14 mapping, which is not
            # a column of the pair tables, so the whole filtered result has to
            # be walked to know how many rows match.
            matched = 0

            for batch in scanner.to_batches():
                for row in batch.to_pylist():
                    enriched = self._enrich(row, mapping)

                    if enriched["relationship"] != active.relationship:
                        continue

                    if start <= matched < end:
                        rows.append(enriched)

                    matched += 1
        else:
            # Count in Arrow rather than in Python, then materialise only the
            # requested page.
            matched = dataset.count_rows(filter=expression)

            seen = 0

            for batch in scanner.to_batches():
                batch_rows = batch.num_rows

                if seen + batch_rows <= start:
                    seen += batch_rows
                    continue

                for row in batch.to_pylist():
                    if seen >= end:
                        break

                    if seen >= start:
                        rows.append(self._enrich(row, mapping))

                    seen += 1

                if seen >= end:
                    break

        return {
            "rows": rows,
            "matched": matched,
            "offset": active.offset,
            "limit": active.limit,
            "filters": active.to_dict(),
            "has_classification_flags": self._classified,
            "has_representative_mapping": bool(mapping),
        }

    def _enrich(
        self,
        row: dict[str, Any],
        mapping: dict[tuple, dict[str, Any]],
    ) -> dict[str, Any]:
        """Add display fields and the Stage-14 relationship, if published."""

        distance_mA = int(row["d_bri_mA"])

        query_key = _chain_key(
            row["query_snapshot"],
            row["query_pdb_id"],
            row["query_model_id"],
            row["query_label_chain_id"],
        )
        subject_key = _chain_key(
            row["subject_snapshot"],
            row["subject_pdb_id"],
            row["subject_model_id"],
            row["subject_label_chain_id"],
        )

        is_exact = (
            bool(row["is_zero_duplicate"])
            if "is_zero_duplicate" in row
            else distance_mA == 0
        )

        enriched: dict[str, Any] = {
            "pdb_id_a": str(row["query_pdb_id"]).lower(),
            "chain_a": row["query_label_chain_id"],
            "model_a": int(row["query_model_id"]),
            "snapshot_a": row["query_snapshot"],
            "pdb_id_b": str(row["subject_pdb_id"]).lower(),
            "chain_b": row["subject_label_chain_id"],
            "model_b": int(row["subject_model_id"]),
            "snapshot_b": row["subject_snapshot"],
            "chain_length": int(row["retained_residue_count"]),
            "d_bri_mA": distance_mA,
            "d_bri_angstrom": (
                float(row["d_bri"])
                if row.get("d_bri") is not None
                else distance_mA / 1000.0
            ),
            "classification": (
                "exact_duplicate" if is_exact else "nonzero_near_duplicate"
            ),
            "relationship": "unknown",
            "representative": None,
            "removed_chain": None,
            "direct_d_bri_mA": None,
            "component_id": None,
            "component_is_clique": None,
            "policy_version": None,
        }

        if not mapping:
            return enriched

        decisions = []

        for label, key in (("a", query_key), ("b", subject_key)):
            entry = mapping.get(key)

            if entry is not None:
                decisions.append((label, entry))

        if not decisions:
            enriched["relationship"] = "unaffected"
            return enriched

        removed = [
            (label, entry)
            for label, entry in decisions
            if entry["action"] == "remove"
        ]

        first = decisions[0][1]
        enriched["component_id"] = first.get("component_id")
        enriched["component_is_clique"] = first.get("component_is_clique")
        enriched["policy_version"] = first.get("policy_version")

        if removed:
            label, entry = removed[0]

            enriched["relationship"] = "removed"
            enriched["removed_chain"] = (
                f"{entry['pdb_id']}:{entry['label_chain_id']}"
            )
            enriched["representative"] = (
                f"{entry['representative_pdb_id']}"
                f":{entry['representative_label_chain_id']}"
            )
            enriched["direct_d_bri_mA"] = entry.get("direct_d_bri_mA")
        else:
            enriched["relationship"] = "retained"

        return enriched


def pair_detail(
    explorer: DuplicateExplorer,
    *,
    pdb_id_a: str,
    chain_a: str,
    pdb_id_b: str,
    chain_b: str,
) -> dict[str, Any] | None:
    """Return the single pair matching two chain identities, or None."""

    filters = DuplicateFilters(pdb_id=pdb_id_a, limit=MAX_PAGE_SIZE)

    result = explorer.query(filters)

    wanted = {
        (pdb_id_a.lower(), chain_a, pdb_id_b.lower(), chain_b),
        (pdb_id_b.lower(), chain_b, pdb_id_a.lower(), chain_a),
    }

    for row in result["rows"]:
        key = (
            row["pdb_id_a"],
            row["chain_a"],
            row["pdb_id_b"],
            row["chain_b"],
        )

        if key in wanted:
            return row

    return None
