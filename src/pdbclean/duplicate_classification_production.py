"""Stage 10 production classification of exact Stage-8 BRI distances."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.duplicate_classification import (
    PAPER_NEAR_DUPLICATE_THRESHOLD_MA,
)
from pdbclean.full_bri_production import (
    CANDIDATE_COMPARISON_SCHEMA,
    M1_COMPARISON_SCHEMA,
)


class ClassificationProductionError(RuntimeError):
    """Raised when Stage-10 publication is invalid."""


CLASSIFICATION_FIELDS = [
    pa.field(
        "is_zero_duplicate",
        pa.bool_(),
        nullable=False,
    ),
    pa.field(
        "is_paper_near_duplicate",
        pa.bool_(),
        nullable=False,
    ),
    pa.field(
        "is_nonzero_near_duplicate",
        pa.bool_(),
        nullable=False,
    ),
]


CANDIDATE_CLASSIFIED_SCHEMA = pa.schema(
    list(CANDIDATE_COMPARISON_SCHEMA)
    + CLASSIFICATION_FIELDS,
    metadata={
        b"schema_name": (
            b"pdbclean_stage10_candidate_classifications"
        ),
        b"schema_version": b"1.0",
        b"classification_basis": (
            b"exact Stage-8 full-BRI integer milliangstrom distance"
        ),
        b"paper_near_duplicate_threshold_mA": b"10",
    },
)


M1_CLASSIFIED_SCHEMA = pa.schema(
    list(M1_COMPARISON_SCHEMA)
    + CLASSIFICATION_FIELDS,
    metadata={
        b"schema_name": (
            b"pdbclean_stage10_m1_classifications"
        ),
        b"schema_version": b"1.0",
        b"classification_basis": (
            b"exact Stage-8 full-BRI integer milliangstrom distance"
        ),
        b"paper_near_duplicate_threshold_mA": b"10",
    },
)


def _read_json(path: Path) -> dict:
    value = json.loads(
        path.read_text(encoding="utf-8")
    )

    if not isinstance(value, dict):
        raise ClassificationProductionError(
            f"JSON artifact is not an object: {path}"
        )

    return value


def _validate_commit(value: str) -> str:
    if (
        len(value) != 40
        or any(
            c not in "0123456789abcdef"
            for c in value
        )
    ):
        raise ClassificationProductionError(
            "Stage-10 commit must be a lowercase 40-character SHA"
        )

    return value


def _write_json_atomic(
    payload: dict,
    path: Path,
) -> Path:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = path.with_suffix(
        path.suffix + ".tmp"
    )

    temporary.write_text(
        json.dumps(
            payload,
            sort_keys=True,
            indent=2,
            ensure_ascii=True,
        )
        + "\n",
        encoding="utf-8",
    )

    temporary.replace(path)

    return path


def _classify_file(
    source: Path,
    destination: Path,
    *,
    input_schema: pa.Schema,
    output_schema: pa.Schema,
) -> dict[str, int]:
    """Classify every Stage-8 row without filtering."""

    source_pf = pq.ParquetFile(
        source
    )

    temporary = destination.with_suffix(
        destination.suffix + ".tmp"
    )

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if temporary.exists():
        temporary.unlink()

    writer = pq.ParquetWriter(
        temporary,
        output_schema,
        compression="zstd",
        version="2.6",
        use_dictionary=True,
    )

    counts = {
        "input": 0,
        "zero": 0,
        "paper_near": 0,
        "nonzero_near": 0,
        "not_near": 0,
    }

    try:
        for rg in range(
            source_pf.metadata.num_row_groups
        ):
            table = source_pf.read_row_group(
                rg
            ).combine_chunks()

            d = np.asarray(
                table[
                    "d_bri_mA"
                ].to_numpy(
                    zero_copy_only=False
                ),
                dtype=np.int64,
            )

            if np.any(d < 0):
                raise ClassificationProductionError(
                    "Negative Stage-8 d_bri_mA encountered"
                )

            zero = d == 0
            paper_near = (
                d
                <= PAPER_NEAR_DUPLICATE_THRESHOLD_MA
            )
            nonzero_near = (
                (d > 0)
                & paper_near
            )

            if np.any(
                zero & ~paper_near
            ):
                raise ClassificationProductionError(
                    "Zero duplicate is not a near duplicate"
                )

            if np.any(
                nonzero_near
                != (paper_near & ~zero)
            ):
                raise ClassificationProductionError(
                    "Nonzero-near classification inconsistency"
                )

            arrays = [
                table[field.name]
                for field in input_schema
            ]

            arrays.extend(
                [
                    pa.array(
                        zero,
                        type=pa.bool_(),
                    ),
                    pa.array(
                        paper_near,
                        type=pa.bool_(),
                    ),
                    pa.array(
                        nonzero_near,
                        type=pa.bool_(),
                    ),
                ]
            )

            output = pa.Table.from_arrays(
                arrays,
                schema=output_schema,
            )

            writer.write_table(
                output
            )

            n = table.num_rows

            counts["input"] += n
            counts["zero"] += int(
                np.count_nonzero(zero)
            )
            counts["paper_near"] += int(
                np.count_nonzero(
                    paper_near
                )
            )
            counts["nonzero_near"] += int(
                np.count_nonzero(
                    nonzero_near
                )
            )
            counts["not_near"] += int(
                np.count_nonzero(
                    ~paper_near
                )
            )

    finally:
        writer.close()

    if (
        counts["zero"]
        + counts["nonzero_near"]
        != counts["paper_near"]
    ):
        temporary.unlink(
            missing_ok=True
        )

        raise ClassificationProductionError(
            "Stage-10 near-duplicate accounting failed"
        )

    if (
        counts["paper_near"]
        + counts["not_near"]
        != counts["input"]
    ):
        temporary.unlink(
            missing_ok=True
        )

        raise ClassificationProductionError(
            "Stage-10 total classification accounting failed"
        )

    if (
        pq.read_metadata(
            temporary
        ).num_rows
        != counts["input"]
    ):
        temporary.unlink(
            missing_ok=True
        )

        raise ClassificationProductionError(
            "Stage-10 Parquet population mismatch"
        )

    temporary.replace(
        destination
    )

    return counts


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--pipeline-git-commit",
        required=True,
    )

    args = parser.parse_args()

    stage10_commit = _validate_commit(
        args.pipeline_git_commit
    )

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
        storage_root = (
            repo / storage_root
        )

    completed = sorted(
        path.parent
        for path in storage_root.glob(
            f"*/{protocol}/full_bri_compare/_SUCCESS"
        )
    )

    if len(completed) != 1:
        raise ClassificationProductionError(
            "Expected exactly one completed Stage-8 publication"
        )

    stage8_root = completed[0]

    success = _read_json(
        stage8_root / "_SUCCESS"
    )
    summary8 = _read_json(
        stage8_root / "global_summary.json"
    )

    if (
        success.get("success_schema_name")
        != "pdbclean_stage8_full_bri_success"
        or success.get("success_schema_version")
        != "1.0"
    ):
        raise ClassificationProductionError(
            "Unexpected Stage-8 _SUCCESS schema"
        )

    if (
        summary8.get("summary_schema_name")
        != "pdbclean_stage8_full_bri_global_summary"
        or summary8.get("summary_schema_version")
        != "1.0"
    ):
        raise ClassificationProductionError(
            "Unexpected Stage-8 summary schema"
        )

    if summary8.get(
        "classification_performed"
    ) is not False:
        raise ClassificationProductionError(
            "Stage-8 input unexpectedly reports classification"
        )

    candidate_input = (
        stage8_root
        / success["candidate_comparisons"]
    )
    m1_input = (
        stage8_root
        / success["m1_comparisons"]
    )

    if (
        pq.read_metadata(
            candidate_input
        ).num_rows
        != summary8[
            "brain_candidate_comparison_count"
        ]
    ):
        raise ClassificationProductionError(
            "Stage-8 candidate population mismatch"
        )

    if (
        pq.read_metadata(
            m1_input
        ).num_rows
        != summary8[
            "m1_direct_comparison_count"
        ]
    ):
        raise ClassificationProductionError(
            "Stage-8 m1 population mismatch"
        )

    output_root = (
        stage8_root.parent
        / "duplicate_classification"
    )

    finalized = (
        output_root / "finalized"
    )

    success_path = (
        output_root / "_SUCCESS"
    )

    if success_path.exists():
        success_path.unlink()

    candidate_counts = _classify_file(
        candidate_input,
        finalized
        / "candidate_classifications.parquet",
        input_schema=(
            CANDIDATE_COMPARISON_SCHEMA
        ),
        output_schema=(
            CANDIDATE_CLASSIFIED_SCHEMA
        ),
    )

    m1_counts = _classify_file(
        m1_input,
        finalized
        / "m1_classifications.parquet",
        input_schema=M1_COMPARISON_SCHEMA,
        output_schema=M1_CLASSIFIED_SCHEMA,
    )

    total = {
        key: (
            candidate_counts[key]
            + m1_counts[key]
        )
        for key in candidate_counts
    }

    if (
        total["input"]
        != summary8[
            "total_full_bri_comparison_count"
        ]
    ):
        raise ClassificationProductionError(
            "Stage-8 -> Stage-10 pair accounting failed"
        )

    provenance = {
        key: value
        for key, value in summary8.items()
        if key.endswith(
            "_pipeline_git_commit"
        )
    }

    global_summary = {
        "summary_schema_name": (
            "pdbclean_stage10_duplicate_classification_global_summary"
        ),
        "summary_schema_version": "1.0",
        "snapshot": summary8["snapshot"],
        "cleaning_protocol": (
            summary8["cleaning_protocol"]
        ),
        **provenance,
        "duplicate_classification_pipeline_git_commit": (
            stage10_commit
        ),
        "classification_basis": (
            "exact_full_bri_integer_milliangstrom"
        ),
        "paper_near_duplicate_threshold_mA": (
            PAPER_NEAR_DUPLICATE_THRESHOLD_MA
        ),
        "paper_near_duplicate_threshold_angstrom": (
            PAPER_NEAR_DUPLICATE_THRESHOLD_MA
            / 1000.0
        ),
        "input_pair_count": total["input"],
        "zero_duplicate_pair_count": (
            total["zero"]
        ),
        "paper_near_duplicate_pair_count": (
            total["paper_near"]
        ),
        "nonzero_near_duplicate_pair_count": (
            total["nonzero_near"]
        ),
        "not_near_duplicate_pair_count": (
            total["not_near"]
        ),
        "brain_candidate_input_pair_count": (
            candidate_counts["input"]
        ),
        "m1_input_pair_count": (
            m1_counts["input"]
        ),
        "brain_candidate_zero_duplicate_pair_count": (
            candidate_counts["zero"]
        ),
        "m1_zero_duplicate_pair_count": (
            m1_counts["zero"]
        ),
        "brain_candidate_paper_near_duplicate_pair_count": (
            candidate_counts["paper_near"]
        ),
        "m1_paper_near_duplicate_pair_count": (
            m1_counts["paper_near"]
        ),
        "processing_error_count": 0,
        "pair_accounting_valid": True,
        "filtering_performed": False,
    }

    _write_json_atomic(
        global_summary,
        output_root
        / "global_summary.json",
    )

    success10 = {
        "success_schema_name": (
            "pdbclean_stage10_duplicate_classification_success"
        ),
        "success_schema_version": "1.0",
        "snapshot": summary8["snapshot"],
        "cleaning_protocol": (
            summary8["cleaning_protocol"]
        ),
        **provenance,
        "duplicate_classification_pipeline_git_commit": (
            stage10_commit
        ),
        "global_summary": (
            "global_summary.json"
        ),
        "candidate_classifications": (
            "finalized/candidate_classifications.parquet"
        ),
        "m1_classifications": (
            "finalized/m1_classifications.parquet"
        ),
    }

    # Completion marker strictly last.
    _write_json_atomic(
        success10,
        success_path,
    )

    print(
        "Input pairs:",
        f"{total['input']:,}",
    )
    print(
        "Zero duplicates:",
        f"{total['zero']:,}",
    )
    print(
        "Paper near-duplicates (<=0.010 A):",
        f"{total['paper_near']:,}",
    )
    print(
        "Nonzero near-duplicates:",
        f"{total['nonzero_near']:,}",
    )
    print(
        "Not near-duplicates:",
        f"{total['not_near']:,}",
    )
    print()
    print(
        "STAGE-10 DUPLICATE CLASSIFICATION PUBLICATION: PASS"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
