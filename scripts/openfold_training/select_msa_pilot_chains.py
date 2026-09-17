"""Select a representative Stage-18 MSA pilot from a run's query population.

The pilot exists to prove the MSA machinery before ~500k production searches
are committed.  It is therefore chosen to hit the cases that actually break
things, not a random sample:

    length 1                  the degenerate minimum a real pipeline must survive
    very short                short enough that profile search behaves unusually
    normal                    the median case
    long                      cost-representative
    near maximum / maximum    the worst-case runtime and memory envelope
    label != auth chain id    the identifier-namespace hazard
    identical sequences       two DIFFERENT chains sharing one exact sequence

That last case is the important one.  Fewer unique sequences exist than chains,
and a naive implementation "optimises" by searching each unique sequence once
and fanning the result out.  This experiment forbids that: every retained chain
gets its own independent search and its own MSA identity.  Including a sequence
-sharing pair in the pilot makes the shortcut observable rather than assumed --
if the two chains come back with a shared MSA identity, the pilot fails.

Nothing here is snapshot-specific.  Every path is an argument, so the same
selector serves any run.  The selection is deterministic: given the same query
population it always returns the same pilot, so a pilot result stays meaningful
evidence for the production run that follows it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def _pick(table: pa.Table, mask: pa.Array, reason: str) -> dict | None:
    """Return the first row matching ``mask``, or ``None``.

    "First" is by ``query_index``, which is stable for a frozen query
    population, so the pilot does not drift between invocations.
    """

    indices = pc.indices_nonzero(mask)

    if len(indices) == 0:
        return None

    position = indices[0].as_py()

    record = {
        name: table.column(name)[position].as_py()
        for name in table.column_names
    }
    record["pilot_reason"] = reason

    return record


def select_pilot(
    *,
    queries_path: Path,
    training_path: Path,
    dedup_path: Path,
) -> dict:
    queries = ds.dataset(queries_path).to_table(
        columns=[
            "query_index",
            "openfold_chain_id",
            "pdb_id",
            "sequence_length",
            "sequence",
        ]
    )

    # Sort once so "first match" is defined by query_index, not row order.
    queries = queries.sort_by([("query_index", "ascending")])

    training = ds.dataset(training_path).to_table(
        columns=[
            "openfold_chain_id",
            "pdb_id",
            "label_chain_id",
            "auth_chain_id",
            "retained_residue_count",
            "terminal_trimmed",
            "dirty_residue_count",
        ]
    )

    lengths = queries.column("sequence_length")
    minimum = pc.min(lengths).as_py()
    maximum = pc.max(lengths).as_py()

    # Median without materialising a sorted copy of every length.
    median = int(
        pc.approximate_median(pc.cast(lengths, pa.float64())).as_py()
    )

    selected: list[dict] = []
    seen: set[str] = set()

    def take(record: dict | None) -> None:
        if record is None:
            return

        key = record["openfold_chain_id"]

        if key in seen:
            return

        seen.add(key)
        selected.append(record)

    # --- length envelope -------------------------------------------------

    take(
        _pick(
            queries,
            pc.equal(lengths, minimum),
            f"minimum length ({minimum})",
        )
    )

    take(
        _pick(
            queries,
            pc.and_(
                pc.greater_equal(lengths, 2),
                pc.less_equal(lengths, 10),
            ),
            "very short (2-10 residues)",
        )
    )

    take(
        _pick(
            queries,
            pc.equal(lengths, median),
            f"median length ({median})",
        )
    )

    take(
        _pick(
            queries,
            pc.and_(
                pc.greater_equal(lengths, 900),
                pc.less_equal(lengths, 1100),
            ),
            "long (900-1100 residues)",
        )
    )

    take(
        _pick(
            queries,
            pc.and_(
                pc.greater_equal(
                    lengths,
                    int(maximum * 0.75),
                ),
                pc.less(lengths, maximum),
            ),
            "near maximum length",
        )
    )

    take(
        _pick(
            queries,
            pc.equal(lengths, maximum),
            f"maximum length ({maximum})",
        )
    )

    # --- identifier-namespace hazard -------------------------------------

    # openfold_chain_id is built from auth_chain_id, while the scientific
    # pipeline computed BRI on label_chain_id.  A chain where the two differ
    # proves the pilot is carrying the right identifier end to end.
    diverges = pc.and_(
        pc.is_valid(training.column("auth_chain_id")),
        pc.not_equal(
            training.column("label_chain_id"),
            training.column("auth_chain_id"),
        ),
    )

    diverging_ids = training.filter(diverges).column("openfold_chain_id")

    take(
        _pick(
            queries,
            pc.is_in(
                queries.column("openfold_chain_id"),
                value_set=diverging_ids,
            ),
            "label_chain_id differs from auth_chain_id",
        )
    )

    # --- sequence-sharing hazard -----------------------------------------

    # Two DIFFERENT retained chains carrying one identical sequence.  Both are
    # searched independently in production; the pilot must show two separate
    # MSA identities, never one reused twice.
    dedup = ds.dataset(dedup_path).to_table(
        columns=[
            "sequence_sha256",
            "openfold_chain_id",
            "query_index",
        ]
    )

    shared_pair: list[dict] = []

    duplicated = dedup.group_by("sequence_sha256").aggregate(
        [("openfold_chain_id", "count")]
    )

    multiplicity = duplicated.column("openfold_chain_id_count")

    # A modest multiplicity keeps the pilot honest without picking the most
    # extreme fan-out in the dataset.
    candidates = duplicated.filter(
        pc.and_(
            pc.greater_equal(multiplicity, 2),
            pc.less_equal(multiplicity, 4),
        )
    )

    if candidates.num_rows > 0:
        target_hash = candidates.column("sequence_sha256")[0].as_py()

        members = dedup.filter(
            pc.equal(dedup.column("sequence_sha256"), target_hash)
        ).sort_by([("query_index", "ascending")])

        member_ids = members.column("openfold_chain_id")[:2].to_pylist()

        for member_id in member_ids:
            record = _pick(
                queries,
                pc.equal(queries.column("openfold_chain_id"), member_id),
                (
                    "identical sequence shared with another retained chain "
                    f"(sha256 {target_hash[:16]}...)"
                ),
            )

            if record is not None:
                record["shared_sequence_sha256"] = target_hash
                shared_pair.append(record)
                take(record)

    # --- provenance -------------------------------------------------------

    manifest = {
        "stage": "openfold_msa_pilot_selection",
        "stage_version": "1.0",
        "selection_is_deterministic": True,
        "query_population_size": queries.num_rows,
        "query_manifest": str(queries_path),
        "query_manifest_sha256": _sha256(queries_path),
        "training_manifest": str(training_path),
        "training_manifest_sha256": _sha256(training_path),
        "sequence_dedup_manifest": str(dedup_path),
        "sequence_dedup_manifest_sha256": _sha256(dedup_path),
        "length_minimum": minimum,
        "length_median": median,
        "length_maximum": maximum,
        "pilot_chain_count": len(selected),
        "sequence_sharing_pair": [
            record["openfold_chain_id"] for record in shared_pair
        ],
        "sequence_sharing_pair_is_present": len(shared_pair) == 2,
        "note": (
            "Sequence-level statistics are audit information only. Production "
            "Stage 18 searches every retained chain independently; identical "
            "sequences never share an MSA identity."
        ),
        "chains": [
            {
                key: value
                for key, value in record.items()
                if key != "sequence"
            }
            for record in selected
        ],
    }

    return {"manifest": manifest, "records": selected}


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Select a deterministic, representative Stage-18 MSA pilot from a "
            "run's existing query population."
        )
    )
    parser.add_argument("--queries", required=True, type=Path)
    parser.add_argument("--training-manifest", required=True, type=Path)
    parser.add_argument("--sequence-dedup", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)

    args = parser.parse_args()

    result = select_pilot(
        queries_path=args.queries,
        training_path=args.training_manifest,
        dedup_path=args.sequence_dedup,
    )

    manifest = result["manifest"]
    records = result["records"]

    args.output_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = args.output_dir / "pilot_queries.fasta"

    with fasta_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record['openfold_chain_id']}\n")
            handle.write(f"{record['sequence']}\n")

    manifest["pilot_fasta"] = str(fasta_path)
    manifest["pilot_fasta_sha256"] = _sha256(fasta_path)

    manifest_path = args.output_dir / "pilot_selection_manifest.json"

    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(manifest, indent=2, sort_keys=True))

    print()
    print(f"pilot FASTA:    {fasta_path}")
    print(f"pilot manifest: {manifest_path}")

    if not manifest["sequence_sharing_pair_is_present"]:
        print()
        print(
            "WARNING: no identical-sequence chain pair was selected; the "
            "pilot cannot prove absence of sequence-sharing reuse."
        )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
