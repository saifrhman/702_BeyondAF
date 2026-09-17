"""Project a Stage-14d (geometry-then-sequence) release into OpenFold inputs.

The surviving chains are a subset of the geometric release's retained chains,
and MSAs are content-addressed by sequence SHA256, so alignment generation is a
re-projection of existing resources rather than a regeneration. This script
asserts that rather than assuming it: if any surviving sequence lacks an MSA it
stops, and says how many.

It produces, under a training root of its own:

    openfold_training_chains.parquet   retained chains + openfold_chain_id
    train_filter.txt                   the chain ids, one per line
    chain_to_sequence.parquet          openfold_chain_id -> sequence_sha256
    alignment_index/<name>.index       built by build_alignment_index.py
    training_view/projection_index.json built by build_projection_index.py
    training_manifest.json             provenance for the whole projection

Read-only with respect to the release it projects and to the MSA store.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

SCRIPTS = Path(__file__).resolve().parent


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)

    return digest.hexdigest()


def run(cmd: list[str], label: str) -> None:
    print(f"\n--- {label}")
    result = subprocess.run([str(c) for c in cmd], text=True)

    if result.returncode != 0:
        raise SystemExit(f"{label} failed ({result.returncode})")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", required=True, type=Path,
                    help="Stage-14d release directory (read-only)")
    ap.add_argument("--msa-store", required=True, type=Path,
                    help="Directory of <sequence_sha256>.a3m")
    ap.add_argument("--output-root", required=True, type=Path)
    ap.add_argument("--python", default=sys.executable)
    args = ap.parse_args()

    retained_path = args.release / "data/retained_chains.parquet"
    table = pq.read_table(retained_path)
    cols = table.to_pydict()
    n = table.num_rows

    print(f"release        {args.release.name}")
    print(f"retained       {n:,} chains")

    # OpenFold keys a chain by AUTHOR chain id, not label_asym_id. Getting this
    # wrong silently mislabels the majority of chains, so it is asserted.
    chain_ids = [f"{p}_{a}" for p, a in zip(cols["pdb_id"], cols["auth_chain_id"])]

    if len(set(chain_ids)) != n:
        dupes = [c for c in set(chain_ids) if chain_ids.count(c) > 1][:5]
        raise SystemExit(
            f"openfold_chain_id is not unique over this population "
            f"({n - len(set(chain_ids))} collisions, e.g. {dupes})"
        )

    args.output_root.mkdir(parents=True, exist_ok=True)

    # ---- MSA coverage, before anything is written ---------------------
    digests = [hashlib.sha256(s.encode()).hexdigest() for s in cols["retained_sequence"]]
    unique = sorted(set(digests))
    missing = [d for d in unique if not (args.msa_store / f"{d}.a3m").is_file()]

    print(f"distinct seqs  {len(unique):,}")
    print(f"MSAs present   {len(unique) - len(missing):,} "
          f"({100 * (len(unique) - len(missing)) / len(unique):.3f}%)")
    print(f"MSAs missing   {len(missing):,}")

    if missing:
        raise SystemExit(
            f"{len(missing)} sequences have no MSA; this population is not a "
            "pure re-projection and alignment generation would be required."
        )

    # ---- training manifest -------------------------------------------
    chains_table = table.append_column("openfold_chain_id", pa.array(chain_ids))
    chains_path = args.output_root / "openfold_training_chains.parquet"
    pq.write_table(chains_table, chains_path, compression="zstd", version="2.6")

    filter_path = args.output_root / "train_filter.txt"
    filter_path.write_text("\n".join(sorted(chain_ids)) + "\n", encoding="utf-8")

    chain_map_path = args.output_root / "chain_to_sequence.parquet"
    pq.write_table(
        pa.table({"openfold_chain_id": pa.array(chain_ids),
                  "sequence_sha256": pa.array(digests)}),
        chain_map_path, compression="zstd", version="2.6",
    )

    # ---- alignment index ----------------------------------------------
    index_dir = args.output_root / "alignment_index"
    index_dir.mkdir(parents=True, exist_ok=True)
    index_path = index_dir / "pdbclean_seqid100_v1.index"

    run([args.python, SCRIPTS / "build_alignment_index.py",
         "--chain-map", chain_map_path, "--msa-store", args.msa_store,
         "--output", index_path,
         "--manifest", index_dir / "alignment_index_manifest.json",
         "--expected-chains", n, "--expected-sequences", len(unique)],
        "alignment index")

    # ---- projection index ----------------------------------------------
    view = args.output_root / "training_view"
    view.mkdir(parents=True, exist_ok=True)

    run([args.python, SCRIPTS / "build_projection_index.py",
         "--training-chains", chains_path,
         "--output", view / "projection_index.json",
         "--manifest", view / "projection_index_manifest.json",
         "--expected-chains", n],
        "projection index")

    # ---- coverage gate --------------------------------------------------
    run([args.python, SCRIPTS / "verify_alignment_index.py",
         "--index", index_path, "--alignment-dir", args.msa_store,
         "--training-chains", chains_path, "--train-filter", filter_path,
         "--verify-chains", "2000", "--seed", "20260101",
         "--report", args.output_root / "alignment_index_verification.json"],
        "alignment index verification")

    manifest = {
        "stage": "openfold_training_materialisation",
        "stage_version": "1.0",
        "snapshot": "20260101",
        "source_release": args.release.name,
        "source_retained_chains_path": str(retained_path),
        "source_retained_chains_sha256": sha256_file(retained_path),
        "openfold_identifier": "<pdb_id>_<auth_chain_id>",
        "training_chain_count": n,
        "unique_pdb_entry_count": len(set(cols["pdb_id"])),
        "unique_sequence_count": len(unique),
        "msas_reused": len(unique),
        "msas_generated": 0,
        "msa_store": str(args.msa_store),
        "openfold_training_manifest": chains_path.name,
        "openfold_training_manifest_sha256": sha256_file(chains_path),
        "train_filter": filter_path.name,
        "train_filter_sha256": sha256_file(filter_path),
        "alignment_index": str(index_path.relative_to(args.output_root)),
        "alignment_index_sha256": sha256_file(index_path),
    }
    (args.output_root / "training_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )

    print(f"\nTRAINING VIEW READY: {n:,} chains, {len(unique):,} MSAs reused, "
          "0 generated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
