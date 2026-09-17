"""Project a release into OpenFold training inputs PLUS a held-out validation set.

The training half is the same projection `build_seqclust_training_view.py`
produces. The validation half is different in kind: OpenFold's validation path
takes a data directory and a per-chain alignment directory, not an alignment
index, so the validation chains get the conventional layout

    val_data_dir/<pdb_id>.cif
    val_alignment_dir/<pdb_id>_<auth_chain_id>/msa.a3m

built from symlinks into the existing mmCIF tree and the content-addressed MSA
store. Nothing is copied: a validation set is a view of data that already
exists, and duplicating 4,000 alignments would only create a second thing to
keep in step.

The split itself is supplied, not computed here -- it is chosen so that no
validation chain shares a geometric near-duplicate edge or an identical
sequence with any training chain, which is a property of the redundancy graph
rather than of this projection.

Read-only with respect to the release, the mmCIF tree and the MSA store.
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
    d = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            d.update(b)
    return d.hexdigest()


def run(cmd: list[str], label: str) -> None:
    print(f"\n--- {label}")
    if subprocess.run([str(c) for c in cmd], text=True).returncode:
        raise SystemExit(f"{label} failed")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", required=True, type=Path)
    ap.add_argument("--split", required=True, type=Path,
                    help="JSON with 'validation' and 'training' chain-key lists")
    ap.add_argument("--msa-store", required=True, type=Path)
    ap.add_argument("--mmcif-dir", required=True, type=Path)
    ap.add_argument("--output-root", required=True, type=Path)
    ap.add_argument("--python", default=sys.executable)
    args = ap.parse_args()

    split = json.loads(args.split.read_text())
    val_keys, train_keys = set(split["validation"]), set(split["training"])

    table = pq.read_table(args.release / "data/retained_chains.parquet")
    cols = table.to_pydict()
    label_keys = [f"{p}_{c}" for p, c in zip(cols["pdb_id"], cols["label_chain_id"])]
    auth_ids = [f"{p}_{a}" for p, a in zip(cols["pdb_id"], cols["auth_chain_id"])]
    digests = [hashlib.sha256(s.encode()).hexdigest() for s in cols["retained_sequence"]]

    if len(set(auth_ids)) != len(auth_ids):
        raise SystemExit("openfold_chain_id is not unique over this release")

    print(f"release   {args.release.name}")
    print(f"train {len(train_keys):,}   validation {len(val_keys):,}")

    # ---- MSA coverage for BOTH halves, before anything is written ------
    missing = [d for d in set(digests)
               if not (args.msa_store / f"{d}.a3m").is_file()]
    if missing:
        raise SystemExit(f"{len(missing)} sequences have no MSA")
    print(f"MSAs present for all {len(set(digests)):,} distinct sequences")

    args.output_root.mkdir(parents=True, exist_ok=True)
    train_idx = [i for i, k in enumerate(label_keys) if k in train_keys]
    val_idx = [i for i, k in enumerate(label_keys) if k in val_keys]

    if len(train_idx) + len(val_idx) != len(label_keys):
        raise SystemExit("split does not partition the release")

    # ---- training half --------------------------------------------------
    tt = table.take(train_idx).append_column(
        "openfold_chain_id", pa.array([auth_ids[i] for i in train_idx]))
    chains_path = args.output_root / "openfold_training_chains.parquet"
    pq.write_table(tt, chains_path, compression="zstd", version="2.6")

    filter_path = args.output_root / "train_filter.txt"
    filter_path.write_text(
        "\n".join(sorted(auth_ids[i] for i in train_idx)) + "\n")

    chain_map = args.output_root / "chain_to_sequence.parquet"
    pq.write_table(
        pa.table({"openfold_chain_id": pa.array([auth_ids[i] for i in train_idx]),
                  "sequence_sha256": pa.array([digests[i] for i in train_idx])}),
        chain_map, compression="zstd", version="2.6")

    index_dir = args.output_root / "alignment_index"
    index_dir.mkdir(parents=True, exist_ok=True)
    index_path = index_dir / "pdbclean_final.index"
    run([args.python, SCRIPTS / "build_alignment_index.py",
         "--chain-map", chain_map, "--msa-store", args.msa_store,
         "--output", index_path,
         "--manifest", index_dir / "alignment_index_manifest.json",
         "--expected-chains", len(train_idx),
         "--expected-sequences", len({digests[i] for i in train_idx})],
        "training alignment index")

    view = args.output_root / "training_view"
    view.mkdir(parents=True, exist_ok=True)
    run([args.python, SCRIPTS / "build_projection_index.py",
         "--training-chains", chains_path,
         "--output", view / "projection_index.json",
         "--manifest", view / "projection_index_manifest.json",
         "--expected-chains", len(train_idx)],
        "projection index")
    filter_view = view / "train_filter_view.txt"
    filter_view.write_text(filter_path.read_text())

    # ---- validation half ------------------------------------------------
    val_data = args.output_root / "val_data_dir"
    val_align = args.output_root / "val_alignment_dir"
    val_data.mkdir(parents=True, exist_ok=True)
    val_align.mkdir(parents=True, exist_ok=True)

    entries, linked = set(), 0
    for i in val_idx:
        pdb = cols["pdb_id"][i]
        entries.add(pdb)
        d = val_align / auth_ids[i]
        d.mkdir(parents=True, exist_ok=True)
        target = args.msa_store / f"{digests[i]}.a3m"
        link = d / "msa.a3m"
        if not link.exists():
            link.symlink_to(target)
        linked += 1

    for pdb in sorted(entries):
        src = args.mmcif_dir / f"{pdb}.cif"
        if not src.is_file():
            raise SystemExit(f"validation entry {pdb} has no mmCIF at {src}")
        dst = val_data / f"{pdb}.cif"
        if not dst.exists():
            dst.symlink_to(src)

    print(f"\nvalidation: {linked:,} alignment dirs, {len(entries):,} mmCIF entries")

    # cache restricted to the validation chains
    cache_src = args.output_root.parent / "pdbclean-dedup-v1/training_view/train_chain_data_cache.json"
    cache = json.loads(cache_src.read_text())
    val_cache = {auth_ids[i]: cache[auth_ids[i]] for i in val_idx
                 if auth_ids[i] in cache}
    if len(val_cache) != len(val_idx):
        raise SystemExit(
            f"chain data cache covers {len(val_cache):,} of {len(val_idx):,} "
            "validation chains")
    val_cache_path = args.output_root / "val_mmcif_data_cache.json"
    val_cache_path.write_text(json.dumps(val_cache))

    train_cache = {auth_ids[i]: cache[auth_ids[i]] for i in train_idx}
    (view / "train_chain_data_cache.json").write_text(json.dumps(train_cache))

    manifest = {
        "stage": "openfold_training_materialisation_with_validation",
        "source_release": args.release.name,
        "openfold_identifier": "<pdb_id>_<auth_chain_id>",
        "training_chain_count": len(train_idx),
        "validation_chain_count": len(val_idx),
        "validation_entry_count": len(entries),
        "unique_sequence_count": len(set(digests)),
        "msas_reused": len(set(digests)),
        "msas_generated": 0,
        "alignment_index_sha256": sha256_file(index_path),
        "train_filter_sha256": sha256_file(filter_path),
        "split_file": str(args.split),
    }
    (args.output_root / "training_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True))

    print(f"\nREADY: {len(train_idx):,} training, {len(val_idx):,} validation, "
          f"0 MSAs generated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
