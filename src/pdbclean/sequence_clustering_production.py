"""Stage 14d -- sequence-redundancy resolution over a completed Gold release.

Consumes the retained chains of the geometric release, clusters their retained
(post-trimming) sequences with MMseqs2, keeps one chain per cluster under the
project's own deterministic ranking, and publishes a new, separately identified
release in the same shape as the geometric one so downstream code needs no
special-casing.

The geometric release is read-only input. This stage never writes into it.

Entry point contract matches the other production stages:

    python -m pdbclean.sequence_clustering_production \
        --config <stage_config.yaml> --pipeline-git-commit <sha>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.sequence_clustering import (
    POLICY_NAME,
    POLICY_VERSION,
    RANKING,
    ChainQuality,
    SequenceClusteringError,
    choose_representatives,
    load_stage14_policy,
    parse_cluster_tsv,
    reconcile,
)
from pdbclean.stage_registry import release_name

SUCCESS_SCHEMA_NAME = "pdbclean_stage14d_sequence_clustering_success"
SUMMARY_SCHEMA_NAME = "pdbclean_stage14d_sequence_clustering_global_summary"
SCHEMA_VERSION = "1.0"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)

    return digest.hexdigest()


def _write_json_atomic(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def _write_parquet_atomic(table: pa.Table, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".parquet.tmp")

    if tmp.exists():
        tmp.unlink()

    pq.write_table(table, tmp, compression="zstd", version="2.6")
    tmp.replace(path)


def _validate_commit(value: str) -> str:
    text = value.strip()

    if len(text) != 40 or any(c not in "0123456789abcdef" for c in text.lower()):
        raise SequenceClusteringError(f"Not a full git commit sha: {value!r}")

    return text.lower()


def _run(cmd: list[str], *, label: str) -> None:
    print(f"  $ {' '.join(cmd[:3])} ... ({label})", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        sys.stderr.write(result.stdout[-4000:])
        sys.stderr.write(result.stderr[-4000:])
        raise SequenceClusteringError(
            f"{label} failed with exit code {result.returncode}"
        )


def _mmseqs_version(binary: str) -> str:
    out = subprocess.run([binary, "version"], capture_output=True, text=True)
    return out.stdout.strip() or "unknown"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--pipeline-git-commit", required=True)
    args = parser.parse_args()

    commit = _validate_commit(args.pipeline_git_commit)
    repo = Path.cwd()
    config = load_config(args.config).data

    section = config.get("sequence_clustering") or {}

    if not section.get("enabled", False):
        print("sequence_clustering.enabled is false; nothing to do.")
        return 0

    dataset = config["release"]["dataset_name"]
    protocol = config["release"]["protocol_version"]
    snapshot = str(config["snapshot"]["snapshot_id"])

    def _root(key: str, fallback: str) -> Path:
        value = Path(config["storage"].get(key, fallback))
        return value if value.is_absolute() else repo / value

    output_root = _root("output_root", "outputs/pdbclean")
    release_root = _root("release_root", "outputs/releases")

    source_release = release_root / release_name(
        dataset_name=dataset, snapshot=snapshot, protocol=protocol,
        suffix=str(section.get("input_release_suffix", "dedup-v1")),
    )
    target_release = release_root / release_name(
        dataset_name=dataset, snapshot=snapshot, protocol=protocol,
        suffix=str(section.get("release_suffix", "dedup-v1-seqclust-v1")),
    )
    stage_root = output_root / snapshot / protocol / "sequence_clustering"

    if target_release.resolve() == source_release.resolve():
        raise SequenceClusteringError(
            "release_suffix would overwrite the input release; the sequence "
            "stage must publish under its own identifier"
        )

    source_success = source_release / "_SUCCESS"

    if not source_success.is_file():
        raise SequenceClusteringError(
            f"Input Gold release is not complete: no {source_success}"
        )

    # ---- scientific parameters ---------------------------------------
    min_seq_id = float(section["min_seq_id"])
    coverage = float(section["coverage"])
    cov_mode = int(section["cov_mode"])
    cluster_mode = int(section["cluster_mode"])
    subcommand = str(section["subcommand"])
    binary = str(section.get("mmseqs_binary") or shutil.which("mmseqs") or "mmseqs")
    threads = int(config["execution"].get("sequence_clustering_threads", 8))

    if subcommand not in ("easy-cluster", "easy-linclust"):
        raise SequenceClusteringError(f"Unsupported subcommand {subcommand!r}")

    version = _mmseqs_version(binary)

    print("=== Stage 14d: sequence-redundancy resolution ===")
    print(f"  input release   {source_release}")
    print(f"  new release     {target_release}")
    print(f"  stage root      {stage_root}")
    print(f"  mmseqs          {binary} ({version})")
    print(f"  parameters      {subcommand} --min-seq-id {min_seq_id} "
          f"-c {coverage} --cov-mode {cov_mode} --cluster-mode {cluster_mode}")
    print()

    # ---- input population --------------------------------------------
    retained_path = source_release / "data/retained_chains.parquet"
    table = pq.read_table(retained_path)
    cols = table.to_pydict()

    keys = [f"{p}_{c}" for p, c in zip(cols["pdb_id"], cols["label_chain_id"])]

    if len(set(keys)) != len(keys):
        raise SequenceClusteringError(
            "pdb_id + label_chain_id is not unique over the input population"
        )

    quality = {
        key: ChainQuality(
            terminal_trimmed=bool(t),
            dirty_residue_count=d,
            pdb_id=p,
        )
        for key, p, t, d in zip(
            keys, cols["pdb_id"], cols["terminal_trimmed"],
            cols["dirty_residue_count"],
        )
    }
    sequences = dict(zip(keys, cols["retained_sequence"]))

    print(f"  input chains    {len(keys):,}")
    print(f"  distinct seqs   {len(set(sequences.values())):,}")

    # ---- deposition metadata, for the resolution term ------------------
    meta_path = (output_root / snapshot / protocol /
                 "downstream_metadata/finalized/entry_metadata.parquet")
    metadata: dict[str, dict[str, Any]] = {}

    if meta_path.is_file():
        mt = pq.read_table(meta_path).to_pylist()
        metadata = {row["pdb_id"]: row for row in mt}

    entries = {q.pdb_id for q in quality.values()}
    print(f"  metadata        {len(entries & set(metadata)):,} of "
          f"{len(entries):,} entries ({100 * len(entries & set(metadata)) / len(entries):.1f}%)")

    # ---- cluster -------------------------------------------------------
    work = stage_root / "work"

    if work.exists():
        shutil.rmtree(work)

    work.mkdir(parents=True, exist_ok=True)
    fasta = work / "input.fasta"

    # Sorted, so the FASTA the clusterer sees is a function of the population
    # alone and not of parquet row order.
    with open(fasta, "w", encoding="utf-8") as handle:
        for key in sorted(keys):
            handle.write(f">{key}\n{sequences[key]}\n")

    started = time.time()
    _run([binary, subcommand, str(fasta), str(work / "clu"), str(work / "tmp"),
          "--min-seq-id", repr(min_seq_id), "-c", repr(coverage),
          "--cov-mode", str(cov_mode), "--cluster-mode", str(cluster_mode),
          "--threads", str(threads), "-v", "1"], label="mmseqs clustering")
    cluster_seconds = time.time() - started

    clusters = parse_cluster_tsv(work / "clu_cluster.tsv")
    assigned = sum(len(v) for v in clusters.values())

    if assigned != len(keys):
        raise SequenceClusteringError(
            f"MMseqs2 assigned {assigned:,} chains, expected {len(keys):,}"
        )

    print(f"  clusters        {len(clusters):,}  "
          f"(clustering took {cluster_seconds:.0f}s)")

    # ---- our representative selection ----------------------------------
    policy = load_stage14_policy(repo)
    rep_of_cluster, rep_of_chain, resolution_clusters = choose_representatives(
        clusters, quality, metadata, policy
    )

    retained_keys = sorted(set(rep_of_cluster.values()))
    removed_keys = sorted(set(keys) - set(retained_keys))

    gates = reconcile(
        input_chains=keys, retained=retained_keys, removed=removed_keys,
        representative_of_chain=rep_of_chain,
    )
    print(f"  retained        {len(retained_keys):,}")
    print(f"  removed         {len(removed_keys):,}")
    print(f"  clusters where the resolution term applied: "
          f"{resolution_clusters:,} / {len(clusters):,}")

    # ---- pairwise identity, removed -> its representative ---------------
    identity, identity_method = _pair_identity(
        binary, work, sequences, removed_keys, rep_of_chain, retained_keys,
        min_seq_id=min_seq_id, coverage=coverage, cov_mode=cov_mode,
        threads=threads,
    )

    # ---- manifests ------------------------------------------------------
    index = {key: i for i, key in enumerate(keys)}
    cluster_of_chain = {
        member: cid for cid, members in clusters.items() for member in members
    }

    retained_idx = [index[k] for k in retained_keys]
    removed_idx = [index[k] for k in removed_keys]

    retained_table = table.take(retained_idx)
    removed_table = table.take(removed_idx)

    removed_table = removed_table.append_column(
        "seqdedup_cluster_id", pa.array([cluster_of_chain[k] for k in removed_keys])
    ).append_column(
        "seqdedup_representative_pdb_id",
        pa.array([rep_of_chain[k].rsplit("_", 1)[0] for k in removed_keys]),
    ).append_column(
        "seqdedup_representative_label_chain_id",
        pa.array([rep_of_chain[k].rsplit("_", 1)[1] for k in removed_keys]),
    ).append_column(
        "seqdedup_sequence_identity",
        pa.array([identity.get(k) for k in removed_keys], type=pa.float64()),
    ).append_column(
        "seqdedup_identity_method",
        pa.array([identity_method.get(k, "unresolved") for k in removed_keys]),
    ).append_column(
        "seqdedup_policy_version", pa.array([POLICY_VERSION] * len(removed_keys))
    )

    _write_parquet_atomic(retained_table, stage_root / "finalized/retained_chains.parquet")
    _write_parquet_atomic(removed_table, stage_root / "finalized/removed_chain_audit.parquet")

    resolved = sum(1 for k in removed_keys if identity.get(k) is not None)
    methods: dict[str, int] = {}

    for k in removed_keys:
        m = identity_method.get(k, "unresolved")
        methods[m] = methods.get(m, 0) + 1

    bands: dict[str, int] = {}

    for k in removed_keys:
        v = identity.get(k)
        band = ("unresolved" if v is None else "100%" if v >= 0.999999
                else ">=90%" if v >= 0.90 else ">=70%" if v >= 0.70
                else ">=50%" if v >= 0.50 else "<50%")
        bands[band] = bands.get(band, 0) + 1

    sizes = sorted(len(v) for v in clusters.values())
    size_hist: dict[str, int] = {}

    for s in sizes:
        band = ("1" if s == 1 else "2" if s == 2 else "3-10" if s <= 10
                else "11-100" if s <= 100 else "101-1000" if s <= 1000 else ">1000")
        size_hist[band] = size_hist.get(band, 0) + 1

    summary = {
        "summary_schema_name": SUMMARY_SCHEMA_NAME,
        "summary_schema_version": SCHEMA_VERSION,
        "snapshot": snapshot,
        "cleaning_protocol": protocol,
        "input_release": source_release.name,
        "release_name": target_release.name,
        "mmseqs_version": version,
        "mmseqs_subcommand": subcommand,
        "min_seq_id": min_seq_id,
        "coverage": coverage,
        "cov_mode": cov_mode,
        "cluster_mode": cluster_mode,
        "policy_name": POLICY_NAME,
        "policy_version": POLICY_VERSION,
        "ranking": list(RANKING),
        "cluster_count": len(clusters),
        "cluster_size_histogram": size_hist,
        "largest_cluster": sizes[-1] if sizes else 0,
        "resolution_term_clusters": resolution_clusters,
        "identity_resolved": resolved,
        "identity_bands": bands,
        "identity_methods": methods,
        "clustering_seconds": round(cluster_seconds, 1),
        **gates,
    }
    _write_json_atomic(summary, stage_root / "global_summary.json")

    # ---- publish the new release ----------------------------------------
    _publish(target_release, stage_root, source_release, summary, commit)

    success = {
        "success_schema_name": SUCCESS_SCHEMA_NAME,
        "success_schema_version": SCHEMA_VERSION,
        "snapshot": snapshot,
        "cleaning_protocol": protocol,
        "release_name": target_release.name,
        "input_release": source_release.name,
        "sequence_clustering_pipeline_git_commit": commit,
        "min_seq_id": min_seq_id,
        "coverage": coverage,
        "cov_mode": cov_mode,
        "cluster_mode": cluster_mode,
        "policy_version": POLICY_VERSION,
        "global_summary": "global_summary.json",
        "retained_chains": "finalized/retained_chains.parquet",
        "removed_chain_audit": "finalized/removed_chain_audit.parquet",
    }
    _write_json_atomic(success, stage_root / "_SUCCESS")

    shutil.rmtree(work, ignore_errors=True)

    print()
    print(f"  cluster sizes   {size_hist}")
    print(f"  identity bands  {bands}")
    print(f"  identity method {methods}")
    print()
    print("STAGE-15 SEQUENCE CLUSTERING PUBLICATION: PASS")
    return 0


def _pair_identity(binary, work, sequences, removed_keys, rep_of_chain,
                   retained_keys, *, min_seq_id, coverage, cov_mode, threads):
    """Alignment identity from each removed chain to its chosen representative.

    Computed with MMseqs2 rather than positionally: 55% of cluster members
    differ in length from their representative under `--cov-mode 0 -c 0.8`, so
    an ungapped position-by-position comparison would be wrong for the
    majority of pairs.
    """

    if not removed_keys:
        return {}

    qf = work / "identity_query.fasta"
    tf = work / "identity_target.fasta"

    with open(qf, "w", encoding="utf-8") as handle:
        for key in removed_keys:
            handle.write(f">{key}\n{sequences[key]}\n")

    with open(tf, "w", encoding="utf-8") as handle:
        for key in retained_keys:
            handle.write(f">{key}\n{sequences[key]}\n")

    out = work / "identity.tsv"
    _run([binary, "easy-search", str(qf), str(tf), str(out), str(work / "tmp_id"),
          "--min-seq-id", repr(min_seq_id), "-c", repr(coverage),
          "--cov-mode", str(cov_mode), "--threads", str(threads),
          "--max-seqs", "300", "--format-output", "query,target,fident",
          "-v", "1"], label="mmseqs identity search")

    wanted = {(k, rep_of_chain[k]) for k in removed_keys}
    identity: dict[str, float] = {}
    method: dict[str, str] = {}

    with open(out, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")

            if len(parts) < 3:
                continue

            q, t, f = parts[0], parts[1], parts[2]

            if (q, t) in wanted and q not in identity:
                identity[q] = float(f)
                method[q] = "mmseqs_alignment"

    # MMseqs2's k-mer prefilter cannot seed very short sequences, so a minority
    # of pairs come back with no alignment -- in the 2026-01-01 population
    # 11,199 of 415,614, median length 10 residues. Leaving those null would
    # break the requirement that every removal records its identity, so they
    # are resolved exactly and deterministically instead of being dropped:
    # an identical pair is 1.0 by inspection, and an equal-length pair in the
    # same cluster is compared position by position. Only a pair that is
    # neither is left unresolved, and the method used is recorded per row so
    # no figure is mistaken for an alignment that did not happen.
    for key in removed_keys:
        if key in identity:
            continue

        query = sequences[key]
        target = sequences[rep_of_chain[key]]

        if query == target:
            identity[key] = 1.0
            method[key] = "exact_match"
            continue

        if len(query) == len(target) and query:
            same = sum(1 for a, b in zip(query, target) if a == b)
            identity[key] = same / len(query)
            method[key] = "ungapped_positional"
            continue

        short, long_ = (query, target) if len(query) <= len(target) else (target, query)

        if not short:
            method[key] = "unresolved_empty_sequence"
            continue

        best = max(
            sum(1 for a, b in zip(short, long_[offset:]) if a == b)
            for offset in range(len(long_) - len(short) + 1)
        )
        identity[key] = best / len(short)
        method[key] = "ungapped_best_offset"

    return identity, method


def _publish(target: Path, stage_root: Path, source: Path,
             summary: dict[str, Any], commit: str) -> None:
    """Write the new release tree, in the same shape as the geometric release."""

    if target.exists():
        shutil.rmtree(target)

    (target / "data").mkdir(parents=True, exist_ok=True)
    (target / "audit").mkdir(parents=True, exist_ok=True)
    (target / "provenance").mkdir(parents=True, exist_ok=True)

    shutil.copy2(stage_root / "finalized/retained_chains.parquet",
                 target / "data/retained_chains.parquet")
    shutil.copy2(stage_root / "finalized/removed_chain_audit.parquet",
                 target / "audit/removed_chain_audit.parquet")
    shutil.copy2(stage_root / "global_summary.json",
                 target / "provenance/sequence_clustering_summary.json")

    # The geometric release this one is derived from, by identity not by copy.
    src_manifest = source / "release_manifest.json"
    derived = {
        "release_name": source.name,
        "release_manifest_sha256": _sha256(src_manifest) if src_manifest.is_file() else None,
        "retained_chains_sha256": _sha256(source / "data/retained_chains.parquet"),
    }
    _write_json_atomic(derived, target / "provenance/derived_from.json")

    files = []

    for path in sorted(target.rglob("*")):
        if path.is_file():
            files.append({
                "path": str(path.relative_to(target)),
                "size_bytes": path.stat().st_size,
                "sha256": _sha256(path),
            })

    manifest = {
        "release_name": target.name,
        "release_version": "1.0",
        "release_kind": "geometric_then_sequence",
        "snapshot": summary["snapshot"],
        "protocol": summary["cleaning_protocol"],
        "derived_from_release": source.name,
        "sequence_clustering": {
            "mmseqs_version": summary["mmseqs_version"],
            "subcommand": summary["mmseqs_subcommand"],
            "min_seq_id": summary["min_seq_id"],
            "coverage": summary["coverage"],
            "cov_mode": summary["cov_mode"],
            "cluster_mode": summary["cluster_mode"],
            "policy_name": summary["policy_name"],
            "policy_version": summary["policy_version"],
            "ranking": summary["ranking"],
        },
        "input_chain_count": summary["input_chain_count"],
        "retained_chain_count": summary["retained_chain_count"],
        "removed_chain_count": summary["removed_chain_count"],
        "pipeline_git_commit": commit,
        "files": files,
    }
    _write_json_atomic(manifest, target / "release_manifest.json")
    _write_json_atomic(
        {"success_schema_name": "pdbclean_stage14d_release_success",
         "success_schema_version": SCHEMA_VERSION,
         "release_name": target.name,
         "retained_chain_count": summary["retained_chain_count"],
         "removed_chain_count": summary["removed_chain_count"]},
        target / "_SUCCESS",
    )


if __name__ == "__main__":
    raise SystemExit(main())
