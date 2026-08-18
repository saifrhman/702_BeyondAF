#!/usr/bin/env python3

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import pickle
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import yaml
from gemmi import cif


CACHE_VERSION = "comp702-nonzero-near-evidence-v1"


# ============================================================
# Basic helpers
# ============================================================

def P(a, b):
    return tuple(sorted((str(a).lower(), str(b).lower())))


def clean(v):
    if v is None:
        return ""

    s = str(v).strip()

    if s in {"?", "."}:
        return ""

    return s.strip("'\"")


def fnum(v):
    try:
        return float(v)
    except Exception:
        return None


def distinct_join(values):
    out = []

    for value in values:
        value = clean(value)

        if value and value not in out:
            out.append(value)

    return "|".join(out)


def values(block, tag):
    col = block.find_values(tag)

    if not col:
        return []

    out = []

    for i in range(len(col)):
        try:
            value = col.str(i)
        except Exception:
            value = col[i]

        value = clean(value)

        if value:
            out.append(value)

    return out


def first(block, *tags):
    for tag in tags:
        vals = values(block, tag)

        if vals:
            return vals[0]

    return ""


def sha_rows(rows):
    h = hashlib.sha256()

    for row in sorted(rows):
        if isinstance(row, (tuple, list)):
            text = "\x1f".join(map(str, row))
        else:
            text = str(row)

        h.update(text.encode("utf-8"))
        h.update(b"\n")

    return h.hexdigest()


def safe_float_delta(a, b):
    fa = fnum(a)
    fb = fnum(b)

    if fa is None or fb is None:
        return None

    return abs(fa - fb)


# ============================================================
# Rigid-body alignment
# ============================================================

def kabsch_metrics(target_xyz, source_xyz):
    """
    Align source onto target using a proper rotation + translation.

    Returns RMSD and maximum per-atom Euclidean displacement.
    Reflection is not permitted.
    """

    a = np.asarray(target_xyz, dtype=np.float64)
    b = np.asarray(source_xyz, dtype=np.float64)

    if a.shape != b.shape:
        raise ValueError(
            f"Kabsch shape mismatch: {a.shape} vs {b.shape}"
        )

    if a.ndim != 2 or a.shape[1] != 3:
        raise ValueError(
            f"Expected Nx3 coordinates, got {a.shape}"
        )

    if len(a) == 0:
        raise ValueError("No coordinates supplied")

    ca = a.mean(axis=0)
    cb = b.mean(axis=0)

    aa = a - ca
    bb = b - cb

    h = bb.T @ aa

    u, _, vt = np.linalg.svd(h)

    r = u @ vt

    if np.linalg.det(r) < 0:
        u[:, -1] *= -1
        r = u @ vt

    aligned = bb @ r + ca

    delta = aligned - a
    distances = np.linalg.norm(delta, axis=1)

    rmsd = math.sqrt(
        float(np.mean(np.sum(delta * delta, axis=1)))
    )

    return {
        "rmsd_A": rmsd,
        "max_atom_deviation_A": float(distances.max()),
        "mean_atom_deviation_A": float(distances.mean()),
    }


# ============================================================
# mmCIF source fetch/cache
# ============================================================

def fetch_bytes(url, expected_size, attempts=3):
    last = None

    for attempt in range(1, attempts + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={
                    "User-Agent":
                        "COMP702-2026-nonzero-near-analysis/1.0"
                },
            )

            with urllib.request.urlopen(
                request,
                timeout=120,
            ) as response:
                data = response.read()

            if len(data) != int(expected_size):
                raise RuntimeError(
                    f"Downloaded size {len(data)} != "
                    f"manifest size {expected_size}"
                )

            return data

        except Exception as exc:
            last = exc

            if attempt == attempts:
                break

    raise RuntimeError(
        f"Failed to fetch {url}"
    ) from last


def cache_path(cache_dir, pdb_id):
    return cache_dir / f"{pdb_id.lower()}.pkl.gz"


def load_cache(path):
    if not path.exists():
        return None

    try:
        with gzip.open(path, "rb") as f:
            payload = pickle.load(f)

        if payload.get("cache_version") != CACHE_VERSION:
            return None

        return payload["entry"]

    except Exception:
        return None


def save_cache(path, entry):
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    payload = {
        "cache_version": CACHE_VERSION,
        "entry": entry,
    }

    fd, tmp_name = tempfile.mkstemp(
        prefix=path.name + ".",
        suffix=".tmp",
        dir=str(path.parent),
    )

    os.close(fd)

    tmp_path = Path(tmp_name)

    try:
        with gzip.open(tmp_path, "wb") as f:
            pickle.dump(
                payload,
                f,
                protocol=pickle.HIGHEST_PROTOCOL,
            )

        tmp_path.replace(path)

    finally:
        if tmp_path.exists():
            tmp_path.unlink()


# ============================================================
# mmCIF parsing
# ============================================================

def residue_sort_key(key):
    label_seq, auth_seq, ins, comp = key

    def numeric(v):
        try:
            return (0, int(v))
        except Exception:
            return (1, str(v))

    return (
        numeric(label_seq),
        numeric(auth_seq),
        str(ins),
        str(comp),
    )


def parse_entry(pdb_id, block):
    """
    Parse model 1 only.

    Stores:
      - deposition-level crystallographic metadata
      - complete model-1 atom signatures
      - per-polymer-chain signatures
      - per-chain N/CA/C backbone coordinates suitable for
        rigid-body alignment when the number of complete
        occupancy-1 residues matches the cleaned retained m.
    """

    asym_cat = block.get_mmcif_category(
        "_struct_asym."
    )

    poly_cat = block.get_mmcif_category(
        "_entity_poly."
    )

    asym_to_entity = {}

    if asym_cat:
        for asym, entity in zip(
            asym_cat.get("id", []),
            asym_cat.get("entity_id", []),
        ):
            asym_to_entity[
                clean(asym)
            ] = clean(entity)

    polymer_entities = set()

    if poly_cat:
        polymer_entities = {
            clean(x)
            for x in poly_cat.get(
                "entity_id",
                [],
            )
        }

    atom = block.get_mmcif_category(
        "_atom_site."
    )

    if not atom:
        raise RuntimeError(
            f"{pdb_id}: no _atom_site category"
        )

    n = len(
        atom.get(
            "Cartn_x",
            [],
        )
    )

    def col(name):
        vals = atom.get(name, [])

        if vals:
            return vals

        return [""] * n

    group = col("group_PDB")
    asym = col("label_asym_id")
    seq = col("label_seq_id")
    auth_seq = col("auth_seq_id")
    ins = col("pdbx_PDB_ins_code")
    comp = col("label_comp_id")
    atom_id = col("label_atom_id")
    alt = col("label_alt_id")
    x = col("Cartn_x")
    y = col("Cartn_y")
    z = col("Cartn_z")
    b = col("B_iso_or_equiv")
    occ = col("occupancy")
    model = col("pdbx_PDB_model_num")

    all_xyz_rows = []
    all_full_rows = []

    polymer_rows = defaultdict(list)

    water_atoms = 0
    nonpolymer_components = set()

    for i in range(n):
        model_id = clean(model[i])

        if model_id not in {"", "1"}:
            continue

        g = clean(group[i])
        chain = clean(asym[i])
        residue = clean(comp[i])
        atom_name = clean(atom_id[i])
        alt_id = clean(alt[i])

        xx = clean(x[i])
        yy = clean(y[i])
        zz = clean(z[i])

        bb = clean(b[i])
        oo = clean(occ[i])

        # Complete-model signatures intentionally ignore
        # chain/residue numbering. Coordinate/chemical content
        # is what matters here.
        model_xyz = (
            g,
            residue,
            atom_name,
            alt_id,
            xx,
            yy,
            zz,
        )

        model_full = (
            *model_xyz,
            bb,
            oo,
        )

        all_xyz_rows.append(
            model_xyz
        )

        all_full_rows.append(
            model_full
        )

        if residue == "HOH":
            water_atoms += 1

        entity = asym_to_entity.get(
            chain,
            "",
        )

        if entity not in polymer_entities:
            if residue and residue != "HOH":
                nonpolymer_components.add(
                    residue
                )

            continue

        polymer_rows[chain].append(
            {
                "label_seq":
                    clean(seq[i]),
                "auth_seq":
                    clean(auth_seq[i]),
                "ins":
                    clean(ins[i]),
                "comp":
                    residue,
                "atom":
                    atom_name,
                "alt":
                    alt_id,
                "x":
                    xx,
                "y":
                    yy,
                "z":
                    zz,
                "b":
                    bb,
                "occ":
                    oo,
            }
        )

    chains = {}

    for chain_id, rows in sorted(
        polymer_rows.items()
    ):
        rows = sorted(
            rows,
            key=lambda r: (
                residue_sort_key(
                    (
                        r["label_seq"],
                        r["auth_seq"],
                        r["ins"],
                        r["comp"],
                    )
                ),
                r["atom"],
                r["alt"],
            ),
        )

        residue_map = defaultdict(list)

        for r in rows:
            key = (
                r["label_seq"],
                r["auth_seq"],
                r["ins"],
                r["comp"],
            )

            residue_map[key].append(r)

        residue_keys = sorted(
            residue_map,
            key=residue_sort_key,
        )

        residue_sequence = [
            key[3]
            for key in residue_keys
        ]

        xyz_rows = []
        full_rows = []

        for r in rows:
            # Chain signatures intentionally ignore residue
            # numbering but preserve residue/atom order.
            xyz_rows.append(
                (
                    r["comp"],
                    r["atom"],
                    r["alt"],
                    r["x"],
                    r["y"],
                    r["z"],
                )
            )

            full_rows.append(
                (
                    r["comp"],
                    r["atom"],
                    r["alt"],
                    r["x"],
                    r["y"],
                    r["z"],
                    r["b"],
                    r["occ"],
                )
            )

        backbone_xyz = []
        backbone_b = []
        backbone_occ = []
        backbone_residue_names = []

        for key in residue_keys:
            rrows = residue_map[key]

            by_atom = defaultdict(list)

            for r in rrows:
                if r["atom"] in {
                    "N",
                    "CA",
                    "C",
                }:
                    by_atom[
                        r["atom"]
                    ].append(r)

            selected = {}

            valid = True

            for atom_name in (
                "N",
                "CA",
                "C",
            ):
                candidates = by_atom.get(
                    atom_name,
                    [],
                )

                # Accepted cleaned chains should generally
                # have a single occupancy-1 backbone atom.
                # Prefer blank altloc, then any occupancy-1
                # candidate.
                candidates = sorted(
                    candidates,
                    key=lambda r: (
                        0 if not r["alt"] else 1,
                        r["alt"],
                    ),
                )

                chosen = None

                for candidate in candidates:
                    occupancy = fnum(
                        candidate["occ"]
                    )

                    if (
                        occupancy is not None
                        and abs(
                            occupancy - 1.0
                        ) <= 1e-9
                    ):
                        chosen = candidate
                        break

                if chosen is None:
                    valid = False
                    break

                selected[
                    atom_name
                ] = chosen

            if not valid:
                continue

            backbone_residue_names.append(
                key[3]
            )

            for atom_name in (
                "N",
                "CA",
                "C",
            ):
                r = selected[
                    atom_name
                ]

                coords = (
                    fnum(r["x"]),
                    fnum(r["y"]),
                    fnum(r["z"]),
                )

                if any(
                    value is None
                    for value in coords
                ):
                    valid = False
                    break

                backbone_xyz.append(
                    coords
                )

                backbone_b.append(
                    fnum(r["b"])
                )

                backbone_occ.append(
                    fnum(r["occ"])
                )

            if not valid:
                # Remove residue if a coordinate conversion
                # failed unexpectedly.
                del backbone_xyz[-3:]
                del backbone_b[-3:]
                del backbone_occ[-3:]
                backbone_residue_names.pop()

        chains[chain_id] = {
            "entity_id":
                asym_to_entity.get(
                    chain_id,
                    "",
                ),
            "residue_count_raw":
                len(residue_keys),
            "atom_count":
                len(rows),
            "sequence":
                residue_sequence,
            "sequence_sha256":
                sha_rows(
                    [
                        tuple(
                            residue_sequence
                        )
                    ]
                ),
            "xyz_sha256":
                sha_rows(
                    xyz_rows
                ),
            "full_sha256":
                sha_rows(
                    full_rows
                ),
            "backbone_complete_residue_count":
                len(
                    backbone_residue_names
                ),
            "backbone_residue_names":
                backbone_residue_names,
            "backbone_xyz":
                backbone_xyz,
            "backbone_b":
                backbone_b,
            "backbone_occ":
                backbone_occ,
        }

    cell = {
        "a":
            first(
                block,
                "_cell.length_a",
            ),
        "b":
            first(
                block,
                "_cell.length_b",
            ),
        "c":
            first(
                block,
                "_cell.length_c",
            ),
        "alpha":
            first(
                block,
                "_cell.angle_alpha",
            ),
        "beta":
            first(
                block,
                "_cell.angle_beta",
            ),
        "gamma":
            first(
                block,
                "_cell.angle_gamma",
            ),
    }

    return {
        "pdb_id":
            pdb_id,
        "experimental_method":
            distinct_join(
                values(
                    block,
                    "_exptl.method",
                )
            ),
        "space_group":
            first(
                block,
                "_space_group.name_H-M_alt",
                "_symmetry.space_group_name_H-M",
            ),
        "cell":
            cell,
        "resolution":
            distinct_join(
                values(
                    block,
                    "_refine.ls_d_res_high",
                )
            ),
        "r_work":
            distinct_join(
                values(
                    block,
                    "_refine.ls_R_factor_R_work",
                )
            ),
        "r_free":
            distinct_join(
                values(
                    block,
                    "_refine.ls_R_factor_R_free",
                )
            ),
        "structure_factor_status":
            first(
                block,
                "_pdbx_database_status.status_code_sf",
            ),
        "model1_atom_count":
            len(
                all_xyz_rows
            ),
        "water_atom_count":
            water_atoms,
        "nonpolymer_components":
            sorted(
                nonpolymer_components
            ),
        "model1_xyz_sha256":
            sha_rows(
                all_xyz_rows
            ),
        "model1_full_sha256":
            sha_rows(
                all_full_rows
            ),
        "chains":
            chains,
    }


def get_entry(
    pdb_id,
    source,
    bucket_url,
    cache_dir,
):
    cp = cache_path(
        cache_dir,
        pdb_id,
    )

    cached = load_cache(cp)

    if cached is not None:
        return cached, True

    src = source[pdb_id]

    url = (
        bucket_url.rstrip("/")
        + "/"
        + str(src["s3_key"])
    )

    compressed = fetch_bytes(
        url,
        src["size_bytes"],
    )

    raw = gzip.decompress(
        compressed
    )

    block = cif.read_string(
        raw
    ).sole_block()

    entry = parse_entry(
        pdb_id,
        block,
    )

    save_cache(
        cp,
        entry,
    )

    return entry, False


# ============================================================
# Multiset comparisons
# ============================================================

def multiset(chain_map, key):
    return Counter(
        chain[key]
        for chain in chain_map.values()
    )


def multiset_intersection_count(a, b):
    return sum(
        (a & b).values()
    )


# ============================================================
# CLI
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Batch evidence analysis for the frozen 2026 "
            "nonzero 1-10 mA BRI detailed-inspection population."
        )
    )

    parser.add_argument(
        "--config",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--review-csv",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--stage11-parquet",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--source-manifest",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--expected-pairs",
        type=int,
        default=None,
    )

    parser.add_argument(
        "--expected-chain-pairs",
        type=int,
        default=None,
    )

    return parser.parse_args()


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    cache_dir = (
        args.output_dir
        / "mmcif_cache"
    )

    cache_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    chain_csv = (
        args.output_dir
        / "chain_pair_evidence.csv"
    )

    chain_parquet = (
        args.output_dir
        / "chain_pair_evidence.parquet"
    )

    pair_csv = (
        args.output_dir
        / "pair_evidence.csv"
    )

    pair_parquet = (
        args.output_dir
        / "pair_evidence.parquet"
    )

    summary_json = (
        args.output_dir
        / "summary.json"
    )

    # --------------------------------------------------------
    # Canonical snapshot URL
    # --------------------------------------------------------

    with args.config.open(
        "r",
        encoding="utf-8",
    ) as f:
        config = yaml.safe_load(f)

    bucket_url = (
        config["snapshot"][
            "bucket_url"
        ].rstrip("/")
    )

    # --------------------------------------------------------
    # Current 2026 review population
    # --------------------------------------------------------

    with args.review_csv.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as f:
        review_rows = list(
            csv.DictReader(f)
        )

    selected = []

    review_by_pair = {}

    for row in review_rows:
        if (
            int(
                row[
                    "exact_bri_chain_pair_count"
                ]
            ) != 0
        ):
            continue

        minimum = int(
            row["minimum_d_bri_mA"]
        )

        if not (
            1 <= minimum <= 10
        ):
            continue

        if (
            row["detailed_review_status"]
            != "unreviewed"
        ):
            continue

        pair = P(
            row["left_pdb_id"],
            row["right_pdb_id"],
        )

        if pair in review_by_pair:
            raise RuntimeError(
                f"Duplicate review pair: {pair}"
            )

        review_by_pair[pair] = row
        selected.append(pair)

    if (
        args.expected_pairs
        is not None
        and len(selected)
        != args.expected_pairs
    ):
        raise RuntimeError(
            f"Expected {args.expected_pairs} pairs, "
            f"found {len(selected)}"
        )

    # --------------------------------------------------------
    # Chain-level Stage-11 evidence
    # --------------------------------------------------------

    stage11 = pq.read_table(
        args.stage11_parquet,
        columns=[
            "query_pdb_id",
            "query_label_chain_id",
            "subject_pdb_id",
            "subject_label_chain_id",
            "retained_residue_count",
            "d_bri_mA",
        ],
    ).to_pylist()

    chain_rows_by_pair = defaultdict(
        list
    )

    for row in stage11:
        d = int(
            row["d_bri_mA"]
        )

        if not (
            1 <= d <= 10
        ):
            continue

        pair = P(
            row["query_pdb_id"],
            row["subject_pdb_id"],
        )

        if pair not in review_by_pair:
            continue

        chain_rows_by_pair[
            pair
        ].append(row)

    missing = (
        set(selected)
        - set(chain_rows_by_pair)
    )

    if missing:
        raise RuntimeError(
            "Missing Stage-11 chain evidence "
            f"for {len(missing)} pairs; "
            f"first={sorted(missing)[:10]}"
        )

    total_chain_pairs = sum(
        len(rows)
        for rows in chain_rows_by_pair.values()
    )

    if (
        args.expected_chain_pairs
        is not None
        and total_chain_pairs
        != args.expected_chain_pairs
    ):
        raise RuntimeError(
            f"Expected {args.expected_chain_pairs} "
            f"chain pairs, found {total_chain_pairs}"
        )

    # --------------------------------------------------------
    # Determine m=1-only versus nontrivial pairs
    # --------------------------------------------------------

    trivial_m1_pairs = set()
    nontrivial_pairs = set()

    for pair, rows in (
        chain_rows_by_pair.items()
    ):
        ms = {
            int(
                row[
                    "retained_residue_count"
                ]
            )
            for row in rows
        }

        if ms == {1}:
            trivial_m1_pairs.add(
                pair
            )
        else:
            nontrivial_pairs.add(
                pair
            )

    required_pdbs = {
        pdb
        for pair in nontrivial_pairs
        for pdb in pair
    }

    # --------------------------------------------------------
    # Source manifest only for nontrivial structures
    # --------------------------------------------------------

    manifest_rows = pq.read_table(
        args.source_manifest,
        columns=[
            "pdb_id",
            "s3_key",
            "size_bytes",
            "etag",
        ],
    ).to_pylist()

    source = {
        str(row["pdb_id"]).lower():
            row
        for row in manifest_rows
        if (
            str(
                row["pdb_id"]
            ).lower()
            in required_pdbs
        )
    }

    unresolved = (
        required_pdbs
        - set(source)
    )

    if unresolved:
        raise RuntimeError(
            "Missing source objects: "
            + repr(
                sorted(
                    unresolved
                )[:20]
            )
        )

    # --------------------------------------------------------
    # Parse each required deposition once
    # --------------------------------------------------------

    entries = {}

    cache_hits = 0
    cache_misses = 0

    print(
        "===== 2026 NONZERO BATCH EVIDENCE =====",
        flush=True,
    )

    print(
        "Deposition pairs:",
        f"{len(selected):,}",
        flush=True,
    )

    print(
        "Chain-pair evidence:",
        f"{total_chain_pairs:,}",
        flush=True,
    )

    print(
        "m=1-only deposition pairs:",
        f"{len(trivial_m1_pairs):,}",
        flush=True,
    )

    print(
        "Nontrivial deposition pairs:",
        f"{len(nontrivial_pairs):,}",
        flush=True,
    )

    print(
        "mmCIF depositions required:",
        f"{len(required_pdbs):,}",
        flush=True,
    )

    for index, pdb_id in enumerate(
        sorted(required_pdbs),
        1,
    ):
        entry, cached = get_entry(
            pdb_id,
            source,
            bucket_url,
            cache_dir,
        )

        entries[pdb_id] = entry

        if cached:
            cache_hits += 1
        else:
            cache_misses += 1

        if (
            index == 1
            or index % 20 == 0
            or index == len(required_pdbs)
        ):
            print(
                f"Parsed/resolved "
                f"{index:,}/"
                f"{len(required_pdbs):,}",
                flush=True,
            )

    # --------------------------------------------------------
    # Chain-pair evidence
    # --------------------------------------------------------

    chain_output = []

    pair_chain_summaries = (
        defaultdict(list)
    )

    for pair in sorted(
        chain_rows_by_pair
    ):
        for row in (
            chain_rows_by_pair[pair]
        ):
            q_pdb = str(
                row[
                    "query_pdb_id"
                ]
            ).lower()

            s_pdb = str(
                row[
                    "subject_pdb_id"
                ]
            ).lower()

            q_chain = str(
                row[
                    "query_label_chain_id"
                ]
            )

            s_chain = str(
                row[
                    "subject_label_chain_id"
                ]
            )

            m = int(
                row[
                    "retained_residue_count"
                ]
            )

            d = int(
                row[
                    "d_bri_mA"
                ]
            )

            record = {
                "left_pdb_id":
                    pair[0],
                "right_pdb_id":
                    pair[1],
                "query_pdb_id":
                    q_pdb,
                "query_label_chain_id":
                    q_chain,
                "subject_pdb_id":
                    s_pdb,
                "subject_label_chain_id":
                    s_chain,
                "retained_residue_count":
                    m,
                "d_bri_mA":
                    d,
                "trivial_m1":
                    m == 1,
                "query_chain_found":
                    None,
                "subject_chain_found":
                    None,
                "raw_xyz_equal":
                    None,
                "raw_full_equal":
                    None,
                "sequence_equal":
                    None,
                "query_backbone_complete_residue_count":
                    None,
                "subject_backbone_complete_residue_count":
                    None,
                "alignment_available":
                    False,
                "alignment_rmsd_A":
                    None,
                "alignment_max_atom_deviation_A":
                    None,
                "alignment_mean_atom_deviation_A":
                    None,
                "backbone_b_equal":
                    None,
                "backbone_b_max_abs_delta":
                    None,
                "backbone_occupancy_equal":
                    None,
            }

            # m=1-only evidence deliberately requires no
            # structure download/alignment.
            if pair in trivial_m1_pairs:
                chain_output.append(
                    record
                )

                pair_chain_summaries[
                    pair
                ].append(record)

                continue

            q_entry = entries[
                q_pdb
            ]

            s_entry = entries[
                s_pdb
            ]

            q = q_entry[
                "chains"
            ].get(q_chain)

            s = s_entry[
                "chains"
            ].get(s_chain)

            record[
                "query_chain_found"
            ] = q is not None

            record[
                "subject_chain_found"
            ] = s is not None

            if (
                q is None
                or s is None
            ):
                chain_output.append(
                    record
                )

                pair_chain_summaries[
                    pair
                ].append(record)

                continue

            record[
                "raw_xyz_equal"
            ] = (
                q["xyz_sha256"]
                == s["xyz_sha256"]
            )

            record[
                "raw_full_equal"
            ] = (
                q["full_sha256"]
                == s["full_sha256"]
            )

            record[
                "sequence_equal"
            ] = (
                q["sequence"]
                == s["sequence"]
            )

            qn = int(
                q[
                    "backbone_complete_residue_count"
                ]
            )

            sn = int(
                s[
                    "backbone_complete_residue_count"
                ]
            )

            record[
                "query_backbone_complete_residue_count"
            ] = qn

            record[
                "subject_backbone_complete_residue_count"
            ] = sn

            # Only align when the raw source reconstruction
            # independently reproduces the cleaned retained m.
            if (
                qn == m
                and sn == m
            ):
                metrics = kabsch_metrics(
                    q[
                        "backbone_xyz"
                    ],
                    s[
                        "backbone_xyz"
                    ],
                )

                record[
                    "alignment_available"
                ] = True

                record[
                    "alignment_rmsd_A"
                ] = metrics[
                    "rmsd_A"
                ]

                record[
                    "alignment_max_atom_deviation_A"
                ] = metrics[
                    "max_atom_deviation_A"
                ]

                record[
                    "alignment_mean_atom_deviation_A"
                ] = metrics[
                    "mean_atom_deviation_A"
                ]

                qb = np.asarray(
                    q[
                        "backbone_b"
                    ],
                    dtype=float,
                )

                sb = np.asarray(
                    s[
                        "backbone_b"
                    ],
                    dtype=float,
                )

                if (
                    qb.shape
                    == sb.shape
                    and np.all(
                        np.isfinite(qb)
                    )
                    and np.all(
                        np.isfinite(sb)
                    )
                ):
                    record[
                        "backbone_b_equal"
                    ] = bool(
                        np.array_equal(
                            qb,
                            sb,
                        )
                    )

                    record[
                        "backbone_b_max_abs_delta"
                    ] = float(
                        np.max(
                            np.abs(
                                qb - sb
                            )
                        )
                    )

                qo = np.asarray(
                    q[
                        "backbone_occ"
                    ],
                    dtype=float,
                )

                so = np.asarray(
                    s[
                        "backbone_occ"
                    ],
                    dtype=float,
                )

                if (
                    qo.shape
                    == so.shape
                    and np.all(
                        np.isfinite(qo)
                    )
                    and np.all(
                        np.isfinite(so)
                    )
                ):
                    record[
                        "backbone_occupancy_equal"
                    ] = bool(
                        np.array_equal(
                            qo,
                            so,
                        )
                    )

            chain_output.append(
                record
            )

            pair_chain_summaries[
                pair
            ].append(record)

    # --------------------------------------------------------
    # Pair-level evidence
    # --------------------------------------------------------

    pair_output = []

    for pair in sorted(selected):
        review = review_by_pair[
            pair
        ]

        chain_records = (
            pair_chain_summaries[
                pair
            ]
        )

        all_m1 = all(
            r["trivial_m1"]
            for r in chain_records
        )

        if all_m1:
            pair_output.append(
                {
                    "deposition_pair_id":
                        review[
                            "deposition_pair_id"
                        ],
                    "left_pdb_id":
                        pair[0],
                    "right_pdb_id":
                        pair[1],
                    "chain_pair_count":
                        len(chain_records),
                    "minimum_d_bri_mA":
                        min(
                            r["d_bri_mA"]
                            for r in chain_records
                        ),
                    "maximum_d_bri_mA":
                        max(
                            r["d_bri_mA"]
                            for r in chain_records
                        ),
                    "retained_lengths":
                        "1",
                    "all_m1":
                        True,
                    "matched_chain_raw_xyz_equal_count":
                        0,
                    "matched_chain_raw_full_equal_count":
                        0,
                    "matched_chain_alignment_available_count":
                        0,
                    "maximum_alignment_rmsd_A":
                        None,
                    "maximum_alignment_atom_deviation_A":
                        None,
                    "left_polymer_chain_count":
                        None,
                    "right_polymer_chain_count":
                        None,
                    "shared_polymer_identity_chain_count":
                        None,
                    "shared_polymer_raw_xyz_chain_count":
                        None,
                    "shared_polymer_full_chain_count":
                        None,
                    "complete_polymer_identity_equal":
                        None,
                    "complete_polymer_raw_xyz_equal":
                        None,
                    "complete_polymer_full_equal":
                        None,
                    "complete_model1_raw_xyz_equal":
                        None,
                    "complete_model1_full_equal":
                        None,
                    "space_group_equal":
                        None,
                    "maximum_cell_length_delta_A":
                        None,
                    "maximum_cell_angle_delta_deg":
                        None,
                    "left_resolution":
                        review[
                            "left_resolution_angstrom"
                        ],
                    "right_resolution":
                        review[
                            "right_resolution_angstrom"
                        ],
                    "left_title":
                        " ".join(
                            review[
                                "left_struct_title"
                            ].split()
                        ),
                    "right_title":
                        " ".join(
                            review[
                                "right_struct_title"
                            ].split()
                        ),
                    "evidence_category":
                        "TRIVIAL_M1_ONLY",
                    "requires_manual_review":
                        False,
                }
            )

            continue

        left = entries[
            pair[0]
        ]

        right = entries[
            pair[1]
        ]

        left_chains = left[
            "chains"
        ]

        right_chains = right[
            "chains"
        ]

        left_identity = multiset(
            left_chains,
            "sequence_sha256",
        )

        right_identity = multiset(
            right_chains,
            "sequence_sha256",
        )

        left_xyz = multiset(
            left_chains,
            "xyz_sha256",
        )

        right_xyz = multiset(
            right_chains,
            "xyz_sha256",
        )

        left_full = multiset(
            left_chains,
            "full_sha256",
        )

        right_full = multiset(
            right_chains,
            "full_sha256",
        )

        shared_identity = (
            multiset_intersection_count(
                left_identity,
                right_identity,
            )
        )

        shared_xyz = (
            multiset_intersection_count(
                left_xyz,
                right_xyz,
            )
        )

        shared_full = (
            multiset_intersection_count(
                left_full,
                right_full,
            )
        )

        complete_identity = (
            left_identity
            == right_identity
        )

        complete_xyz = (
            left_xyz
            == right_xyz
        )

        complete_full = (
            left_full
            == right_full
        )

        complete_model_xyz = (
            left[
                "model1_xyz_sha256"
            ]
            == right[
                "model1_xyz_sha256"
            ]
        )

        complete_model_full = (
            left[
                "model1_full_sha256"
            ]
            == right[
                "model1_full_sha256"
            ]
        )

        raw_xyz_count = sum(
            r["raw_xyz_equal"] is True
            for r in chain_records
        )

        raw_full_count = sum(
            r["raw_full_equal"] is True
            for r in chain_records
        )

        aligned = [
            r
            for r in chain_records
            if r[
                "alignment_available"
            ]
        ]

        max_rmsd = (
            max(
                r[
                    "alignment_rmsd_A"
                ]
                for r in aligned
            )
            if aligned
            else None
        )

        max_atom_dev = (
            max(
                r[
                    "alignment_max_atom_deviation_A"
                ]
                for r in aligned
            )
            if aligned
            else None
        )

        cell_length_deltas = [
            safe_float_delta(
                left["cell"]["a"],
                right["cell"]["a"],
            ),
            safe_float_delta(
                left["cell"]["b"],
                right["cell"]["b"],
            ),
            safe_float_delta(
                left["cell"]["c"],
                right["cell"]["c"],
            ),
        ]

        cell_angle_deltas = [
            safe_float_delta(
                left["cell"]["alpha"],
                right["cell"]["alpha"],
            ),
            safe_float_delta(
                left["cell"]["beta"],
                right["cell"]["beta"],
            ),
            safe_float_delta(
                left["cell"]["gamma"],
                right["cell"]["gamma"],
            ),
        ]

        valid_length_deltas = [
            x
            for x in cell_length_deltas
            if x is not None
        ]

        valid_angle_deltas = [
            x
            for x in cell_angle_deltas
            if x is not None
        ]

        # Evidence categories only.
        # These are NOT scientific duplicate decisions.
        if complete_model_full:
            category = (
                "COMPLETE_MODEL1_RAW_REUSE"
            )
            manual = False

        elif complete_full:
            category = (
                "COMPLETE_POLYMER_RAW_REUSE"
            )
            manual = False

        elif complete_xyz:
            category = (
                "COMPLETE_POLYMER_XYZ_REUSE"
            )
            manual = False

        elif shared_full > 0:
            category = (
                "PARTIAL_POLYMER_RAW_REUSE"
            )
            manual = False

        elif shared_xyz > 0:
            category = (
                "PARTIAL_POLYMER_XYZ_REUSE"
            )
            manual = False

        elif (
            len(aligned)
            == len(chain_records)
        ):
            category = (
                "NEAR_BRI_GEOMETRY_ONLY"
            )
            manual = True

        else:
            category = (
                "AMBIGUOUS_REQUIRES_MANUAL_REVIEW"
            )
            manual = True

        retained_lengths = sorted(
            {
                int(
                    r[
                        "retained_residue_count"
                    ]
                )
                for r in chain_records
            }
        )

        pair_output.append(
            {
                "deposition_pair_id":
                    review[
                        "deposition_pair_id"
                    ],
                "left_pdb_id":
                    pair[0],
                "right_pdb_id":
                    pair[1],
                "chain_pair_count":
                    len(chain_records),
                "minimum_d_bri_mA":
                    min(
                        r["d_bri_mA"]
                        for r in chain_records
                    ),
                "maximum_d_bri_mA":
                    max(
                        r["d_bri_mA"]
                        for r in chain_records
                    ),
                "retained_lengths":
                    ";".join(
                        map(
                            str,
                            retained_lengths,
                        )
                    ),
                "all_m1":
                    False,
                "matched_chain_raw_xyz_equal_count":
                    raw_xyz_count,
                "matched_chain_raw_full_equal_count":
                    raw_full_count,
                "matched_chain_alignment_available_count":
                    len(aligned),
                "maximum_alignment_rmsd_A":
                    max_rmsd,
                "maximum_alignment_atom_deviation_A":
                    max_atom_dev,
                "left_polymer_chain_count":
                    len(
                        left_chains
                    ),
                "right_polymer_chain_count":
                    len(
                        right_chains
                    ),
                "shared_polymer_identity_chain_count":
                    shared_identity,
                "shared_polymer_raw_xyz_chain_count":
                    shared_xyz,
                "shared_polymer_full_chain_count":
                    shared_full,
                "complete_polymer_identity_equal":
                    complete_identity,
                "complete_polymer_raw_xyz_equal":
                    complete_xyz,
                "complete_polymer_full_equal":
                    complete_full,
                "complete_model1_raw_xyz_equal":
                    complete_model_xyz,
                "complete_model1_full_equal":
                    complete_model_full,
                "space_group_equal":
                    (
                        left[
                            "space_group"
                        ]
                        == right[
                            "space_group"
                        ]
                    ),
                "maximum_cell_length_delta_A":
                    (
                        max(
                            valid_length_deltas
                        )
                        if (
                            len(
                                valid_length_deltas
                            )
                            == 3
                        )
                        else None
                    ),
                "maximum_cell_angle_delta_deg":
                    (
                        max(
                            valid_angle_deltas
                        )
                        if (
                            len(
                                valid_angle_deltas
                            )
                            == 3
                        )
                        else None
                    ),
                "left_resolution":
                    left[
                        "resolution"
                    ],
                "right_resolution":
                    right[
                        "resolution"
                    ],
                "left_title":
                    " ".join(
                        review[
                            "left_struct_title"
                        ].split()
                    ),
                "right_title":
                    " ".join(
                        review[
                            "right_struct_title"
                        ].split()
                    ),
                "evidence_category":
                    category,
                "requires_manual_review":
                    manual,
            }
        )

    # --------------------------------------------------------
    # Write chain evidence
    # --------------------------------------------------------

    chain_fields = list(
        chain_output[0]
    )

    with chain_csv.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=chain_fields,
        )

        writer.writeheader()
        writer.writerows(
            chain_output
        )

    pq.write_table(
        pa.Table.from_pylist(
            chain_output
        ),
        chain_parquet,
        compression="zstd",
    )

    # --------------------------------------------------------
    # Write pair evidence
    # --------------------------------------------------------

    pair_fields = list(
        pair_output[0]
    )

    with pair_csv.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=pair_fields,
        )

        writer.writeheader()
        writer.writerows(
            pair_output
        )

    pq.write_table(
        pa.Table.from_pylist(
            pair_output
        ),
        pair_parquet,
        compression="zstd",
    )

    # --------------------------------------------------------
    # Final audit and summary
    # --------------------------------------------------------

    if len(pair_output) != len(selected):
        raise RuntimeError(
            "Pair output count mismatch"
        )

    if len(chain_output) != total_chain_pairs:
        raise RuntimeError(
            "Chain output count mismatch"
        )

    category_counts = Counter(
        row[
            "evidence_category"
        ]
        for row in pair_output
    )

    manual_count = sum(
        bool(
            row[
                "requires_manual_review"
            ]
        )
        for row in pair_output
    )

    alignment_available = sum(
        row[
            "alignment_available"
        ] is True
        for row in chain_output
    )

    alignment_missing = sum(
        (
            not row["trivial_m1"]
            and row[
                "alignment_available"
            ] is not True
        )
        for row in chain_output
    )

    summary = {
        "cache_version":
            CACHE_VERSION,
        "deposition_pair_count":
            len(
                pair_output
            ),
        "chain_pair_count":
            len(
                chain_output
            ),
        "m1_only_pair_count":
            len(
                trivial_m1_pairs
            ),
        "nontrivial_pair_count":
            len(
                nontrivial_pairs
            ),
        "required_mmcif_deposition_count":
            len(
                required_pdbs
            ),
        "cache_hits":
            cache_hits,
        "cache_misses":
            cache_misses,
        "alignment_available_chain_pair_count":
            alignment_available,
        "alignment_missing_nontrivial_chain_pair_count":
            alignment_missing,
        "evidence_category_counts":
            dict(
                sorted(
                    category_counts.items()
                )
            ),
        "manual_review_pair_count":
            manual_count,
        "scientific_duplicate_decisions_made":
            0,
        "review_csv_modified":
            False,
        "old_snapshot_comparison_used":
            False,
    }

    with summary_json.open(
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(
            summary,
            f,
            indent=2,
            sort_keys=True,
        )

        f.write("\n")

    print()
    print(
        "===== EVIDENCE CATEGORY COUNTS ====="
    )

    for category, count in sorted(
        category_counts.items()
    ):
        print(
            f"{category}: {count:,}"
        )

    print()
    print(
        "Pairs requiring manual review:",
        f"{manual_count:,}",
    )

    print(
        "Pairs handled as m=1-only evidence:",
        f"{len(trivial_m1_pairs):,}",
    )

    print(
        "Nontrivial chain alignments available:",
        f"{alignment_available:,}",
    )

    print(
        "Nontrivial chain alignments unavailable:",
        f"{alignment_missing:,}",
    )

    print()
    print(
        "Chain evidence:",
        chain_parquet,
    )

    print(
        "Pair evidence:",
        pair_parquet,
    )

    print(
        "Summary:",
        summary_json,
    )

    print()
    print(
        "Scientific duplicate decisions made: 0"
    )

    print(
        "Review CSV modified: NO"
    )

    print(
        "Old-snapshot comparison used: NO"
    )

    print(
        "Evidence source: frozen 2026 snapshot only"
    )

    print(
        "2026 NONZERO BATCH EVIDENCE: PASS"
    )


if __name__ == "__main__":
    main()
