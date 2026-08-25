"""The alignment index is a bridge between two namespaces, and bridges rot quietly.

OpenFold addresses alignments by chain; the MSA corpus is addressed by sequence
digest.  The index joins them.  Every failure mode of that join is silent at
run time -- a chain resolves to *an* MSA either way, and training proceeds on
the wrong evolutionary signal without raising anything.

These tests therefore do two things the build itself cannot.  They pin the
index format against OpenFold's read contract (``open(join(dir, db))``,
``seek(start)``, ``read(size)``), including the claim that addressing a whole
file with ``start=0`` is equivalent to addressing a slice of a packed database.
And they corrupt a known-good index in each of the ways it could plausibly be
wrong, asserting the verifier fails -- because a verifier that cannot fail is
not evidence.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("pyarrow")

import pyarrow as pa
import pyarrow.parquet as pq


REPO_ROOT = Path(__file__).resolve().parents[2]

BUILDER = REPO_ROOT / "scripts/openfold_training/build_alignment_index.py"
VERIFIER = REPO_ROOT / "scripts/openfold_training/verify_alignment_index.py"


def load_script(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


builder = load_script(BUILDER, "build_alignment_index")
verifier = load_script(VERIFIER, "verify_alignment_index")


# --------------------------------------------------------------------------
# a miniature corpus with the same shape as the real one
# --------------------------------------------------------------------------

SEQUENCES = {
    "alpha": "MVLSEGEWQLVLHVWAKVEAD",
    "beta": "MNIFEMLRIDEGLRLKIYKDT",
    "gamma": "MKTAYIAKQRQISFVKSHFSR",
}

# beta is shared by three chains, exactly as real duplicates are.
CHAINS = {
    "101m_A": "alpha",
    "102l_A": "beta",
    "102l_B": "beta",
    "9zzz_c": "beta",
    "3abc_A": "gamma",
}


def digest(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def a3m_text(sequence: str, depth: int = 3) -> str:
    """A minimal but structurally real a3m: query first, then homologues."""

    lines = [f">query_{digest(sequence)[:8]}", sequence]

    for i in range(depth):
        # Lower case marks an insertion relative to the query, as in a real a3m.
        lines += [f">homologue_{i}", sequence[:-1] + "a"]

    return "\n".join(lines) + "\n"


@pytest.fixture
def corpus(tmp_path: Path):
    """Build a chain map and an MSA store, and return their paths."""

    store = tmp_path / "msas"
    store.mkdir()

    for sequence in SEQUENCES.values():
        (store / f"{digest(sequence)}.a3m").write_text(a3m_text(sequence))

    chains = list(CHAINS)
    shas = [digest(SEQUENCES[CHAINS[c]]) for c in chains]

    chain_map = tmp_path / "chain_to_sequence.parquet"
    pq.write_table(
        pa.table({"openfold_chain_id": chains, "sequence_sha256": shas}),
        chain_map,
    )

    training = tmp_path / "training_chains.parquet"
    pq.write_table(
        pa.table(
            {
                "openfold_chain_id": chains,
                "retained_sequence": [SEQUENCES[CHAINS[c]] for c in chains],
            }
        ),
        training,
    )

    return {
        "root": tmp_path,
        "store": store,
        "chain_map": chain_map,
        "training": training,
    }


def run_builder(corpus, **extra) -> subprocess.CompletedProcess:
    output = corpus["root"] / "out.index"

    argv = [
        sys.executable,
        str(BUILDER),
        "--chain-map", str(corpus["chain_map"]),
        "--msa-store", str(corpus["store"]),
        "--output", str(output),
    ]

    for key, value in extra.items():
        argv += [f"--{key.replace('_', '-')}", str(value)]

    return subprocess.run(argv, capture_output=True, text=True)


def run_verifier(corpus, index: Path, **extra) -> subprocess.CompletedProcess:
    argv = [
        sys.executable,
        str(VERIFIER),
        "--index", str(index),
        "--alignment-dir", str(corpus["store"]),
        "--training-chains", str(corpus["training"]),
        "--verify-chains", "-1",
    ]

    for key, value in extra.items():
        argv += [f"--{key.replace('_', '-')}", str(value)]

    return subprocess.run(argv, capture_output=True, text=True)


@pytest.fixture
def built(corpus):
    result = run_builder(corpus)
    assert result.returncode == 0, result.stderr

    index_path = corpus["root"] / "out.index"
    return corpus, index_path, json.loads(index_path.read_text())


# --------------------------------------------------------------------------
# the format OpenFold requires
# --------------------------------------------------------------------------

def test_every_chain_is_indexed(built):
    _, _, index = built

    assert set(index) == set(CHAINS)


def test_entries_match_openfolds_read_contract(built):
    _, _, index = built

    for entry in index.values():
        assert set(entry) == {"db", "files"}
        assert isinstance(entry["db"], str)
        assert len(entry["files"]) == 1

        name, start, size = entry["files"][0]

        assert name.endswith(".a3m")
        assert start == 0
        assert isinstance(size, int) and size > 0


def test_chain_keys_split_into_file_id_and_chain_id(built):
    """OpenFold does name.rsplit('_', 1) to find the structure file."""

    _, _, index = built

    for chain in index:
        file_id, chain_id = chain.rsplit("_", 1)

        assert file_id
        assert chain_id


def test_reading_an_entry_yields_that_chains_sequence(built):
    corpus, _, index = built

    for chain, key in CHAINS.items():
        entry = index[chain]
        _, start, size = entry["files"][0]

        observed = verifier.query_sequence(
            corpus["store"] / entry["db"], start, size
        )

        assert observed == SEQUENCES[key]


# --------------------------------------------------------------------------
# sequence reuse must cost keys, not bytes
# --------------------------------------------------------------------------

def test_chains_sharing_a_sequence_share_one_entry(built):
    _, _, index = built

    beta = [c for c, k in CHAINS.items() if k == "beta"]

    assert len(beta) == 3

    entries = [index[c] for c in beta]

    assert all(e == entries[0] for e in entries)


def test_distinct_msas_equals_distinct_sequences(built):
    _, _, index = built

    assert len({e["db"] for e in index.values()}) == len(SEQUENCES)


def test_no_msa_bytes_are_duplicated(built):
    """The whole point: 5 chains, 3 files, nothing copied."""

    corpus, _, index = built

    stored = sorted(p.name for p in corpus["store"].glob("*.a3m"))

    assert stored == sorted({e["db"] for e in index.values()})


# --------------------------------------------------------------------------
# the equivalence the builder's design rests on
# --------------------------------------------------------------------------

def test_whole_file_addressing_equals_packed_slice_addressing(corpus):
    """A zero-copy entry and a packed-db entry must read the same bytes.

    This is the claim that lets the index skip duplicating 161 GiB.  It is
    asserted at the byte level, which is the level OpenFold reads at.
    """

    sequence = SEQUENCES["alpha"]
    source = corpus["store"] / f"{digest(sequence)}.a3m"
    raw = source.read_bytes()

    packed = corpus["root"] / "packed_0.db"
    padding = b"\x00" * 4096
    packed.write_bytes(padding + raw + b"trailing junk")

    with packed.open("rb") as handle:
        handle.seek(len(padding))
        from_packed = handle.read(len(raw))

    with source.open("rb") as handle:
        handle.seek(0)
        from_zero_copy = handle.read(len(raw))

    assert from_packed == from_zero_copy == raw


def test_size_must_be_the_whole_file_not_a_prefix(built):
    """A truncated size silently drops homologues from the MSA."""

    corpus, _, index = built

    for entry in index.values():
        _, _, size = entry["files"][0]

        assert size == (corpus["store"] / entry["db"]).stat().st_size


# --------------------------------------------------------------------------
# the builder must refuse to produce a broken index
# --------------------------------------------------------------------------

def test_builder_fails_when_an_msa_is_missing(corpus):
    (corpus["store"] / f"{digest(SEQUENCES['gamma'])}.a3m").unlink()

    result = run_builder(corpus)

    assert result.returncode == 1
    assert "no MSA" in result.stderr


def test_builder_fails_on_an_empty_msa(corpus):
    (corpus["store"] / f"{digest(SEQUENCES['gamma'])}.a3m").write_text("")

    result = run_builder(corpus)

    assert result.returncode == 1
    assert "empty" in result.stderr


def test_builder_enforces_the_expected_chain_count(corpus):
    result = run_builder(corpus, expected_chains=len(CHAINS) + 1)

    assert result.returncode == 1
    assert "expected" in result.stderr


def test_builder_enforces_the_expected_sequence_count(corpus):
    result = run_builder(corpus, expected_sequences=len(SEQUENCES) + 1)

    assert result.returncode == 1
    assert "expected" in result.stderr


def test_builder_accepts_the_correct_expectations(corpus):
    result = run_builder(
        corpus,
        expected_chains=len(CHAINS),
        expected_sequences=len(SEQUENCES),
    )

    assert result.returncode == 0, result.stderr


# --------------------------------------------------------------------------
# the verifier must actually be able to fail
# --------------------------------------------------------------------------

def test_verifier_passes_on_a_good_index(built):
    corpus, index_path, _ = built

    result = run_verifier(corpus, index_path)

    assert result.returncode == 0, result.stdout + result.stderr
    assert "STATUS: PASS" in result.stdout


def test_verifier_catches_a_chain_pointed_at_the_wrong_msa(built):
    """The failure no count-based check can see."""

    corpus, index_path, index = built

    index["101m_A"] = index["102l_A"]          # alpha chain, beta MSA
    index_path.write_text(json.dumps(index))

    result = run_verifier(corpus, index_path)

    assert result.returncode == 1
    assert "not that chain's sequence" in result.stdout


def test_verifier_catches_a_missing_chain(built):
    corpus, index_path, index = built

    del index["3abc_A"]
    index_path.write_text(json.dumps(index))

    result = run_verifier(corpus, index_path)

    assert result.returncode == 1
    assert "no index entry" in result.stdout


def test_verifier_catches_an_entry_for_an_unknown_chain(built):
    corpus, index_path, index = built

    index["9xyz_Z"] = index["101m_A"]
    index_path.write_text(json.dumps(index))

    result = run_verifier(corpus, index_path)

    assert result.returncode == 1
    assert "no training chain" in result.stdout


def test_verifier_catches_a_truncated_byte_range(built):
    """Truncation still parses -- it just yields a shorter query row."""

    corpus, index_path, index = built

    name, start, size = index["101m_A"]["files"][0]
    index["101m_A"] = {"db": index["101m_A"]["db"], "files": [[name, start, 12]]}
    index_path.write_text(json.dumps(index))

    result = run_verifier(corpus, index_path)

    assert result.returncode == 1


def test_verifier_catches_a_dangling_db_reference(built):
    corpus, index_path, index = built

    index["101m_A"] = {
        "db": "does_not_exist.a3m",
        "files": index["101m_A"]["files"],
    }
    index_path.write_text(json.dumps(index))

    result = run_verifier(corpus, index_path)

    assert result.returncode == 1
    assert "could not be read" in result.stdout


def test_verifier_catches_chains_missing_from_the_train_filter(built):
    corpus, index_path, index = built

    filter_path = corpus["root"] / "train_filter.txt"
    filter_path.write_text("\n".join(list(CHAINS) + ["4zzz_Q"]) + "\n")

    result = run_verifier(corpus, index_path, train_filter=filter_path)

    assert result.returncode == 1
    assert "train filter" in result.stdout


def test_verifier_accepts_a_consistent_train_filter(built):
    corpus, index_path, _ = built

    filter_path = corpus["root"] / "train_filter.txt"
    filter_path.write_text("\n".join(CHAINS) + "\n")

    result = run_verifier(corpus, index_path, train_filter=filter_path)

    assert result.returncode == 0, result.stdout + result.stderr


# --------------------------------------------------------------------------
# provenance
# --------------------------------------------------------------------------

def test_manifest_records_what_was_built(corpus):
    manifest = corpus["root"] / "manifest.json"

    result = run_builder(corpus, manifest=manifest)

    assert result.returncode == 0, result.stderr

    recorded = json.loads(manifest.read_text())

    assert recorded["chains_indexed"] == len(CHAINS)
    assert recorded["distinct_msas_referenced"] == len(SEQUENCES)
    assert recorded["chains_reusing_a_shared_msa"] == len(CHAINS) - len(SEQUENCES)
    assert recorded["alignment_filename"] == builder.ALIGNMENT_FILENAME
    assert len(recorded["alignment_index_sha256"]) == 64


def test_manifest_digest_matches_the_index_on_disk(corpus):
    manifest = corpus["root"] / "manifest.json"

    run_builder(corpus, manifest=manifest)

    recorded = json.loads(manifest.read_text())
    index_path = corpus["root"] / "out.index"

    assert recorded["alignment_index_sha256"] == hashlib.sha256(
        index_path.read_bytes()
    ).hexdigest()
