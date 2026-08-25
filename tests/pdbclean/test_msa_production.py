"""Production MSAs must be produced by the protocol the pilot validated.

The pilot is the control: it proves the mmseqs command sequence on nine
representative chains, including a deliberately duplicated sequence pair, and
it is the only evidence that the protocol is correct.  Production inherits that
evidence *only* while it runs the identical protocol.  A parameter changed in
one and not the other silently invalidates the control, so these tests compare
the two scripts directly.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]

PILOT = REPO_ROOT / "scripts/openfold_training/run_msa_pilot.sbatch"
PRODUCTION = REPO_ROOT / "scripts/openfold_training/run_msa_production.sbatch"


def mmseqs_calls(path: Path) -> list[tuple[str, tuple[str, ...]]]:
    """Every mmseqs invocation in a script, as (verb, flags).

    Positional arguments are deliberately excluded: the two scripts read
    different query sets and write to different places, which is the whole
    point of one being production.  What must not differ is the verb sequence
    and every scientific parameter.
    """

    source = path.read_text(encoding="utf-8")

    calls: list[tuple[str, tuple[str, ...]]] = []

    for verb, body in re.findall(
        r'"\$MMSEQS" (\w+)((?:[^\n]*\\\n)*[^\n]*)', source
    ):
        flags = re.findall(
            r'(--[a-z0-9-]+|(?<!\w)-[ae](?!\w))\s+("?[^\s"\\]*"?)?', body
        )

        calls.append(
            (verb, tuple(f"{k} {v}".strip() for k, v in flags if k))
        )

    return calls


@pytest.fixture(scope="module")
def pilot_calls():
    return mmseqs_calls(PILOT)


@pytest.fixture(scope="module")
def production_calls():
    return mmseqs_calls(PRODUCTION)


def test_both_scripts_exist():
    assert PILOT.is_file()
    assert PRODUCTION.is_file()


def test_the_command_sequence_is_identical(pilot_calls, production_calls):
    assert [verb for verb, _ in pilot_calls] == [
        verb for verb, _ in production_calls
    ]


def test_every_scientific_parameter_is_identical(
    pilot_calls, production_calls
):
    for (pilot_verb, pilot_flags), (prod_verb, prod_flags) in zip(
        pilot_calls, production_calls
    ):
        assert pilot_verb == prod_verb
        assert pilot_flags == prod_flags, (
            f"{pilot_verb} differs between pilot and production:\n"
            f"  pilot:      {pilot_flags}\n"
            f"  production: {prod_flags}"
        )


def test_the_search_keeps_three_profile_iterations(production_calls):
    """The single parameter most likely to be 'optimised' into wrongness."""

    search = next(
        flags for verb, flags in production_calls if verb == "search"
    )

    assert "--num-iterations 3" in search


def test_production_installs_the_result_lookup_before_unpacking():
    """The defect the pilot caught, and the reason it must not recur.

    Without a lookup on the *result* database, --unpack-name-mode 1 falls back
    to numeric DB keys and every MSA loses the identity of the sequence it
    describes.
    """

    source = PRODUCTION.read_text(encoding="utf-8")

    lookup = source.index('cp "$WORK/qdb.lookup" "$WORK/uniref.a3m.lookup"')
    unpack = source.index('"$MMSEQS" unpackdb')

    assert lookup < unpack, "the lookup must be installed before unpacking"
    assert "--unpack-name-mode 1" in source


def test_production_refuses_mismatched_keys():
    """Copying qdb.lookup is only correct while the keys correspond."""

    source = PRODUCTION.read_text(encoding="utf-8")

    assert 'RESULT_KEYS=' in source
    assert 'QUERY_KEYS=' in source
    assert "result database keys do not match" in source


def test_production_verifies_every_output_is_traceable(pilot_calls):
    """No orphan outputs, no missing MSAs -- both pilot failure modes."""

    source = PRODUCTION.read_text(encoding="utf-8")

    assert "orphan output not traceable to a query" in source
    assert "missing MSA for query" in source


def test_production_publishes_only_after_verification():
    """A killed job must not leave a half-written MSA in the shared store."""

    source = PRODUCTION.read_text(encoding="utf-8")

    verify = source.index("identity check:")
    publish = source.index('mv -f "$UNPACK"/*.a3m')

    assert verify < publish


def test_production_is_resumable():
    source = PRODUCTION.read_text(encoding="utf-8")

    assert "_SHARD_SUCCESS" in source
    assert "Shard already complete" in source


def test_production_refuses_an_unverified_database():
    source = PRODUCTION.read_text(encoding="utf-8")

    assert "_UNIREF30_DB_SUCCESS" in source
    assert "unverified database" in source


def test_sharding_is_disjoint_and_covers_everything():
    """Every sequence is searched exactly once across the array.

    Re-implements the shard predicate the script uses, because an overlap
    would duplicate expensive searches and a gap would silently omit chains
    from the training set.
    """

    source = PRODUCTION.read_text(encoding="utf-8")

    assert "(index % shard_count) == shard_id" in source

    total, shards = 142_056, 64

    seen: list[int] = []

    for shard in range(shards):
        seen.extend(i for i in range(total) if i % shards == shard)

    assert len(seen) == total
    assert len(set(seen)) == total
