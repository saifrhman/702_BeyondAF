"""Sequence-redundancy resolution over the geometric survivors.

Geometry runs first and is not recomputed: this stage consumes the retained
chains of a completed Gold release and asks a separate question -- which of
those chains are redundant *by sequence*.

Representative selection is deliberately NOT delegated to MMseqs2. Its
``--cluster-mode 0`` set cover is a greedy, order-dependent heuristic; on this
very population it was observed to place byte-identical sequences under
different representatives (for example ``2n0k_A`` and ``2n0k_B``, the same
89-residue sequence, landed under ``2wj7_E`` and ``4jut_C``). MMseqs2 is
therefore used for one thing only -- deciding which chains belong together --
and the survivor within each cluster is chosen by the project's own
deterministic ranking, the same rule Stage 14b applies to geometric components:

    1. untrimmed preferred over terminal-trimmed
    2. fewer defective (dirty) residues preferred
    3. better nominal resolution preferred, but only where the cluster's
       experimental methods are comparable
    4. canonical chain key, as a total deterministic tie-break

Rule 3 carries a caveat worth stating plainly: the resolution term is only
consulted when every member of a cluster shares one experimental-method tuple,
and deposition metadata exists for a minority of entries (the
``downstream_metadata`` stage only fetches the depositions that participate in
*geometric* near-duplicate pairs). Where metadata is absent the method tuples
cannot agree, so rule 3 stands down and selection falls through to rules 1, 2
and 4 -- which are defined for every chain and are themselves deterministic.

The ranking helpers are imported from the Stage 14b entry point rather than
reimplemented, so the two stages cannot drift apart in what they consider a
good representative.
"""

from __future__ import annotations

import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

#: The ranking, in the vocabulary the resolved configuration records.
RANKING: tuple[str, ...] = (
    "terminal_trimmed_false_preferred",
    "lower_dirty_residue_count_preferred",
    "comparable_method_resolution_preferred",
    "canonical_chain_key_tiebreak",
)

POLICY_NAME = "comp702_sequence_representative_selection"
POLICY_VERSION = "1.0"

_STAGE14_ENTRY_POINT = "scripts/select_stage14_representatives.py"


class SequenceClusteringError(RuntimeError):
    """Raised when sequence-redundancy resolution cannot proceed safely."""


def load_stage14_policy(repo_root: Path):
    """Import the Stage 14b ranking helpers from their entry point.

    Loading the script rather than copying its functions is deliberate: the
    requirement is that geometric and sequence representative selection agree,
    and the only way to guarantee that is to execute the same code. The script
    guards its ``main()`` behind ``if __name__ == "__main__"``, so importing it
    runs nothing.
    """

    path = repo_root / _STAGE14_ENTRY_POINT

    if not path.is_file():
        raise SequenceClusteringError(
            f"Stage 14b entry point not found at {path}; the sequence stage "
            "reuses its ranking helpers and will not substitute a copy."
        )

    spec = importlib.util.spec_from_file_location(
        "pdbclean._stage14_policy", path
    )

    if spec is None or spec.loader is None:
        raise SequenceClusteringError(f"Cannot load {path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)

    for required in ("methods_tuple", "method_resolution"):
        if not hasattr(module, required):
            raise SequenceClusteringError(
                f"Stage 14b entry point has no {required}(); the ranking "
                "contract has changed and this stage must be revisited."
            )

    return module


@dataclass(frozen=True)
class ChainQuality:
    """The per-chain facts the ranking consumes."""

    terminal_trimmed: bool
    dirty_residue_count: int
    pdb_id: str


def parse_cluster_tsv(path: Path) -> dict[str, list[str]]:
    """Read an MMseqs2 ``*_cluster.tsv`` into {mmseqs_rep: [members]}.

    The MMseqs2 representative is retained only as the cluster's identity; it
    carries no authority over which chain survives.
    """

    clusters: dict[str, list[str]] = {}

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue

            rep, _, member = line.rstrip("\n").partition("\t")

            if not member:
                raise SequenceClusteringError(
                    f"Malformed cluster line in {path}: {line!r}"
                )

            clusters.setdefault(rep, []).append(member)

    return clusters


def make_rank(
    policy,
    members: Sequence[str],
    quality: dict[str, ChainQuality],
    metadata: dict[str, dict[str, Any]],
) -> tuple[Callable[[str], tuple], bool, tuple | None]:
    """Build the Stage 14b ranking closure for one cluster.

    Returns ``(rank, comparable_method, common_method)`` exactly as Stage 14b's
    ``component_ranker`` does.
    """

    method_tuples = {
        policy.methods_tuple(
            (metadata.get(quality[key].pdb_id) or {}).get("experimental_methods")
        )
        for key in members
    }

    comparable_method = len(method_tuples) == 1
    common_method = next(iter(method_tuples)) if comparable_method else None

    # A cluster whose members all lack metadata shares the empty method tuple,
    # which is "comparable" but carries no resolution. Treat it as
    # incomparable so the empty case cannot masquerade as an agreed method.
    if comparable_method and not common_method:
        comparable_method = False
        common_method = None

    def rank(key: str) -> tuple:
        chain = quality[key]

        dirty = chain.dirty_residue_count
        dirty = 10**9 if dirty is None else int(dirty)

        resolution = None

        if comparable_method:
            meta = metadata.get(chain.pdb_id)

            if meta is not None:
                resolution = policy.method_resolution(meta, common_method)

        resolution_missing = resolution is None
        resolution_value = float("inf") if resolution is None else float(resolution)

        return (
            1 if chain.terminal_trimmed else 0,
            dirty,
            1 if resolution_missing else 0,
            resolution_value,
            key,
        )

    return rank, comparable_method, common_method


def choose_representatives(
    clusters: dict[str, list[str]],
    quality: dict[str, ChainQuality],
    metadata: dict[str, dict[str, Any]],
    policy,
) -> tuple[dict[str, str], dict[str, str], int]:
    """Pick one survivor per cluster under the project ranking.

    Returns ``(representative_of_cluster, representative_of_chain, n_resolution_used)``
    where ``representative_of_chain`` maps every input chain -- survivor and
    removed alike -- to the chain that represents it.
    """

    representative_of_cluster: dict[str, str] = {}
    representative_of_chain: dict[str, str] = {}
    resolution_used = 0

    # Sorting the cluster ids makes the walk itself order-independent, so the
    # output cannot depend on the order MMseqs2 happened to emit.
    for cluster_id in sorted(clusters):
        members = sorted(clusters[cluster_id])

        missing = [key for key in members if key not in quality]

        if missing:
            raise SequenceClusteringError(
                f"Cluster {cluster_id} contains chains absent from the input "
                f"population: {missing[:5]}"
            )

        rank, comparable, _ = make_rank(policy, members, quality, metadata)

        if comparable:
            resolution_used += 1

        winner = min(members, key=rank)

        representative_of_cluster[cluster_id] = winner

        for member in members:
            representative_of_chain[member] = winner

    return representative_of_cluster, representative_of_chain, resolution_used


def reconcile(
    *,
    input_chains: Iterable[str],
    retained: Iterable[str],
    removed: Iterable[str],
    representative_of_chain: dict[str, str],
) -> dict[str, Any]:
    """Run the stage's accounting and attribution gates.

    Raises on any failure; returns the numbers it verified.
    """

    input_set = set(input_chains)
    retained_set = set(retained)
    removed_set = set(removed)

    if retained_set & removed_set:
        raise SequenceClusteringError(
            f"{len(retained_set & removed_set)} chains are both retained and "
            "removed"
        )

    if retained_set | removed_set != input_set:
        missing = len(input_set - (retained_set | removed_set))
        extra = len((retained_set | removed_set) - input_set)
        raise SequenceClusteringError(
            f"Chain accounting does not reconcile: {missing} input chains "
            f"unaccounted for, {extra} outputs not in the input"
        )

    if len(retained_set) + len(removed_set) != len(input_set):
        raise SequenceClusteringError("Chain counts do not sum to the input")

    for chain in removed_set:
        rep = representative_of_chain.get(chain)

        if rep is None:
            raise SequenceClusteringError(
                f"Removed chain {chain} has no representative"
            )

        if rep == chain:
            raise SequenceClusteringError(
                f"Removed chain {chain} is its own representative"
            )

        if rep not in retained_set:
            raise SequenceClusteringError(
                f"Removed chain {chain} points at representative {rep}, which "
                "is not retained"
            )

    for chain in retained_set:
        if representative_of_chain.get(chain) != chain:
            raise SequenceClusteringError(
                f"Retained chain {chain} is not its own representative"
            )

    return {
        "input_chain_count": len(input_set),
        "retained_chain_count": len(retained_set),
        "removed_chain_count": len(removed_set),
        "every_removal_attributed": True,
        "no_representative_removed": True,
        "accounting_reconciles": True,
    }
