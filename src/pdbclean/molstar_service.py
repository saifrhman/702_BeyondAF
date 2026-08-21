"""On-demand Mol* scene generation for detected duplicate pairs.

The Duplicate Explorer must be able to visualise *any* detected pair whose
source structures are actually available, not only the handful of pairs a
preparation script happened to pre-generate.

This module resolves the two structures for a pair, builds a scene with the
validated logic in :mod:`pdbclean.molstar_scenes`, and caches the result.

Three rules govern everything here:

**Mol\\* is inspection only.**  A scene is built *from* the recorded
complete-BRI classification; it never produces or revises one.  Nothing in this
module writes to a stage output, a release, or run provenance.

**Snapshot correctness.**  A pair detected in a run is visualised with the
structures belonging to *that run's resolved snapshot*.  Fetching whatever the
PDB serves today would silently show a different structure than the one the
science was computed on, so remote fetching is not implemented: a structure is
either available from a preserved/materialised local source or the pair is
reported unavailable, with the reason.

**The cache is disposable.**  Generated scenes are derived data keyed by run,
snapshot, pair and view.  Deleting the cache loses nothing scientific, and
generating a scene never changes any scientific hash.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import gemmi

from pdbclean import molstar_scenes as scenes


SCENE_SCHEMA_VERSION = "1.0"

#: Colours used by the frozen example scenes, kept for visual consistency.
COLOUR_A = "#2563eb"
COLOUR_B = "#f97316"


# ---------------------------------------------------------------------------
# Availability
# ---------------------------------------------------------------------------

AVAILABLE = "available"
CACHED = "cached"

#: The source structure is not held locally for this run's snapshot.
STRUCTURE_UNAVAILABLE = "structure_unavailable"

#: The run's snapshot has not been preserved or materialised, so its
#: structures cannot be resolved without risking a different PDB version.
SNAPSHOT_NOT_MATERIALISED = "snapshot_not_materialised"

#: The chain named by the pair is not present in the resolved structure.
CHAIN_NOT_FOUND = "chain_not_found"

#: The two chains do not share a comparable backbone, so a superposed view
#: cannot be built (a side-by-side view usually still can).
BACKBONE_MISMATCH = "backbone_mismatch"

#: The source object identity is known from the Bronze manifest but no local
#: copy exists yet. Recoverable: it can be materialised on demand.
SOURCE_NOT_MATERIALISED = "source_not_materialised"

#: The run's Bronze source manifest is absent, so the source object identity
#: for this snapshot cannot be established at all.
SOURCE_MANIFEST_MISSING = "source_manifest_missing"

#: The entry is not listed in this snapshot's Bronze manifest.
NOT_IN_SNAPSHOT = "not_in_snapshot"

#: Retrieval or ETag verification failed.
SOURCE_FETCH_FAILED = "source_fetch_failed"

UNSUPPORTED_SOURCE = "unsupported_historical_source"


REASON_TEXT = {
    STRUCTURE_UNAVAILABLE: (
        "Source structure not available locally for this run's snapshot. "
        "Preserve or materialise the snapshot to enable inspection."
    ),
    SNAPSHOT_NOT_MATERIALISED: (
        "This run's snapshot has not been preserved or materialised. "
        "Structures are not fetched from the current PDB, because a newer "
        "entry may differ from the one this result was computed on."
    ),
    CHAIN_NOT_FOUND: (
        "The named chain is not present in the resolved structure file."
    ),
    BACKBONE_MISMATCH: (
        "The two chains do not expose a comparable backbone atom set, so this "
        "view cannot be constructed."
    ),
    UNSUPPORTED_SOURCE: (
        "The source format for this snapshot is not supported by the viewer."
    ),
    SOURCE_NOT_MATERIALISED: (
        "Source object identified but not yet held locally. It can be "
        "materialised on demand from the run's own snapshot."
    ),
    SOURCE_MANIFEST_MISSING: (
        "This run has no Bronze source manifest, so the snapshot's source "
        "object identity cannot be established."
    ),
    NOT_IN_SNAPSHOT: (
        "This entry is not listed in the snapshot's Bronze source manifest."
    ),
    SOURCE_FETCH_FAILED: (
        "The source object could not be retrieved or failed ETag "
        "verification against the run's manifest."
    ),
}


# ---------------------------------------------------------------------------
# Views
# ---------------------------------------------------------------------------

VIEW_SIDE_BY_SIDE = "side_by_side"
VIEW_SUPERPOSED = "superposed"
VIEW_DEPOSITED = "deposited"
VIEW_CHAINS_ONLY = "chains_only"

VIEWS: tuple[tuple[str, str, str], ...] = (
    (
        VIEW_SIDE_BY_SIDE,
        "Side by side",
        "Both chains in their deposited orientation, translated apart. No "
        "rotation is applied.",
    ),
    (
        VIEW_SUPERPOSED,
        "Superposed",
        "The second chain is rigidly superposed onto the first (Kabsch, on "
        "paired backbone atoms) for visual comparison only.",
    ),
    (
        VIEW_CHAINS_ONLY,
        "Chains only",
        "Only the two chains of the pair, without the rest of each "
        "deposited entry.",
    ),
    (
        VIEW_DEPOSITED,
        "Deposited context",
        "Each chain shown within its full deposited structure.",
    ),
)


class MolstarServiceError(RuntimeError):
    """Raised when a scene cannot be produced for a reason worth surfacing."""


# ---------------------------------------------------------------------------
# Structure resolution
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StructureSource:
    """Where one structure was resolved from, and under which snapshot."""

    pdb_id: str
    path: Path
    origin: str
    snapshot_id: str | None

    #: The Bronze source-object identity this structure corresponds to.
    source_object: Any = None

    def to_dict(self) -> dict[str, Any]:
        payload = {
            "pdb_id": self.pdb_id,
            "path": str(self.path),
            "origin": self.origin,
            "snapshot_id": self.snapshot_id,
        }

        if self.source_object is not None:
            payload["source_object"] = self.source_object.to_dict()

        return payload


@dataclass(frozen=True)
class StructureLocator:
    """Resolves a PDB id to a structure file for one run's snapshot.

    Resolution is anchored to the run's **source** layer, never to whether the
    chain survived Stage 14. A chain removed from the Gold retained population
    is still a deposited structure in the snapshot, and stays inspectable.

    Search order, strongest provenance first:

    1. the hot materialisation of *this run's* snapshot;
    2. the durable preserved store for *this run's* snapshot;
    3. prepared example structures shipped with the repository;
    4. the run's **Bronze source manifest**, which gives the exact
       snapshot-scoped object key and ETag -- materialised on demand into (1)
       when ``allow_materialisation`` is set.

    Step 4 fetches the *snapshot* object, e.g.
    ``20260101/pub/.../7acj.cif.gz``, and verifies its ETag against the run's
    manifest. It is therefore the object the pipeline actually parsed, not
    whatever the PDB serves today.
    """

    snapshot_id: str | None = None
    hot_root: Path | None = None
    durable_root: Path | None = None
    example_root: Path | None = None

    #: The run's Bronze source manifest index, if one is available.
    source_index: Any = None

    #: The run's chain-namespace index (label_asym_id <-> auth_asym_id).
    chain_index: Any = None

    #: Bucket the snapshot objects are served from.
    bucket_url: str | None = None

    #: Whether a single missing object may be fetched on demand.
    allow_materialisation: bool = True

    def candidate_paths(self, pdb_id: str) -> list[tuple[Path, str]]:
        lowered = pdb_id.lower()
        candidates: list[tuple[Path, str]] = []

        if self.hot_root is not None and self.snapshot_id:
            base = self.hot_root / self.snapshot_id
            for suffix in (".cif", ".cif.gz", ".bcif"):
                candidates.append(
                    (base / f"{lowered}{suffix}", "hot_cache")
                )
                candidates.append(
                    (base / lowered[1:3] / f"{lowered}{suffix}", "hot_cache")
                )

        if self.durable_root is not None and self.snapshot_id:
            base = self.durable_root / "materialised" / self.snapshot_id
            for suffix in (".cif", ".cif.gz"):
                candidates.append(
                    (base / f"{lowered}{suffix}", "durable_store")
                )

        if self.example_root is not None:
            candidates.append(
                (self.example_root / f"{lowered}.cif", "prepared_example")
            )

        return candidates

    def source_object(self, pdb_id: str):
        """The Bronze source-object identity for this entry, if known."""

        if self.source_index is None:
            return None

        return self.source_index.lookup(pdb_id)

    def resolve(
        self,
        pdb_id: str,
        *,
        materialise_if_missing: bool | None = None,
    ) -> StructureSource | None:
        """Resolve one structure, materialising from source if permitted."""

        for path, origin in self.candidate_paths(pdb_id):
            if path.is_file():
                return StructureSource(
                    pdb_id=pdb_id.lower(),
                    path=path,
                    origin=origin,
                    snapshot_id=(
                        self.snapshot_id
                        if origin != "prepared_example"
                        else None
                    ),
                    source_object=self.source_object(pdb_id),
                )

        allowed = (
            self.allow_materialisation
            if materialise_if_missing is None
            else materialise_if_missing
        )

        if not allowed:
            return None

        return self._materialise(pdb_id)

    def _materialise(self, pdb_id: str) -> StructureSource | None:
        """Fetch exactly this snapshot's object for one entry."""

        from pdbclean.source_index import (
            SourceIndexError,
            cache_path_for,
            materialise,
        )

        found = self.source_object(pdb_id)

        if found is None or not self.bucket_url or self.hot_root is None:
            return None

        destination = cache_path_for(found, hot_root=self.hot_root)

        try:
            path = materialise(
                found,
                bucket_url=self.bucket_url,
                destination=destination,
            )
        except SourceIndexError:
            return None

        return StructureSource(
            pdb_id=found.pdb_id,
            path=path,
            origin="materialised_from_snapshot",
            snapshot_id=found.snapshot_id,
            source_object=found,
        )


def locator_for_run(
    *,
    repo_root: str | Path,
    snapshot_id: str | None,
    hot_cache_root: str | Path | None = None,
    durable_root: str | Path | None = None,
    output_root: str | Path | None = None,
    bucket_url: str | None = None,
    protocol: str | None = None,
    allow_materialisation: bool = True,
) -> StructureLocator:
    """Build a locator scoped to one run's snapshot.

    The Bronze source manifest for that snapshot is attached when present, so
    a structure that is not held locally can still be identified -- and
    materialised on demand -- from the run's own immutable source record.
    """

    root = Path(repo_root)

    index = None

    if output_root is not None and snapshot_id:
        from pdbclean.source_index import BronzeSourceIndex

        candidate = Path(output_root)

        if not candidate.is_absolute():
            candidate = root / candidate

        index = BronzeSourceIndex.for_snapshot(
            output_root=candidate, snapshot_id=snapshot_id
        )

    chains = None

    if output_root is not None and snapshot_id and protocol:
        from pdbclean.source_index import ChainNameIndex

        candidate = Path(output_root)

        if not candidate.is_absolute():
            candidate = root / candidate

        chains = ChainNameIndex.for_protocol(
            output_root=candidate,
            snapshot_id=snapshot_id,
            protocol=protocol,
        )

    return StructureLocator(
        snapshot_id=snapshot_id,
        hot_root=Path(hot_cache_root) if hot_cache_root else None,
        durable_root=Path(durable_root) if durable_root else None,
        example_root=root / "reports" / "molstar_exact_duplicate_examples",
        source_index=index,
        chain_index=chains,
        bucket_url=bucket_url,
        allow_materialisation=allow_materialisation,
    )


# ---------------------------------------------------------------------------
# Pair identity and cache keys
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PairRequest:
    """One duplicate pair, as recorded by the pipeline."""

    pdb_id_a: str
    chain_a: str
    pdb_id_b: str
    chain_b: str
    snapshot_id: str | None = None
    run_id: str | None = None
    model_id: int = 1
    chain_length: int | None = None
    d_bri_units: int | None = None
    classification: str | None = None
    relationship: str | None = None
    representative: str | None = None

    @property
    def label(self) -> str:
        return (
            f"{self.pdb_id_a}:{self.chain_a} <-> "
            f"{self.pdb_id_b}:{self.chain_b}"
        )

    def cache_key(self, view: str) -> str:
        """Deterministic key over run, snapshot, pair and view.

        Two requests that are scientifically the same pair, in the same run,
        for the same view, reuse the same scene.
        """

        material = "|".join(
            str(part)
            for part in (
                SCENE_SCHEMA_VERSION,
                self.run_id or "-",
                self.snapshot_id or "-",
                self.pdb_id_a.lower(),
                self.chain_a,
                self.pdb_id_b.lower(),
                self.chain_b,
                self.model_id,
                view,
            )
        )

        return hashlib.sha256(material.encode("utf-8")).hexdigest()[:32]

    def to_dict(self) -> dict[str, Any]:
        return {
            "pdb_id_a": self.pdb_id_a,
            "chain_a": self.chain_a,
            "pdb_id_b": self.pdb_id_b,
            "chain_b": self.chain_b,
            "snapshot_id": self.snapshot_id,
            "run_id": self.run_id,
            "model_id": self.model_id,
            "chain_length": self.chain_length,
            "d_bri_units": self.d_bri_units,
            "classification": self.classification,
            "relationship": self.relationship,
            "representative": self.representative,
        }


# ---------------------------------------------------------------------------
# Scene construction
# ---------------------------------------------------------------------------



#: A cartoon needs a secondary-structure run to draw anything. Chains shorter
#: than this render as ball-and-stick so their atoms are actually visible --
#: many detected duplicate chains are only two or three residues long.
MINIMUM_CARTOON_RESIDUES = 12


def _representation_for(chain_length: int | None) -> str:
    """Pick a representation that is actually visible for this chain."""

    if chain_length is not None and chain_length < MINIMUM_CARTOON_RESIDUES:
        return "ball_and_stick"

    return "cartoon"


def _chain_length(path: Path, chain_id: str) -> int | None:
    """Residue count of one label_asym_id chain, from the structure itself."""

    try:
        keys, _ = scenes.backbone_coordinates(path, chain_id)
    except Exception:  # pragma: no cover - unreadable chain
        return None

    return len(keys) // 3 if keys else None


def _chain_selector(chain_id: str) -> list[dict[str, Any]]:
    """MolViewSpec selector for one chain, in the frozen scenes' format."""

    return [{"label_asym_id": chain_id}]


def _all_chains_selector(path: Path) -> list[dict[str, Any]]:
    """Selector covering every polymer chain in a structure.

    Reads ``_atom_site.label_asym_id`` directly. gemmi's ``chain.name`` is the
    **auth** identifier, which is a different namespace and would silently
    select nothing when used against a ``label_asym_id`` selector.
    """

    document = gemmi.cif.read(str(path))
    block = document.sole_block()

    labels = block.find_loop("_atom_site.label_asym_id")

    chains = sorted({value for value in labels if value})

    return [{"label_asym_id": chain} for chain in chains] or [{}]


def _build_side_by_side(
    source_a: StructureSource,
    source_b: StructureSource,
    request: PairRequest,
) -> dict[str, Any]:
    xyz_a = scenes.selected_atom_coordinates(source_a.path, [request.chain_a])
    xyz_b = scenes.selected_atom_coordinates(source_b.path, [request.chain_b])

    left, right, _separation = scenes.side_by_side_translations(
        xyz_a, xyz_b
    )

    return scenes.scene(
        f"{request.label} — side by side",
        [
            scenes.structure_branch(
                f"./{source_a.pdb_id}.cif",
                _chain_selector(request.chain_a),
                COLOUR_A,
                rep_type=_representation_for(
                    _chain_length(source_a.path, request.chain_a)
                ),
                translation=left,
            ),
            scenes.structure_branch(
                f"./{source_b.pdb_id}.cif",
                _chain_selector(request.chain_b),
                COLOUR_B,
                rep_type=_representation_for(
                    _chain_length(source_b.path, request.chain_b)
                ),
                translation=right,
            ),
        ],
    )


def _build_superposed(
    source_a: StructureSource,
    source_b: StructureSource,
    request: PairRequest,
) -> tuple[dict[str, Any], dict[str, Any]]:
    keys_a, _ = scenes.backbone_coordinates(source_a.path, request.chain_a)
    keys_b, _ = scenes.backbone_coordinates(source_b.path, request.chain_b)

    common = len(set(keys_a) & set(keys_b))

    if common < 3:
        raise MolstarServiceError(BACKBONE_MISMATCH)

    _, ref, mov = scenes.paired_backbones(
        source_a.path,
        request.chain_a,
        source_b.path,
        request.chain_b,
        expected_residues=common // 3,
    )

    rotation, translation = scenes.kabsch_reference_from_moving(ref, mov)
    aligned = scenes.apply_transform(mov, rotation, translation)

    # Reported for the viewer only. These are visual-agreement figures; the
    # authoritative duplicate distance is the recorded complete-BRI L-infinity.
    metrics = {
        "paired_backbone_atom_count": int(ref.shape[0]),
        "aligned_backbone_rmsd_angstrom": scenes.rmsd(ref, aligned),
        "aligned_backbone_max_distance_angstrom": (
            scenes.max_atom_distance(ref, aligned)
        ),
        "note": (
            "Visual superposition figures only. The authoritative duplicate "
            "distance is the recorded complete-BRI L-infinity value."
        ),
    }

    document = scenes.scene(
        f"{request.label} — superposed",
        [
            scenes.structure_branch(
                f"./{source_a.pdb_id}.cif",
                _chain_selector(request.chain_a),
                COLOUR_A,
                rep_type=_representation_for(
                    _chain_length(source_a.path, request.chain_a)
                ),
            ),
            scenes.structure_branch(
                f"./{source_b.pdb_id}.cif",
                _chain_selector(request.chain_b),
                COLOUR_B,
                rep_type=_representation_for(
                    _chain_length(source_b.path, request.chain_b)
                ),
                matrix=scenes.mvs_matrix(rotation, translation),
            ),
        ],
    )

    return document, metrics


def _build_chains_only(
    source_a: StructureSource,
    source_b: StructureSource,
    request: PairRequest,
) -> dict[str, Any]:
    return scenes.scene(
        f"{request.label} — chains only",
        [
            scenes.structure_branch(
                f"./{source_a.pdb_id}.cif",
                _chain_selector(request.chain_a),
                COLOUR_A,
                rep_type=_representation_for(
                    _chain_length(source_a.path, request.chain_a)
                ),
            ),
            scenes.structure_branch(
                f"./{source_b.pdb_id}.cif",
                _chain_selector(request.chain_b),
                COLOUR_B,
                rep_type=_representation_for(
                    _chain_length(source_b.path, request.chain_b)
                ),
            ),
        ],
    )


def _build_deposited(
    source_a: StructureSource,
    source_b: StructureSource,
    request: PairRequest,
) -> dict[str, Any]:
    return scenes.scene(
        f"{request.label} — deposited context",
        [
            scenes.structure_branch(
                f"./{source_a.pdb_id}.cif",
                _all_chains_selector(source_a.path),
                COLOUR_A,
            ),
            scenes.structure_branch(
                f"./{source_b.pdb_id}.cif",
                _all_chains_selector(source_b.path),
                COLOUR_B,
            ),
        ],
    )


# ---------------------------------------------------------------------------
# Service
# ---------------------------------------------------------------------------


class MolstarSceneService:
    """Resolves structures, builds scenes on demand and caches them."""

    def __init__(
        self,
        *,
        locator: StructureLocator,
        cache_root: str | Path,
    ) -> None:
        self.locator = locator
        self.cache_root = Path(cache_root)

    # -- availability ---------------------------------------------------

    def availability(self, request: PairRequest) -> dict[str, Any]:
        """Report whether this pair can be visualised, and why not if not.

        Availability depends only on **source** availability. A chain removed
        by Stage 14 is still a deposited structure and remains inspectable;
        its retained/removed status is reported as metadata, never as a
        reason to refuse.

        Never returns a bare "no prepared scene": either the pair is
        visualisable, or the specific blocking reason is named.
        """

        entries: dict[str, dict[str, Any]] = {}

        for pdb_id in (request.pdb_id_a, request.pdb_id_b):
            # Do not fetch during a mere availability check.
            local = self.locator.resolve(pdb_id, materialise_if_missing=False)

            if local is not None:
                entries[pdb_id] = {
                    "state": "available",
                    "source": local.to_dict(),
                }
                continue

            index = self.locator.source_index

            if index is None or not getattr(index, "available", False):
                entries[pdb_id] = {
                    "state": SOURCE_MANIFEST_MISSING,
                    "source": None,
                }
                continue

            found = self.locator.source_object(pdb_id)

            if found is None:
                entries[pdb_id] = {"state": NOT_IN_SNAPSHOT, "source": None}
                continue

            # Identity known, bytes not held yet -- recoverable on demand.
            entries[pdb_id] = {
                "state": SOURCE_NOT_MATERIALISED,
                "source": found.to_dict(),
            }

        states = {entry["state"] for entry in entries.values()}

        recoverable = states <= {"available", SOURCE_NOT_MATERIALISED}

        if recoverable:
            return {
                "available": True,
                "reason": None,
                "reason_text": None,
                "requires_materialisation": (
                    SOURCE_NOT_MATERIALISED in states
                ),
                "missing_structures": [],
                "snapshot_id": self.locator.snapshot_id,
                "structures": entries,
                "sources": [
                    entry["source"]
                    for entry in entries.values()
                    if entry["source"]
                ],
                "views": [
                    {"key": key, "label": label, "description": description}
                    for key, label, description in VIEWS
                ],
            }

        blocking = next(
            state for state in states if state != "available"
        )

        return {
            "available": False,
            "reason": blocking,
            "reason_text": REASON_TEXT.get(blocking, blocking),
            "requires_materialisation": False,
            "missing_structures": [
                pdb_id
                for pdb_id, entry in entries.items()
                if entry["state"] != "available"
            ],
            "snapshot_id": self.locator.snapshot_id,
            "structures": entries,
            "views": [],
        }

    # -- scene generation -----------------------------------------------

    def _cache_path(self, request: PairRequest, view: str) -> Path:
        key = request.cache_key(view)

        return self.cache_root / key[:2] / f"{key}.mvsj"

    def scene(
        self,
        request: PairRequest,
        view: str = VIEW_SIDE_BY_SIDE,
        *,
        use_cache: bool = True,
    ) -> dict[str, Any]:
        """Return a Mol* scene for one pair, generating it if necessary."""

        if view not in {key for key, _, _ in VIEWS}:
            raise MolstarServiceError(f"Unknown view: {view}")

        cache_path = self._cache_path(request, view)

        if use_cache and cache_path.is_file():
            payload = json.loads(cache_path.read_text(encoding="utf-8"))
            payload["cache"] = {"hit": True, "path": str(cache_path)}

            # The cached artefact is the geometry. The recorded classification,
            # distance and Stage-14 relationship belong to this request, so
            # they are refreshed rather than served stale.
            payload["pair"] = request.to_dict()
            payload.setdefault("provenance", {}).update(
                {
                    "run_id": request.run_id,
                    "classification": request.classification,
                    "d_bri_units": request.d_bri_units,
                    "relationship": request.relationship,
                    "representative": request.representative,
                }
            )

            return payload

        source_a = self.locator.resolve(request.pdb_id_a)
        source_b = self.locator.resolve(request.pdb_id_b)

        if source_a is None or source_b is None:
            availability = self.availability(request)

            raise MolstarServiceError(
                availability["reason_text"] or STRUCTURE_UNAVAILABLE
            )

        metrics: dict[str, Any] = {}

        if view == VIEW_SIDE_BY_SIDE:
            document = _build_side_by_side(source_a, source_b, request)
        elif view == VIEW_SUPERPOSED:
            document, metrics = _build_superposed(source_a, source_b, request)
        elif view == VIEW_CHAINS_ONLY:
            document = _build_chains_only(source_a, source_b, request)
        else:
            document = _build_deposited(source_a, source_b, request)

        payload = {
            "schema_version": SCENE_SCHEMA_VERSION,
            "generated_at": datetime.now(timezone.utc).strftime(
                "%Y-%m-%dT%H:%M:%SZ"
            ),
            "view": view,
            "pair": request.to_dict(),
            # Provenance the viewer displays, so it is always obvious which
            # run and snapshot the structures on screen belong to.
            "provenance": {
                "run_id": request.run_id,
                "snapshot_id": self.locator.snapshot_id,
                "structure_a": source_a.to_dict(),
                "structure_b": source_b.to_dict(),
                # Both chain namespaces, so it is explicit which identifier
                # selects the deposited chain. They frequently differ.
                "chain_a": self._chain_names(
                    request.pdb_id_a, request.chain_a
                ),
                "chain_b": self._chain_names(
                    request.pdb_id_b, request.chain_b
                ),
                "classification": request.classification,
                "d_bri_units": request.d_bri_units,
                "relationship": request.relationship,
                "representative": request.representative,
                "authority": (
                    "Classification shown is the recorded complete-BRI "
                    "L-infinity result. Mol* is visual inspection only and "
                    "never determines whether two chains are duplicates."
                ),
            },
            "visual_metrics": metrics,
            "scene": document,
            "cache": {"hit": False, "path": str(cache_path)},
        }

        # The cache is derived data. Failing to write it must never fail the
        # request, and it is never part of any scientific artefact.
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        except OSError:
            payload["cache"] = {"hit": False, "path": None}

        return payload

    def _chain_names(self, pdb_id: str, label_asym_id: str) -> dict[str, Any]:
        """Both chain identifiers, from the run's own cleaning output."""

        index = self.locator.chain_index

        if index is not None and getattr(index, "available", False):
            found = index.lookup(pdb_id, label_asym_id)

            if found is not None:
                return found.to_dict()

        # Never guess that auth == label; report only what is known.
        return {
            "pdb_id": pdb_id.lower(),
            "label_asym_id": label_asym_id,
            "auth_asym_id": None,
            "diverges": None,
            "selected_on": "label_asym_id",
        }

    def structure_paths(self, request: PairRequest) -> list[Path]:
        """Local structure files a viewer needs for this pair."""

        return [
            source.path
            for source in (
                self.locator.resolve(request.pdb_id_a),
                self.locator.resolve(request.pdb_id_b),
            )
            if source is not None
        ]
