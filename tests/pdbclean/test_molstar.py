"""Mol* on demand: reuse, snapshot correctness and inspection-only status.

Mol* must be usable for the pairs the pipeline actually detected, must show
structures belonging to the run's own snapshot, and must never influence a
classification.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from pdbclean import molstar_scenes as scenes
from pdbclean.molstar_service import (
    NOT_IN_SNAPSHOT,
    SNAPSHOT_NOT_MATERIALISED,
    SOURCE_MANIFEST_MISSING,
    STRUCTURE_UNAVAILABLE,
    VIEW_SIDE_BY_SIDE,
    VIEW_SUPERPOSED,
    VIEWS,
    MolstarSceneService,
    MolstarServiceError,
    PairRequest,
    StructureLocator,
    locator_for_run,
)


REPO = Path(__file__).resolve().parents[2]
EXAMPLES = REPO / "reports" / "molstar_exact_duplicate_examples"


def _examples() -> Path:
    if not (EXAMPLES / "1a0t.cif").is_file():
        pytest.skip("prepared example structures are not present")

    return EXAMPLES


@pytest.fixture()
def service(tmp_path):
    _examples()

    return MolstarSceneService(
        locator=locator_for_run(repo_root=REPO, snapshot_id="20260101"),
        cache_root=tmp_path / "scenes",
    )


@pytest.fixture()
def known_pair():
    return PairRequest(
        pdb_id_a="1a0t",
        chain_a="B",
        pdb_id_b="1oh2",
        chain_b="B",
        snapshot_id="20260101",
        run_id="run-test",
        chain_length=413,
        d_bri_units=0,
        classification="exact_duplicate",
    )


# --------------------------------------------------------------------------
# The extracted core still reproduces the frozen scenes
# --------------------------------------------------------------------------


def test_extracted_core_reproduces_the_frozen_superposition():
    """Proof the scene logic was reused, not reimplemented."""

    directory = _examples()

    _, reference, moving = scenes.paired_backbones(
        directory / "1a0t.cif",
        "B",
        directory / "1oh2.cif",
        "B",
        expected_residues=413,
    )

    rotation, translation = scenes.kabsch_reference_from_moving(
        reference, moving
    )
    matrix = scenes.mvs_matrix(rotation, translation)

    frozen = json.loads(
        (directory / "pair1_chainB_superposed.mvsj").read_text(
            encoding="utf-8"
        )
    )

    def _matrices(node):
        if isinstance(node, dict):
            if node.get("kind") == "transform":
                found = node.get("params", {}).get("matrix")

                if found:
                    yield found

            for value in node.values():
                yield from _matrices(value)
        elif isinstance(node, list):
            for value in node:
                yield from _matrices(value)

    assert any(
        np.allclose(np.array(candidate), np.array(matrix), atol=1e-9)
        for candidate in _matrices(frozen)
    )


def test_extracted_core_reproduces_the_frozen_metrics():
    directory = _examples()

    metrics = json.loads(
        (directory / "metrics.json").read_text(encoding="utf-8")
    )["pair1"]

    _, reference, moving = scenes.paired_backbones(
        directory / "1a0t.cif",
        "B",
        directory / "1oh2.cif",
        "B",
        expected_residues=metrics["retained_residue_count"],
    )

    rotation, translation = scenes.kabsch_reference_from_moving(
        reference, moving
    )
    aligned = scenes.apply_transform(moving, rotation, translation)

    assert reference.shape[0] == metrics["backbone_atom_count"]
    assert scenes.rmsd(reference, aligned) == pytest.approx(
        metrics["aligned_backbone_rmsd_A"], abs=1e-9
    )


def test_scene_module_has_no_import_side_effects(tmp_path, monkeypatch):
    """Importing the core must not write scenes or print a report."""

    monkeypatch.chdir(tmp_path)

    import importlib

    importlib.reload(scenes)

    assert not list(tmp_path.iterdir())


# --------------------------------------------------------------------------
# Availability, and what replaced "no prepared scene"
# --------------------------------------------------------------------------


def test_pair_with_local_structures_is_visualisable(service, known_pair):
    availability = service.availability(known_pair)

    assert availability["available"] is True
    assert availability["reason"] is None
    assert {view["key"] for view in availability["views"]} == {
        key for key, _, _ in VIEWS
    }


def test_unavailable_pair_names_a_specific_reason(service):
    request = PairRequest(
        pdb_id_a="7acj",
        chain_a="Z",
        pdb_id_b="7acr",
        chain_b="Z",
        snapshot_id="20260101",
    )

    availability = service.availability(request)

    assert availability["available"] is False
    # The reason names what is actually missing, rather than a generic
    # "unavailable".
    assert availability["reason"] in {
        STRUCTURE_UNAVAILABLE,
        SNAPSHOT_NOT_MATERIALISED,
        SOURCE_MANIFEST_MISSING,
        NOT_IN_SNAPSHOT,
    }
    assert set(availability["missing_structures"]) == {"7acj", "7acr"}

    # Never the old ambiguous message.
    assert "no prepared scene" not in availability["reason_text"].lower()
    assert len(availability["reason_text"]) > 40


def test_reason_names_the_absent_source_layer(tmp_path):
    """Without a Bronze manifest the viewer says so, rather than guessing."""

    from pdbclean.molstar_service import SOURCE_MANIFEST_MISSING

    for locator in (
        StructureLocator(
            snapshot_id="20260101",
            hot_root=tmp_path / "hot",
            durable_root=tmp_path / "durable",
        ),
        StructureLocator(snapshot_id="20260101"),
    ):
        service = MolstarSceneService(
            locator=locator, cache_root=tmp_path / "scenes"
        )

        request = PairRequest(
            "1abc", "A", "2abc", "B", snapshot_id="20260101"
        )

        assert service.availability(request)["reason"] == (
            SOURCE_MANIFEST_MISSING
        )


# --------------------------------------------------------------------------
# Scene generation and caching
# --------------------------------------------------------------------------


@pytest.mark.parametrize("view", [key for key, _, _ in VIEWS])
def test_every_view_can_be_generated(service, known_pair, view):
    payload = service.scene(known_pair, view)

    assert payload["view"] == view
    assert payload["scene"]["root"]["children"]
    assert payload["pair"]["pdb_id_a"] == "1a0t"


def test_scene_is_cached_and_reused(service, known_pair):
    first = service.scene(known_pair, VIEW_SIDE_BY_SIDE)
    second = service.scene(known_pair, VIEW_SIDE_BY_SIDE)

    assert first["cache"]["hit"] is False
    assert second["cache"]["hit"] is True
    assert Path(first["cache"]["path"]).is_file()


def test_cache_key_covers_run_snapshot_pair_and_view(known_pair):
    other_view = known_pair.cache_key(VIEW_SUPERPOSED)
    other_run = PairRequest(
        **{**known_pair.to_dict(), "run_id": "run-other"}
    ).cache_key(VIEW_SIDE_BY_SIDE)
    other_snapshot = PairRequest(
        **{**known_pair.to_dict(), "snapshot_id": "20260415"}
    ).cache_key(VIEW_SIDE_BY_SIDE)

    base = known_pair.cache_key(VIEW_SIDE_BY_SIDE)

    assert len({base, other_view, other_run, other_snapshot}) == 4


def test_cache_is_regenerable(service, known_pair):
    first = service.scene(known_pair, VIEW_SIDE_BY_SIDE)

    Path(first["cache"]["path"]).unlink()

    regenerated = service.scene(known_pair, VIEW_SIDE_BY_SIDE)

    assert regenerated["cache"]["hit"] is False

    # Structurally identical. Only the generation timestamp in the scene
    # metadata differs, which is what makes the cache safe to discard.
    assert regenerated["scene"]["root"] == first["scene"]["root"]


def test_cache_is_not_a_scientific_artefact(service, known_pair, tmp_path):
    """Generating a scene must not touch any release or configuration hash."""

    from pdbclean.runconfig import resolve_run_config

    before = resolve_run_config(
        config_path=REPO / "config/pdbclean/profiles/comp702_frozen_20260101.yaml"
    )

    service.scene(known_pair, VIEW_SUPERPOSED)

    after = resolve_run_config(
        config_path=REPO / "config/pdbclean/profiles/comp702_frozen_20260101.yaml"
    )

    assert after.scientific_sha256 == before.scientific_sha256
    assert after.sha256 == before.sha256

    # And the cache lives outside every release directory.
    assert "releases" not in str(service.cache_root)


def test_generating_a_scene_does_not_modify_the_release(service, known_pair):
    manifest = (
        REPO
        / "outputs"
        / "releases"
        / "PDBClean-20260101-protocol3.2-comp702-v1-dedup-v1"
        / "release_manifest.json"
    )

    if not manifest.is_file():
        pytest.skip("frozen release not present")

    before = hashlib.sha256(manifest.read_bytes()).hexdigest()

    service.scene(known_pair, VIEW_SUPERPOSED)

    assert hashlib.sha256(manifest.read_bytes()).hexdigest() == before


def test_generating_a_scene_does_not_modify_the_source_structures(
    service, known_pair
):
    directory = _examples()
    target = directory / "1a0t.cif"

    before = hashlib.sha256(target.read_bytes()).hexdigest()

    for view, _, _ in VIEWS:
        service.scene(known_pair, view)

    assert hashlib.sha256(target.read_bytes()).hexdigest() == before


# --------------------------------------------------------------------------
# Snapshot correctness and provenance
# --------------------------------------------------------------------------


def test_scene_records_the_run_and_snapshot(service, known_pair):
    payload = service.scene(known_pair, VIEW_SIDE_BY_SIDE)

    provenance = payload["provenance"]

    assert provenance["run_id"] == "run-test"
    assert provenance["snapshot_id"] == "20260101"
    assert provenance["structure_a"]["pdb_id"] == "1a0t"
    assert provenance["structure_b"]["pdb_id"] == "1oh2"


def test_scene_shows_the_recorded_classification_not_a_new_one(
    service, known_pair
):
    payload = service.scene(known_pair, VIEW_SUPERPOSED)

    provenance = payload["provenance"]

    assert provenance["classification"] == "exact_duplicate"
    assert provenance["d_bri_units"] == 0
    assert "never determines" in provenance["authority"]


def test_visual_metrics_are_labelled_as_non_authoritative(
    service, known_pair
):
    payload = service.scene(known_pair, VIEW_SUPERPOSED)

    metrics = payload["visual_metrics"]

    assert "authoritative" in metrics["note"]
    assert metrics["aligned_backbone_rmsd_angstrom"] < 1e-6


def test_structures_resolve_from_the_runs_snapshot_first(tmp_path):
    """A materialised snapshot copy wins over the prepared example."""

    hot = tmp_path / "hot" / "20260101"
    hot.mkdir(parents=True)
    (hot / "1a0t.cif").write_text("data_fake\n", encoding="utf-8")

    locator = StructureLocator(
        snapshot_id="20260101",
        hot_root=tmp_path / "hot",
        example_root=EXAMPLES,
    )

    source = locator.resolve("1a0t")

    assert source is not None
    assert source.origin == "hot_cache"
    assert source.snapshot_id == "20260101"


def test_no_remote_fetch_is_attempted(service):
    """An absent structure fails locally; it is never downloaded."""

    request = PairRequest("9zzz", "A", "8zzz", "B", snapshot_id="20260101")

    with pytest.raises(MolstarServiceError):
        service.scene(request, VIEW_SIDE_BY_SIDE)


def test_locator_never_lists_a_remote_source():
    locator = locator_for_run(repo_root=REPO, snapshot_id="20260101")

    for path, _origin in locator.candidate_paths("1abc"):
        assert not str(path).startswith("http")


# --------------------------------------------------------------------------
# Prepared example scenes still work
# --------------------------------------------------------------------------


def test_prepared_example_scenes_are_still_present_and_valid():
    directory = _examples()

    prepared = sorted(directory.glob("*.mvsj"))

    assert prepared

    for scene_file in prepared:
        document = json.loads(scene_file.read_text(encoding="utf-8"))

        assert document.get("root")
        assert document["root"].get("children")


def test_unknown_view_is_rejected(service, known_pair):
    with pytest.raises(MolstarServiceError, match="Unknown view"):
        service.scene(known_pair, "hologram")


# --------------------------------------------------------------------------
# Determinism: only non-scientific metadata may vary between generations
# --------------------------------------------------------------------------


def _strip_volatile(payload):
    """Remove the fields that are legitimately generation-time metadata.

    Everything else must be byte-identical between two generations of the
    same pair/view, because it describes what is shown scientifically.
    """

    import copy

    stripped = copy.deepcopy(payload)

    stripped.pop("generated_at", None)
    stripped.pop("cache", None)
    stripped.get("scene", {}).get("metadata", {}).pop("timestamp", None)

    return stripped


def test_only_the_timestamp_varies_between_generations(service, known_pair):
    """Pin exactly which fields are volatile, so none can quietly join them."""

    first = service.scene(known_pair, VIEW_SUPERPOSED, use_cache=False)
    second = service.scene(known_pair, VIEW_SUPERPOSED, use_cache=False)

    # The two differ only in generation metadata ...
    assert first != second
    assert (
        first["scene"]["metadata"]["timestamp"]
        != second["scene"]["metadata"]["timestamp"]
    )

    # ... and are otherwise identical.
    assert _strip_volatile(first) == _strip_volatile(second)


def test_volatile_fields_are_not_scientific(service, known_pair):
    """The volatile metadata carries no scientific content."""

    payload = service.scene(known_pair, VIEW_SUPERPOSED, use_cache=False)

    volatile = {
        "generated_at": payload["generated_at"],
        "timestamp": payload["scene"]["metadata"]["timestamp"],
    }

    for value in volatile.values():
        assert isinstance(value, str)
        # A timestamp, not a distance, classification or identity.
        assert value.startswith("20")

    # None of the scientific facts live in the volatile block.
    stripped = _strip_volatile(payload)

    assert stripped["provenance"]["classification"] == "exact_duplicate"
    assert stripped["provenance"]["d_bri_units"] == 0
    assert stripped["provenance"]["snapshot_id"] == "20260101"


def test_superposition_transform_is_deterministic(service, known_pair):
    """The matrix on screen must not drift between generations."""

    def _matrices(node):
        if isinstance(node, dict):
            if node.get("kind") == "transform":
                found = node.get("params", {}).get("matrix")

                if found:
                    yield tuple(found)

            for value in node.values():
                yield from _matrices(value)
        elif isinstance(node, list):
            for value in node:
                yield from _matrices(value)

    first = list(
        _matrices(service.scene(known_pair, VIEW_SUPERPOSED, use_cache=False))
    )
    second = list(
        _matrices(service.scene(known_pair, VIEW_SUPERPOSED, use_cache=False))
    )

    assert first
    assert first == second


def test_visual_metrics_are_deterministic(service, known_pair):
    first = service.scene(known_pair, VIEW_SUPERPOSED, use_cache=False)
    second = service.scene(known_pair, VIEW_SUPERPOSED, use_cache=False)

    assert first["visual_metrics"] == second["visual_metrics"]


def test_cached_scene_matches_a_freshly_generated_one(service, known_pair):
    """A cache hit must not show anything different from a fresh build."""

    fresh = service.scene(known_pair, VIEW_SIDE_BY_SIDE, use_cache=False)
    cached = service.scene(known_pair, VIEW_SIDE_BY_SIDE)

    assert cached["cache"]["hit"] is True
    assert cached["scene"]["root"] == fresh["scene"]["root"]
    assert cached["provenance"] == fresh["provenance"]


def test_scene_identity_follows_the_snapshot(service, known_pair):
    """The same pair under a different snapshot is a different scene."""

    other = PairRequest(
        **{**known_pair.to_dict(), "snapshot_id": "20260415"}
    )

    assert known_pair.cache_key(VIEW_SIDE_BY_SIDE) != other.cache_key(
        VIEW_SIDE_BY_SIDE
    )


# ==========================================================================
# Source resolution: Gold retention must never gate visualisation
# ==========================================================================

from pdbclean.source_index import (          # noqa: E402
    BronzeSourceIndex,
    ChainNameIndex,
    SourceIndexError,
    SourceObject,
    cache_path_for,
    materialise,
)


OUTPUT_ROOT = REPO / "outputs" / "pdbclean"
SNAPSHOT = "20260101"
PROTOCOL = "protocol3.2-comp702-v1"
BUCKET = "https://pdbsnapshots.s3.us-west-2.amazonaws.com"


def _bronze():
    index = BronzeSourceIndex.for_snapshot(
        output_root=OUTPUT_ROOT, snapshot_id=SNAPSHOT
    )

    if not index.available:
        pytest.skip("frozen Bronze source manifest is not present")

    return index


def _chains():
    index = ChainNameIndex.for_protocol(
        output_root=OUTPUT_ROOT, snapshot_id=SNAPSHOT, protocol=PROTOCOL
    )

    if not index.available:
        pytest.skip("frozen accepted-chain table is not present")

    return index


@pytest.fixture()
def source_locator(tmp_path):
    _bronze()

    return locator_for_run(
        repo_root=REPO,
        snapshot_id=SNAPSHOT,
        hot_cache_root=tmp_path / "hot",
        durable_root=tmp_path / "durable",
        output_root=OUTPUT_ROOT,
        protocol=PROTOCOL,
        bucket_url=BUCKET,
        allow_materialisation=False,
    )


# --------------------------------------------------------------------------
# D. The resolver uses the source layer, not retained_chains
# --------------------------------------------------------------------------


def test_resolver_never_consults_the_gold_retained_set(source_locator):
    """Retention is metadata; it must not appear in resolution at all."""

    for path, _origin in source_locator.candidate_paths("7acr"):
        rendered = str(path)

        assert "retained_chains" not in rendered
        assert "releases" not in rendered


def test_source_identity_is_known_for_a_removed_chain(source_locator):
    """7acr:Z was removed by Stage 14; its source object is still identified."""

    found = source_locator.source_object("7acr")

    assert found is not None
    assert found.snapshot_id == SNAPSHOT
    assert found.s3_key.startswith(f"{SNAPSHOT}/pub/")
    assert found.s3_key.endswith("7acr.cif.gz")
    assert found.etag


def test_removed_chain_is_reported_available_not_missing(tmp_path):
    """The regression this fix addresses."""

    service = MolstarSceneService(
        locator=locator_for_run(
            repo_root=REPO,
            snapshot_id=SNAPSHOT,
            hot_cache_root=tmp_path / "hot",
            durable_root=tmp_path / "durable",
            output_root=OUTPUT_ROOT,
            protocol=PROTOCOL,
            bucket_url=BUCKET,
        ),
        cache_root=tmp_path / "scenes",
    )

    _bronze()

    request = PairRequest(
        pdb_id_a="7acj",
        chain_a="Z",
        pdb_id_b="7acr",
        chain_b="Z",
        snapshot_id=SNAPSHOT,
        d_bri_units=0,
        classification="exact_duplicate",
        relationship="removed",
        representative="7acj:Z",
    )

    availability = service.availability(request)

    assert availability["available"] is True
    assert availability["missing_structures"] == []
    assert availability["requires_materialisation"] is True


def test_retained_and_removed_members_resolve_identically(source_locator):
    """A retained chain and a removed chain resolve the same way."""

    retained = source_locator.source_object("7acj")   # representative
    removed = source_locator.source_object("7acr")    # removed by Stage 14

    assert retained is not None and removed is not None
    assert retained.snapshot_id == removed.snapshot_id == SNAPSHOT
    assert bool(retained.etag) and bool(removed.etag)


# --------------------------------------------------------------------------
# E. Historical snapshot fidelity
# --------------------------------------------------------------------------


def test_source_keys_are_snapshot_scoped(source_locator):
    """The key names the snapshot, so it cannot be the current PDB."""

    for pdb_id in ("7acj", "7acr", "4qby"):
        found = source_locator.source_object(pdb_id)

        assert found is not None
        assert found.s3_key.startswith(f"{SNAPSHOT}/")


def test_url_is_built_from_the_snapshot_key():
    source = SourceObject(
        pdb_id="7acj",
        snapshot_id=SNAPSHOT,
        s3_key=f"{SNAPSHOT}/pub/pdb/data/structures/divided/mmCIF/ac/7acj.cif.gz",
        etag="abc",
    )

    url = source.url(BUCKET)

    assert url.startswith(BUCKET)
    assert f"/{SNAPSHOT}/" in url
    # Never an undated "current entry" endpoint.
    assert "files.rcsb.org" not in url


def test_etag_mismatch_is_refused(tmp_path):
    """A different revision must never be displayed as the run's structure."""

    import io

    source = SourceObject(
        pdb_id="1abc",
        snapshot_id=SNAPSHOT,
        s3_key=f"{SNAPSHOT}/pub/1abc.cif",
        etag="expected-etag",
        size_bytes=10,
    )

    class _Response:
        headers = {"ETag": '"a-different-etag"'}

        def read(self):
            return b"data_1ABC\n"

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    with pytest.raises(SourceIndexError, match="does not match the ETag"):
        materialise(
            source,
            bucket_url=BUCKET,
            destination=tmp_path / "1abc.cif",
            opener=lambda url, timeout=None: _Response(),
        )


def test_matching_etag_is_accepted(tmp_path):
    source = SourceObject(
        pdb_id="1abc",
        snapshot_id=SNAPSHOT,
        s3_key=f"{SNAPSHOT}/pub/1abc.cif",
        etag="good-etag",
        size_bytes=10,
    )

    class _Response:
        headers = {"ETag": '"good-etag"'}

        def read(self):
            return b"data_1ABC\n"

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    written = materialise(
        source,
        bucket_url=BUCKET,
        destination=tmp_path / "1abc.cif",
        opener=lambda url, timeout=None: _Response(),
    )

    assert written.read_text() == "data_1ABC\n"


def test_cache_path_is_scoped_by_snapshot(tmp_path):
    first = cache_path_for(
        SourceObject("7acj", "20260101", "k"), hot_root=tmp_path
    )
    second = cache_path_for(
        SourceObject("7acj", "20260415", "k"), hot_root=tmp_path
    )

    assert first != second
    assert "20260101" in str(first)
    assert "20260415" in str(second)


def test_oversized_object_is_refused(tmp_path):
    from pdbclean.source_index import MAX_OBJECT_BYTES

    source = SourceObject(
        pdb_id="huge",
        snapshot_id=SNAPSHOT,
        s3_key="k",
        size_bytes=MAX_OBJECT_BYTES + 1,
    )

    with pytest.raises(SourceIndexError, match="materialisation limit"):
        materialise(
            source, bucket_url=BUCKET, destination=tmp_path / "huge.cif"
        )


# --------------------------------------------------------------------------
# F. Chain namespace
# --------------------------------------------------------------------------


def test_label_and_auth_chain_identifiers_are_resolved():
    """They differ for most removed chains; never assume they match."""

    index = _chains()

    found = index.lookup("7acr", "Z")

    assert found is not None
    assert found.label_asym_id == "Z"
    assert found.auth_asym_id == "W"
    assert found.diverges is True


def test_divergent_chain_case_is_real_and_selected_correctly():
    index = _chains()

    for pdb_id, label, auth in (("7acj", "Z", "W"), ("1a0e", "B", "D")):
        found = index.lookup(pdb_id, label)

        assert found is not None, (pdb_id, label)
        assert found.auth_asym_id == auth
        assert found.to_dict()["selected_on"] == "label_asym_id"


def test_scene_reports_both_chain_namespaces(tmp_path):
    _chains()

    service = MolstarSceneService(
        locator=locator_for_run(
            repo_root=REPO,
            snapshot_id=SNAPSHOT,
            hot_cache_root=tmp_path / "hot",
            output_root=OUTPUT_ROOT,
            protocol=PROTOCOL,
            bucket_url=BUCKET,
            allow_materialisation=False,
        ),
        cache_root=tmp_path / "scenes",
    )

    names = service._chain_names("7acr", "Z")

    assert names["label_asym_id"] == "Z"
    assert names["auth_asym_id"] == "W"
    assert names["selected_on"] == "label_asym_id"


def test_unknown_chain_reports_null_auth_rather_than_guessing():
    index = _chains()

    assert index.lookup("zzzz", "Q") is None


# --------------------------------------------------------------------------
# I. Meaningful failure reasons
# --------------------------------------------------------------------------


def test_missing_manifest_is_distinguished_from_missing_object(tmp_path):
    from pdbclean.molstar_service import (
        NOT_IN_SNAPSHOT,
        SOURCE_MANIFEST_MISSING,
    )

    no_manifest = MolstarSceneService(
        locator=StructureLocator(snapshot_id=SNAPSHOT),
        cache_root=tmp_path / "a",
    )

    request = PairRequest("1abc", "A", "2abc", "B", snapshot_id=SNAPSHOT)

    assert no_manifest.availability(request)["reason"] == (
        SOURCE_MANIFEST_MISSING
    )

    with_manifest = MolstarSceneService(
        locator=locator_for_run(
            repo_root=REPO,
            snapshot_id=SNAPSHOT,
            hot_cache_root=tmp_path / "hot",
            output_root=OUTPUT_ROOT,
            protocol=PROTOCOL,
            bucket_url=BUCKET,
            allow_materialisation=False,
        ),
        cache_root=tmp_path / "b",
    )

    _bronze()

    assert with_manifest.availability(
        PairRequest("zzzz", "A", "yyyy", "B", snapshot_id=SNAPSHOT)
    )["reason"] == NOT_IN_SNAPSHOT


def test_every_unavailable_reason_has_readable_text():
    from pdbclean.molstar_service import REASON_TEXT

    for reason, text in REASON_TEXT.items():
        assert len(text) > 30, reason
        assert "no prepared scene" not in text.lower()


# --------------------------------------------------------------------------
# G. Scientific immutability
# --------------------------------------------------------------------------


def test_source_resolution_does_not_touch_the_release(source_locator):
    release = (
        REPO
        / "outputs"
        / "releases"
        / "PDBClean-20260101-protocol3.2-comp702-v1-dedup-v1"
    )

    if not release.is_dir():
        pytest.skip("frozen release not present")

    manifest = release / "release_manifest.json"
    retained = release / "data" / "retained_chains.parquet"

    before = (
        hashlib.sha256(manifest.read_bytes()).hexdigest(),
        retained.stat().st_size,
    )

    source_locator.source_object("7acr")
    source_locator.resolve("7acr", materialise_if_missing=False)

    assert (
        hashlib.sha256(manifest.read_bytes()).hexdigest(),
        retained.stat().st_size,
    ) == before


def test_materialisation_writes_only_into_the_cache(tmp_path):
    """Fetched structures land in the disposable cache, nowhere else."""

    source = SourceObject(
        pdb_id="1abc",
        snapshot_id=SNAPSHOT,
        s3_key=f"{SNAPSHOT}/pub/1abc.cif",
        etag="e",
        size_bytes=10,
    )

    class _Response:
        headers = {"ETag": '"e"'}

        def read(self):
            return b"data_1ABC\n"

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    hot = tmp_path / "hot"

    written = materialise(
        source,
        bucket_url=BUCKET,
        destination=cache_path_for(source, hot_root=hot),
        opener=lambda url, timeout=None: _Response(),
    )

    assert hot in written.parents
    assert "releases" not in str(written)
    assert "bronze" not in str(written)


# --------------------------------------------------------------------------
# J. Performance: only the requested pair is prepared
# --------------------------------------------------------------------------


def test_only_requested_entries_are_materialised(tmp_path):
    fetched: list[str] = []

    class _Response:
        def __init__(self, url):
            self.url = url
            self.headers = {}

        def read(self):
            return b"data_X\n"

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    def _opener(url, timeout=None):
        fetched.append(url)

        return _Response(url)

    source = SourceObject("1abc", SNAPSHOT, f"{SNAPSHOT}/pub/1abc.cif")

    materialise(
        source,
        bucket_url=BUCKET,
        destination=tmp_path / "1abc.cif",
        opener=_opener,
    )

    assert len(fetched) == 1


def test_already_cached_object_is_not_refetched(tmp_path):
    calls: list[str] = []

    source = SourceObject("1abc", SNAPSHOT, f"{SNAPSHOT}/pub/1abc.cif")
    destination = tmp_path / "1abc.cif"
    destination.write_text("cached\n", encoding="utf-8")

    def _opener(url, timeout=None):  # pragma: no cover - must not run
        calls.append(url)

        raise AssertionError("should not refetch a cached object")

    written = materialise(
        source, bucket_url=BUCKET, destination=destination, opener=_opener
    )

    assert written.read_text() == "cached\n"
    assert calls == []


def test_bronze_lookup_does_not_load_the_whole_manifest():
    """246k rows must not be materialised to answer one lookup."""

    index = _bronze()

    found = index.lookup("7acj")

    assert found is not None
    # The index caches per-entry, not the whole table.
    assert set(index._cache) == {"7acj"}


# ==========================================================================
# Rendering preconditions
# ==========================================================================
#
# A scene that loads but draws nothing is a failure. These pin the conditions
# that made the viewer render nothing before: a collapsed viewport, a cartoon
# representation on a dipeptide, and an error handler that removed the host.


def test_viewer_layout_gives_the_molstar_host_the_flexible_row():
    """A sibling being shown must not collapse the viewport to zero height."""

    html = (
        REPO / "src" / "pdbclean" / "ui" / "static" / "viewer.html"
    ).read_text(encoding="utf-8")

    body_css = html.split("body {")[1].split("}")[0]
    rows = [
        line for line in body_css.splitlines()
        if "grid-template-rows" in line
    ]

    assert rows, "viewer body must declare its grid rows"

    tracks = rows[0].split(":")[1].strip().rstrip(";").split()

    # One track per top-level child, with exactly one flexible track.
    children = html.count('\n<div id="') + html.count('\n<details id="')

    assert len(tracks) == children, (tracks, children)
    assert tracks.count("1fr") == 1
    assert tracks[-1] == "1fr", "the viewer host must own the flexible row"


def test_molstar_mounts_into_its_own_host():
    """Clearing the mount must not destroy the status overlay."""

    html = (
        REPO / "src" / "pdbclean" / "ui" / "static" / "viewer.html"
    ).read_text(encoding="utf-8")

    assert 'id="molstar-host"' in html
    assert 'Viewer.create("molstar-host"' in html
    assert 'getElementById("molstar-host").innerHTML = ""' in html

    # The host container itself is never wiped.
    assert 'getElementById("app").innerHTML' not in html


def test_load_failure_never_removes_the_viewer():
    """The regression: a failed load used to leave a blank area."""

    html = (
        REPO / "src" / "pdbclean" / "ui" / "static" / "viewer.html"
    ).read_text(encoding="utf-8")

    assert "showError(" in html
    assert "Structure could not be rendered" in html
    assert "Retry" in html

    # No handler replaces the viewer host.
    assert '#app").innerHTML =' not in html


def test_loading_states_are_reported():
    html = (
        REPO / "src" / "pdbclean" / "ui" / "static" / "viewer.html"
    ).read_text(encoding="utf-8")

    for stage in (
        "Resolving snapshot source",
        "Preparing Mol* viewer",
        "Loading structures",
        "Rendering",
    ):
        assert stage in html, stage


def test_short_chains_use_a_visible_representation():
    """Cartoon draws nothing for a dipeptide; many duplicates are dipeptides."""

    from pdbclean.molstar_service import (
        MINIMUM_CARTOON_RESIDUES,
        _representation_for,
    )

    assert _representation_for(2) == "ball_and_stick"
    assert _representation_for(MINIMUM_CARTOON_RESIDUES - 1) == "ball_and_stick"
    assert _representation_for(MINIMUM_CARTOON_RESIDUES) == "cartoon"
    assert _representation_for(413) == "cartoon"
    assert _representation_for(None) == "cartoon"


def test_generated_scene_requests_a_representation(service, known_pair):
    """A structure with no representation node is invisible."""

    for view, _label, _description in VIEWS:
        payload = service.scene(known_pair, view, use_cache=False)

        found: list[str] = []

        def walk(node):
            if isinstance(node, dict):
                if node.get("kind") == "representation":
                    found.append(node["params"]["type"])

                for child in node.get("children", []) or []:
                    walk(child)
            elif isinstance(node, list):
                for child in node:
                    walk(child)

        walk(payload["scene"]["root"])

        assert found, f"{view} produced no representation"
        assert all(
            kind in {"cartoon", "ball_and_stick", "spacefill", "surface"}
            for kind in found
        ), (view, found)


def test_deposited_view_selects_label_asym_ids_that_exist():
    """gemmi chain.name is the auth namespace and would select nothing."""

    from pdbclean.molstar_service import _all_chains_selector

    cached = (
        REPO
        / "outputs"
        / "snapshot_cache"
        / "20260101"
        / "qb"
        / "4qby.cif"
    )

    if not cached.is_file():
        pytest.skip("4qby not materialised in this checkout")

    import gemmi

    block = gemmi.cif.read(str(cached)).sole_block()
    present = {v for v in block.find_loop("_atom_site.label_asym_id") if v}

    selectors = {
        entry["label_asym_id"] for entry in _all_chains_selector(cached)
    }

    assert selectors
    assert selectors <= present
    assert {"DA", "Z"} <= selectors


def test_cached_scene_reports_this_requests_metadata(service, known_pair):
    """Geometry may be cached; the recorded result must not go stale."""

    service.scene(known_pair, VIEW_SIDE_BY_SIDE)

    relabelled = PairRequest(
        **{
            **known_pair.to_dict(),
            "classification": "nonzero_near_duplicate",
            "d_bri_units": 7,
            "relationship": "retained",
        }
    )

    payload = service.scene(relabelled, VIEW_SIDE_BY_SIDE)

    assert payload["cache"]["hit"] is True
    assert payload["provenance"]["classification"] == "nonzero_near_duplicate"
    assert payload["provenance"]["d_bri_units"] == 7
    assert payload["provenance"]["relationship"] == "retained"
