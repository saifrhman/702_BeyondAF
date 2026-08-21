"""Durable snapshot preservation and disposable hot materialisation.

A completed run must stay reproducible after a temporary cache expires, and
preservation state must never be able to change what a run scientifically was.
"""

from __future__ import annotations

import pytest

from pdbclean.runconfig import resolve_run_config
from pdbclean.snapshot import ResolvedSnapshot
from pdbclean.snapshot_selection import resolve_snapshot_for_run
from pdbclean.snapshot_store import (
    HOT,
    MATERIALISED,
    PRESERVED,
    REMOTE_AVAILABLE,
    UNKNOWN,
    ObjectIdentity,
    SnapshotManifest,
    SnapshotStoreError,
    SnapshotStoreLayout,
    plan_preservation,
    provenance_block,
    read_snapshot_manifest,
    snapshot_status,
    write_snapshot_manifest,
)


@pytest.fixture()
def layout(tmp_path):
    return SnapshotStoreLayout(
        durable_root=tmp_path / "store", hot_root=tmp_path / "cache"
    )


def _objects(count, *, prefix="20260101", start=0):
    return [
        ObjectIdentity(
            pdb_id=f"1ab{index}",
            source_key=f"{prefix}/pub/1ab{index}.cif.gz",
            size_bytes=1000 + index,
            etag=f"etag{index}",
        )
        for index in range(start, start + count)
    ]


# --------------------------------------------------------------------------
# F. Snapshot pinning happens once
# --------------------------------------------------------------------------


def _resolver(snapshot_config, **_kwargs):
    snapshot_id = snapshot_config.get("snapshot_id") or "20260415"

    return ResolvedSnapshot(
        snapshot_id=snapshot_id,
        layout="pub/pdb/data/structures/divided/mmCIF",
        source_prefix=f"{snapshot_id}/pub",
        sample_mmcif_key=f"{snapshot_id}/pub/1abc.cif.gz",
        selection_mode=snapshot_config.get("mode", "fixed"),
    )


def test_latest_complete_resolves_to_a_concrete_snapshot():
    resolved, snapshot = resolve_snapshot_for_run(
        resolve_run_config(), resolver=_resolver
    )

    assert resolved.get("snapshot.mode") == "latest_complete"
    assert resolved.get("snapshot.snapshot_id") == "20260415"
    assert snapshot.snapshot_id == "20260415"


def test_provenance_retains_both_the_mode_and_the_resolved_identity():
    resolved, _ = resolve_snapshot_for_run(
        resolve_run_config(), resolver=_resolver
    )

    # The selection *mode* is preserved for the audit record ...
    assert resolved.get("snapshot.resolved_selection_mode") == "latest_complete"
    # ... and the concrete identity is what the run is pinned to.
    assert resolved.get("snapshot.snapshot_id") == "20260415"


def test_a_resumed_run_never_re_resolves_latest(tmp_path):
    """Re-loading a pinned run must not drift to a newer snapshot."""

    from pdbclean.run_provenance import RunProvenance

    pinned, _ = resolve_snapshot_for_run(
        resolve_run_config(), resolver=_resolver
    )

    run = RunProvenance.create(
        resolved=pinned,
        run_root=tmp_path / "runs",
        repo_root=tmp_path,
        snapshot={"snapshot_id": "20260415", "selection_mode": "latest_complete"},
    )

    # The archive has since moved on; the run must not follow it.
    def _later(snapshot_config, **_kwargs):
        return ResolvedSnapshot(
            snapshot_id=snapshot_config.get("snapshot_id") or "20260901",
            layout="pub",
            source_prefix="x",
            sample_mmcif_key="y",
            selection_mode="latest_complete",
        )

    reloaded = RunProvenance.load(run.run_dir).resolved_config

    assert reloaded.get("snapshot.snapshot_id") == "20260415"

    again, snapshot = resolve_snapshot_for_run(reloaded, resolver=_later)

    assert again.get("snapshot.snapshot_id") == "20260415"
    assert snapshot.snapshot_id == "20260415"


def test_pinned_snapshot_identity_is_part_of_the_scientific_hash():
    first, _ = resolve_snapshot_for_run(
        resolve_run_config(overrides=["snapshot.snapshot_id=20260101"]),
        resolver=_resolver,
    )
    second, _ = resolve_snapshot_for_run(
        resolve_run_config(overrides=["snapshot.snapshot_id=20260415"]),
        resolver=_resolver,
    )

    assert first.scientific_sha256 != second.scientific_sha256


# --------------------------------------------------------------------------
# G. Preservation state never alters scientific identity
# --------------------------------------------------------------------------


def test_availability_states_are_reported_accurately(layout):
    assert snapshot_status("20260101", layout)["availability"] == UNKNOWN
    assert snapshot_status(
        "20260101", layout, remote_available=True
    )["availability"] == REMOTE_AVAILABLE

    layout.hot_path("20260101").mkdir(parents=True)

    assert snapshot_status("20260101", layout)["availability"] == HOT

    write_snapshot_manifest(
        SnapshotManifest("20260101", "https://x", "20260101/pub", _objects(2)),
        layout,
    )

    assert snapshot_status("20260101", layout)["availability"] == MATERIALISED


def test_preserved_without_a_cache_is_still_reproducible(layout):
    write_snapshot_manifest(
        SnapshotManifest("20260101", "https://x", "20260101/pub", _objects(3)),
        layout,
    )

    status = snapshot_status("20260101", layout)

    assert status["availability"] == PRESERVED
    assert status["reproducible_without_cache"] is True
    assert status["preserved_object_count"] == 3


def test_a_deleted_hot_cache_does_not_lose_reproducibility(layout):
    write_snapshot_manifest(
        SnapshotManifest("20260101", "https://x", "20260101/pub", _objects(3)),
        layout,
    )

    hot = layout.hot_path("20260101")
    hot.mkdir(parents=True)

    assert snapshot_status("20260101", layout)["availability"] == MATERIALISED

    hot.rmdir()   # the cache expires

    status = snapshot_status("20260101", layout)

    assert status["availability"] == PRESERVED
    assert status["reproducible_without_cache"] is True


def test_preservation_status_is_not_in_the_scientific_configuration():
    resolved = resolve_run_config(overrides=["snapshot.snapshot_id=20260101"])

    from pdbclean.runconfig import scientific_projection

    projection = scientific_projection(resolved.data)

    assert set(projection["snapshot"]) == {"snapshot_id", "bucket_url"}
    assert "availability" not in projection["snapshot"]
    assert "preservation_status" not in projection["snapshot"]


def test_storage_roots_do_not_change_the_scientific_identity(tmp_path):
    baseline = resolve_run_config(overrides=["snapshot.snapshot_id=20260101"])
    relocated = resolve_run_config(
        overrides={
            "snapshot.snapshot_id": "20260101",
            "storage.durable_snapshot_root": str(tmp_path / "elsewhere"),
            "storage.hot_cache_root": str(tmp_path / "scratch"),
        }
    )

    assert relocated.scientific_sha256 == baseline.scientific_sha256


def test_provenance_block_records_availability_without_claiming_identity(layout):
    write_snapshot_manifest(
        SnapshotManifest("20260101", "https://x", "20260101/pub", _objects(1)),
        layout,
    )

    block = provenance_block(
        "20260101", layout, selection_mode="latest_complete"
    )

    assert block["snapshot_id"] == "20260101"
    assert block["selection_mode"] == "latest_complete"
    assert block["availability"] == PRESERVED
    assert "scientific identity" in block["note"]


# --------------------------------------------------------------------------
# Storage roots are configurable, never hard-coded
# --------------------------------------------------------------------------


def test_storage_roots_come_from_configuration(tmp_path):
    resolved = resolve_run_config(
        overrides={
            "storage.durable_snapshot_root": str(tmp_path / "durable"),
            "storage.hot_cache_root": str(tmp_path / "hot"),
        }
    )

    layout = SnapshotStoreLayout.from_config(resolved, repo_root=tmp_path)

    assert layout.durable_root == tmp_path / "durable"
    assert layout.hot_root == tmp_path / "hot"


def test_relative_roots_resolve_against_the_repository(tmp_path):
    layout = SnapshotStoreLayout.from_config(
        resolve_run_config(), repo_root=tmp_path
    )

    assert layout.durable_root == tmp_path / "outputs" / "snapshot_store"
    assert layout.hot_root == tmp_path / "outputs" / "snapshot_cache"


# --------------------------------------------------------------------------
# Content-addressed identity and incremental preservation
# --------------------------------------------------------------------------


def test_content_hash_is_preferred_over_the_etag():
    obj = ObjectIdentity(
        "1abc", "k", 10, etag="abc", content_sha256="cd" * 32
    )

    assert obj.content_id.startswith("sha256:")
    assert obj.verified is True


def test_etag_identity_is_accepted_when_no_hash_exists():
    obj = ObjectIdentity("1abc", "k", 10, etag='"abc123"')

    assert obj.content_id == "etag:abc123:10"
    assert obj.verified is False


def test_filename_only_identity_is_refused():
    with pytest.raises(SnapshotStoreError, match="filename-only"):
        ObjectIdentity("1abc", "k", 10).content_id


def test_unchanged_objects_are_preserved_once_across_snapshots(layout):
    shared = _objects(5)
    added = _objects(1, prefix="20260415", start=99)

    first = plan_preservation(shared, layout)

    assert first["to_transfer_count"] == 5

    # Preserve them.
    write_snapshot_manifest(
        SnapshotManifest("20260101", "https://x", "20260101/pub", shared),
        layout,
    )

    for obj in shared:
        path = layout.object_path(obj.content_id)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"preserved")

    second = plan_preservation(shared + added, layout)

    assert second["already_preserved_count"] == 5
    assert second["to_transfer_count"] == 1
    assert second["bytes_saved_by_deduplication"] > 0
    assert second["to_transfer_bytes"] < second[
        "total_bytes_if_copied_naively"
    ]


def test_repeated_objects_within_one_snapshot_count_once(layout):
    obj = _objects(1)[0]

    plan = plan_preservation([obj, obj, obj], layout)

    assert plan["object_count"] == 3
    assert plan["distinct_object_count"] == 1
    assert plan["duplicate_references_within_snapshot"] == 2
    assert plan["to_transfer_count"] == 1


def test_planning_moves_no_data(layout):
    plan_preservation(_objects(10), layout)

    assert not layout.objects_root.exists()
    assert not layout.manifests_root.exists()


def test_object_paths_are_sharded(layout):
    path = layout.object_path("sha256:" + "ab" * 32)

    assert path.parent.parent.parent.name == "sha256"
    assert len(path.parent.parent.name) == 2
    assert len(path.parent.name) == 2


# --------------------------------------------------------------------------
# Manifests are immutable
# --------------------------------------------------------------------------


def test_manifests_round_trip(layout):
    manifest = SnapshotManifest(
        "20260101", "https://x", "20260101/pub", _objects(4)
    )

    write_snapshot_manifest(manifest, layout)

    reloaded = read_snapshot_manifest("20260101", layout)

    assert reloaded is not None
    assert reloaded.snapshot_id == "20260101"
    assert len(reloaded.objects) == 4
    assert reloaded.objects[0].content_id == manifest.objects[0].content_id


def test_an_existing_manifest_is_never_mutated(layout):
    manifest = SnapshotManifest(
        "20260101", "https://x", "20260101/pub", _objects(2)
    )

    write_snapshot_manifest(manifest, layout)

    with pytest.raises(SnapshotStoreError, match="immutable"):
        write_snapshot_manifest(
            SnapshotManifest(
                "20260101", "https://x", "20260101/pub", _objects(9)
            ),
            layout,
        )

    assert len(read_snapshot_manifest("20260101", layout).objects) == 2


def test_missing_manifest_reads_as_none(layout):
    assert read_snapshot_manifest("29991231", layout) is None


def test_manifest_records_the_documented_fields(layout):
    manifest = SnapshotManifest(
        "20260101", "https://x", "20260101/pub", _objects(2)
    )

    payload = manifest.to_dict()

    assert payload["snapshot_id"] == "20260101"
    assert payload["object_count"] == 2
    assert payload["total_bytes"] > 0

    entry = payload["objects"][0]

    for field in (
        "pdb_id",
        "source_key",
        "size_bytes",
        "etag",
        "content_sha256",
        "content_id",
        "identity_verified",
    ):
        assert field in entry, field
