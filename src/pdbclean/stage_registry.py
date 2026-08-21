"""Declarative description of the COMP702 PDBClean pipeline stages.

This registry is metadata only.  It describes *where* each existing stage
writes, *which* resolved configuration values it consumes, and *how* its output
is recognised and validated.  It contains no scientific computation and never
reimplements a stage: the orchestrator in :mod:`pdbclean.pipeline` drives the
existing stage entry points.

Layers follow the Bronze / Silver / Gold lifecycle documented in
``docs/pdbclean/pipeline_spec.md`` section 3:

``bronze``
    Immutable source inventory and provenance.  No scientific filtering.

``silver``
    The deterministic parsed representation.  Intentionally *not persisted*:
    it is reconstructed from the immutable Bronze object identity (S3 key, byte
    size, ETag) by the versioned parser.

``gold``
    Scientifically derived outputs, from cleaning decisions through to the
    final retained-chain release.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


STAGE_REGISTRY_VERSION = "1.1"

#: The registry's ``ordinal`` is an *execution order index*, not a scientific
#: stage number.  The project's canonical scientific vocabulary is Stage 1
#: through Stage 14, fixed by ``docs/PDBCLEAN_2026_FINDINGS_AND_DECISIONS.md``
#: and by the ``pdbclean_stage<N>_*`` schema names embedded in the frozen
#: artefacts.  Every registry entry therefore carries an explicit
#: :attr:`StageSpec.canonical_stage` label, and that label -- never the
#: ordinal -- is what any interface shows.
#:
#: This matters concretely: the near-duplicate graph is the first step of
#: canonical **Stage 14**, and it must never be presented as canonical
#: "Stage 13", which is the manual detailed-investigation subset and is
#: explicitly *not* the global deletion relation.
CANONICAL_PREREQUISITE = "prerequisite"

#: Canonical Stage 1-14 steps that are deliberately not orchestrated by the
#: runner because they were investigation and review passes feeding the paper,
#: not steps on the release path.  Recorded here so the registry never
#: silently omits part of the canonical vocabulary.
NON_ORCHESTRATED_CANONICAL_STAGES: tuple[dict[str, str], ...] = (
    {
        "canonical_stage": "Stage 11",
        "title": "Acta-style downstream investigation",
        "role": "investigation",
        "frozen_output": "acta_downstream_investigation_v2",
        "note": (
            "Review pass over detected duplicates. Not on the release path "
            "and not a deletion relation."
        ),
    },
    {
        "canonical_stage": "Stage 12",
        "title": "Scientific validation gates",
        "role": "validation",
        "frozen_output": "acta_manual_review_manifest_v2",
        "note": (
            "Validation evidence for the review. Not on the release path and "
            "not a deletion relation."
        ),
    },
    {
        "canonical_stage": "Stage 13",
        "title": "Detailed investigation",
        "role": "investigation",
        "frozen_output": "acta_detailed_review_v2",
        "note": (
            "Manual detailed-review subset. Explicitly NOT the global "
            "Stage-14 deletion relation, and never used as one."
        ),
    },
)

LAYER_SNAPSHOT = "snapshot"
LAYER_BRONZE = "bronze"
LAYER_SILVER = "silver"
LAYER_GOLD = "gold"

#: Ordered lifecycle layers, for UI grouping.
LAYERS: tuple[tuple[str, str], ...] = (
    (LAYER_SNAPSHOT, "Snapshot"),
    (LAYER_BRONZE, "Bronze"),
    (LAYER_SILVER, "Silver"),
    (LAYER_GOLD, "Gold"),
)


@dataclass(frozen=True)
class StageSpec:
    """Static description of one pipeline stage."""

    stage_id: str
    ordinal: int
    title: str
    layer: str
    purpose: str

    #: The canonical scientific stage this entry realises, e.g. ``"Stage 7"``,
    #: ``"Stage 14a"``, or :data:`CANONICAL_PREREQUISITE`.  This -- not
    #: :attr:`ordinal` -- is the label every interface displays.
    canonical_stage: str = CANONICAL_PREREQUISITE

    #: Directory holding the stage output, as a template over
    #: ``output_root``, ``snapshot``, ``protocol``, ``release_root`` and
    #: ``release_name``.  ``None`` means the stage persists nothing.
    directory: str | None = None

    #: Success marker relative to :attr:`directory`.
    success_marker: str | None = "_SUCCESS"

    #: Machine-readable stage summary relative to :attr:`directory`.
    summary: str | None = "global_summary.json"

    #: Primary output artefact relative to :attr:`directory`.
    primary_output: str | None = None

    #: Summary keys reported as the stage's input and output counts.
    input_count_keys: tuple[str, ...] = ()
    output_count_keys: tuple[str, ...] = ()

    #: Dotted resolved-config keys whose values this stage consumes.
    scientific_parameters: tuple[str, ...] = ()

    #: Summary keys compared against the resolved configuration to decide
    #: whether an existing output may be reused.  Each entry maps a summary key
    #: to the dotted resolved-config key it must agree with.
    compatibility: dict[str, str] = field(default_factory=dict)

    #: Stage ids that must reach validation PASS before this stage may start.
    depends_on: tuple[str, ...] = ()

    #: What the validation gate checks, in one human sentence.
    validation: str = ""

    #: How the stage is executed.  Informational for the UI and the plan.
    entry_point: str | None = None

    persisted: bool = True


def _protocol_dir(name: str) -> str:
    return "{output_root}/{snapshot}/{protocol}/" + name


STAGES: tuple[StageSpec, ...] = (
    StageSpec(
        stage_id="snapshot",
        canonical_stage="Prerequisite A",
        ordinal=1,
        title="Snapshot resolution",
        layer=LAYER_SNAPSHOT,
        purpose=(
            "Resolve the requested snapshot to a concrete immutable identity "
            "and verify coordinate availability before any processing begins."
        ),
        directory=None,
        success_marker=None,
        summary=None,
        persisted=False,
        scientific_parameters=(
            "snapshot.mode",
            "snapshot.snapshot_id",
            "snapshot.bucket_url",
        ),
        validation=(
            "The snapshot identity is pinned to YYYYMMDD and a sample "
            "coordinate mmCIF object is confirmed to exist."
        ),
        entry_point="pdbclean.snapshot_selection.resolve_snapshot_for_run",
    ),
    StageSpec(
        stage_id="bronze_source_manifest",
        canonical_stage="Prerequisite B",
        ordinal=2,
        title="Bronze source manifest",
        layer=LAYER_BRONZE,
        purpose=(
            "Immutable inventory of every source object in the snapshot: PDB "
            "ID, S3 key, compressed size and ETag. No scientific filtering."
        ),
        directory="{output_root}/{snapshot}/bronze",
        success_marker=None,
        summary="source_manifest_summary.json",
        primary_output="source_manifest.parquet",
        output_count_keys=("row_count", "object_count", "manifest_rows"),
        scientific_parameters=(
            "snapshot.snapshot_id",
            "snapshot.bucket_url",
        ),
        compatibility={"resolved_snapshot": "snapshot.snapshot_id"},
        depends_on=("snapshot",),
        validation=(
            "Manifest rows all carry the resolved snapshot, and the recorded "
            "sizes and ETags match the listed source objects."
        ),
        entry_point="scripts/pdbclean/create_manifest.py",
    ),
    StageSpec(
        stage_id="silver_parse",
        canonical_stage="Prerequisite C",
        ordinal=3,
        title="Silver parsed representation",
        layer=LAYER_SILVER,
        purpose=(
            "Deterministic parsed structure representation. Not persisted by "
            "design: reconstructed from the immutable Bronze object identity "
            "by the versioned parser, so no second copy of the archive is "
            "stored."
        ),
        directory=None,
        success_marker=None,
        summary=None,
        persisted=False,
        scientific_parameters=(
            "selection.models.policy",
            "selection.models.model_id",
            "selection.chain_namespace.canonical",
        ),
        depends_on=("bronze_source_manifest",),
        validation=(
            "Verified transitively: every Gold chain record carries the "
            "source key and ETag it was parsed from."
        ),
        entry_point="pdbclean.mmcif_parser",
    ),
    StageSpec(
        stage_id="structural_cleaning",
        canonical_stage="Stage 1",
        ordinal=4,
        title="Structural cleaning (Protocol 3.2)",
        layer=LAYER_GOLD,
        purpose=(
            "Apply quality rules Q001-Q006 and the terminal-residue exception; "
            "record every chain as accepted, rejected or failed."
        ),
        directory=_protocol_dir("quality"),
        summary="global_summary.json",
        primary_output="merged/accepted.parquet",
        input_count_keys=("candidate_chain_count", "input_chain_count"),
        output_count_keys=("accepted_chain_count",),
        scientific_parameters=(
            "selection.models.model_id",
            "quality_rules.entry_protein.enabled",
            "quality_rules.disorder.enabled",
            "quality_rules.residue_continuity.enabled",
            "quality_rules.backbone_atoms.required_atoms",
            "quality_rules.backbone_distance.minimum_distance_angstrom",
            "quality_rules.amino_acids.enabled",
        ),
        compatibility={
            "snapshot": "snapshot.snapshot_id",
            "cleaning_protocol": "release.protocol_version",
        },
        depends_on=("silver_parse",),
        validation=(
            "Chain accounting closes: accepted + rejected + non-candidate + "
            "errors equals the candidate population."
        ),
        entry_point="scripts/pdbclean/submit_quality_pipeline.sh",
    ),
    StageSpec(
        stage_id="geometric_validation",
        canonical_stage="Stage 2",
        ordinal=5,
        title="Geometric validation",
        layer=LAYER_GOLD,
        purpose=(
            "Quarantine cleaned chains whose backbone geometry is degenerate, "
            "producing the canonical eligible population."
        ),
        directory=_protocol_dir("geometric_validation"),
        primary_output="finalized/eligible.parquet",
        input_count_keys=("input_accepted_chain_count",),
        output_count_keys=("eligible_chain_count",),
        scientific_parameters=(
            "quality_rules.backbone_distance.minimum_distance_angstrom",
            "post_cleaning_geometric_validation.minimum_triangle_angle_degrees",
        ),
        compatibility={
            "snapshot": "snapshot.snapshot_id",
            "cleaning_protocol": "release.protocol_version",
            "configured_minimum_triangle_angle_degrees": (
                "post_cleaning_geometric_validation."
                "minimum_triangle_angle_degrees"
            ),
        },
        depends_on=("structural_cleaning",),
        validation=(
            "Every accepted chain is either eligible or quarantined with an "
            "explicit geometric reason; the configured thresholds recorded in "
            "every task summary agree."
        ),
        entry_point="scripts/pdbclean/finalize_geometric_validation.py",
    ),
    StageSpec(
        stage_id="complete_bri",
        canonical_stage="Stage 3-4",
        ordinal=6,
        title="Complete BRI",
        layer=LAYER_GOLD,
        purpose=(
            "Compute the complete backbone rigid invariant for every eligible "
            "chain. Complete BRI is the final geometric representation."
        ),
        directory=_protocol_dir("bri"),
        primary_output="finalized/bri.parquet",
        input_count_keys=("input_eligible_chain_count",),
        output_count_keys=("bri_chain_count",),
        scientific_parameters=(
            "bri.implementation",
            "bri.implementation_version",
            "bri.representation_precision_angstrom",
            "bri.vector_partition_key",
        ),
        compatibility={
            "snapshot": "snapshot.snapshot_id",
            "cleaning_protocol": "release.protocol_version",
        },
        depends_on=("geometric_validation",),
        validation=(
            "Chain accounting closes against the eligible population and the "
            "differential gate against pinned BRI v1.2.2 has passed."
        ),
        entry_point="scripts/pdbclean/finalize_bri.py",
    ),
    StageSpec(
        stage_id="brain",
        canonical_stage="Stage 5",
        ordinal=7,
        title="Brain (9-D average BRI)",
        layer=LAYER_GOLD,
        purpose=(
            "Compute the 9-dimensional average BRI vector. Brain is the "
            "filtering and indexing layer only; it never classifies duplicates."
        ),
        directory=_protocol_dir("brain"),
        primary_output="finalized/brain.parquet",
        input_count_keys=("input_bri_chain_count",),
        output_count_keys=("brain_chain_count",),
        scientific_parameters=("brain.dimension", "brain.definition"),
        compatibility={
            "snapshot": "snapshot.snapshot_id",
            "cleaning_protocol": "release.protocol_version",
        },
        depends_on=("complete_bri",),
        validation=(
            "Every chain with m >= 2 has a defined Brain vector; m = 1 chains "
            "are explicitly recorded as undefined."
        ),
        entry_point="pdbclean.brain_finalize_cli",
    ),
    StageSpec(
        stage_id="length_buckets",
        canonical_stage="Stage 6",
        ordinal=8,
        title="Exact chain-length grouping",
        layer=LAYER_GOLD,
        purpose=(
            "Partition chains into exact chain-length buckets. Comparison only "
            "ever happens inside one bucket, because BRI vectors of different "
            "lengths are not comparable."
        ),
        directory=_protocol_dir("length_buckets"),
        primary_output="finalized/bucket_index.parquet",
        input_count_keys=("input_chain_count",),
        output_count_keys=(
            "distinct_length_bucket_count",
            "length_bucket_count",
        ),
        scientific_parameters=("bri.vector_partition_key",),
        compatibility={
            "snapshot": "snapshot.snapshot_id",
            "cleaning_protocol": "release.protocol_version",
        },
        depends_on=("brain",),
        validation=(
            "Bucket membership partitions the population exactly once with no "
            "chain lost or duplicated."
        ),
        entry_point="pdbclean.length_buckets_cli",
    ),
    StageSpec(
        stage_id="candidate_filtering",
        canonical_stage="Stage 7",
        ordinal=9,
        title="Brain candidate filtering",
        layer=LAYER_GOLD,
        purpose=(
            "Lossless same-length Brain prefilter using scipy cKDTree plus an "
            "exact integer post-filter. Produces candidate pairs only."
        ),
        directory=_protocol_dir("brain_prefilter"),
        primary_output="finalized/candidates.parquet",
        input_count_keys=("input_chain_count", "brain_defined_chain_count"),
        output_count_keys=("candidate_pair_count",),
        scientific_parameters=(
            "brain_filter.threshold_angstrom",
            "brain_filter.operator",
            "brain_filter.engine",
            "brain_filter.exact_integer_post_filter",
        ),
        compatibility={
            "snapshot": "snapshot.snapshot_id",
            "cleaning_protocol": "release.protocol_version",
            "brain_threshold_angstrom": "brain_filter.threshold_angstrom",
        },
        depends_on=("length_buckets",),
        validation=(
            "No emitted pair exceeds the exact integer Brain radius, and the "
            "population accounting over buckets closes."
        ),
        entry_point="pdbclean.brain_prefilter_production",
    ),
    StageSpec(
        stage_id="complete_bri_nn",
        canonical_stage="Stage 8-9",
        ordinal=10,
        title="Complete-BRI nearest-neighbour search",
        layer=LAYER_GOLD,
        purpose=(
            "Exact complete-BRI L-infinity radius search over the Brain "
            "candidates using the compressed cover tree."
        ),
        directory=_protocol_dir("full_bri_nn"),
        primary_output="finalized/candidate_near_duplicates.parquet",
        input_count_keys=("brain_candidate_pair_count",),
        output_count_keys=("total_near_duplicate_count",),
        scientific_parameters=(
            "duplicate_search.near_duplicate_threshold_angstrom",
            "duplicate_search.operator",
            "duplicate_search.metric",
            "duplicate_search.search_engine",
            "duplicate_search.representation",
        ),
        compatibility={
            "snapshot": "snapshot.snapshot_id",
            "cleaning_protocol": "release.protocol_version",
            "near_duplicate_threshold_angstrom": (
                "duplicate_search.near_duplicate_threshold_angstrom"
            ),
        },
        depends_on=("candidate_filtering",),
        validation=(
            "Pair accounting closes, no emitted pair exceeds the radius, and "
            "no complete-BRI hit was missed by the Brain prefilter."
        ),
        entry_point="pdbclean.full_bri_nn_production",
    ),
    StageSpec(
        stage_id="duplicate_classification",
        canonical_stage="Stage 10",
        ordinal=11,
        title="Duplicate classification",
        layer=LAYER_GOLD,
        purpose=(
            "Classify each pair from its exact integer-milliangstrom "
            "complete-BRI distance into exact duplicate and near duplicate."
        ),
        directory=_protocol_dir("duplicate_classification"),
        primary_output="finalized/candidate_classifications.parquet",
        input_count_keys=("input_pair_count",),
        output_count_keys=("paper_near_duplicate_pair_count",),
        scientific_parameters=(
            "duplicate_search.near_duplicate_threshold_angstrom",
            "duplicate_search.operator",
            "duplicate_search.exact_duplicate_criterion",
            "duplicate_search.final_classification_basis",
        ),
        compatibility={
            "snapshot": "snapshot.snapshot_id",
            "cleaning_protocol": "release.protocol_version",
            "paper_near_duplicate_threshold_angstrom": (
                "duplicate_search.near_duplicate_threshold_angstrom"
            ),
        },
        depends_on=("complete_bri_nn",),
        validation=(
            "Every zero-distance pair is also a near duplicate and no "
            "classification contradicts the exact integer distance."
        ),
        entry_point="pdbclean.duplicate_classification_production",
    ),
    StageSpec(
        stage_id="downstream_metadata",
        canonical_stage="Stage 14 input",
        ordinal=12,
        title="Deposition metadata",
        layer=LAYER_GOLD,
        purpose=(
            "Snapshot-consistent experimental method and resolution metadata "
            "used by representative ranking. Not a duplicate criterion."
        ),
        directory=_protocol_dir("downstream_metadata"),
        primary_output="finalized/entry_metadata.parquet",
        output_count_keys=(
            "participating_deposition_count",
            "entry_count",
            "deposition_count",
        ),
        scientific_parameters=("snapshot.snapshot_id",),
        compatibility={
            "snapshot": "snapshot.snapshot_id",
            "cleaning_protocol": "release.protocol_version",
        },
        depends_on=("structural_cleaning",),
        validation=(
            "Metadata is drawn only from the resolved snapshot; no "
            "old-snapshot values are admitted."
        ),
        entry_point="pdbclean.downstream_metadata_finalize",
    ),
    StageSpec(
        stage_id="redundancy_graph",
        canonical_stage="Stage 14a",
        ordinal=13,
        title="Near-duplicate graph",
        layer=LAYER_GOLD,
        purpose=(
            "Build the near-duplicate graph over the classified pairs. "
            "Connected components are recorded but are NOT duplicate "
            "equivalence classes."
        ),
        directory=_protocol_dir("stage14_geometric_graph"),
        summary="global_summary.json",
        primary_output="edge_node_components.parquet",
        input_count_keys=("source_edge_count",),
        output_count_keys=("edge_component_count",),
        scientific_parameters=(
            "duplicate_search.near_duplicate_threshold_angstrom",
            "graph.connected_component_is_duplicate_equivalence",
            "graph.require_direct_edge_for_removal",
            "representative_selection.minimum_deduplicated_chain_length",
        ),
        compatibility={
            "threshold": "__near_duplicate_threshold_text__",
        },
        depends_on=("duplicate_classification",),
        validation=(
            "No edge exceeds the threshold, m = 1 edges are excluded from the "
            "deduplication graph, and no removal decision is taken here."
        ),
        entry_point="scripts/build_stage14_geometric_graph.py",
    ),
    StageSpec(
        stage_id="representative_selection",
        canonical_stage="Stage 14b",
        ordinal=14,
        title="Redundancy resolution",
        layer=LAYER_GOLD,
        purpose=(
            "Select a representative for each component under the COMP702 "
            "policy. Every removed chain must have its own direct edge to its "
            "assigned representative."
        ),
        directory=_protocol_dir("stage14_representative_selection_v1"),
        primary_output="representative_mapping.parquet",
        input_count_keys=("canonical_input_chain_count",),
        output_count_keys=("final_retained_chain_count",),
        scientific_parameters=(
            "representative_selection.policy_name",
            "representative_selection.policy_version",
            "representative_selection.ranking",
            "representative_selection.m1_policy",
            "representative_selection.nonclique_algorithm",
            "graph.require_direct_edge_for_removal",
            "graph.automatic_transitive_removal",
        ),
        compatibility={
            "duplicate_threshold": "__near_duplicate_threshold_text__",
            "policy_version": "representative_selection.policy_version",
        },
        depends_on=("redundancy_graph", "downstream_metadata"),
        validation=(
            "Every removed chain has a direct edge to its representative, no "
            "transitive removal occurred, and all m = 1 chains are retained."
        ),
        entry_point="scripts/select_stage14_representatives.py",
    ),
    StageSpec(
        stage_id="gold_release",
        canonical_stage="Stage 14c",
        ordinal=15,
        title="Final Gold release",
        layer=LAYER_GOLD,
        purpose=(
            "Publish the retained-chain manifest, removed-chain audit, "
            "representative mapping and release provenance."
        ),
        directory="{release_root}/{release_name}",
        summary="release_manifest.json",
        primary_output="data/retained_chains.parquet",
        input_count_keys=("canonical_input_chain_count",),
        output_count_keys=("retained_chain_count",),
        scientific_parameters=(
            "duplicate_search.near_duplicate_threshold_angstrom",
            "representative_selection.policy_name",
            "graph.require_direct_edge_for_removal",
        ),
        compatibility={
            "snapshot": "snapshot.snapshot_id",
            "protocol": "release.protocol_version",
            "near_duplicate_threshold": "__near_duplicate_threshold_text__",
        },
        depends_on=("representative_selection",),
        validation=(
            "Release artefact hashes match, counts reconcile with the "
            "representative mapping, and every configured expectation gate "
            "passes."
        ),
        entry_point="scripts/build_stage14_final_release.py",
    ),
)


STAGES_BY_ID: dict[str, StageSpec] = {stage.stage_id: stage for stage in STAGES}


def stage_order() -> list[str]:
    """Return stage ids in pipeline order."""

    return [stage.stage_id for stage in sorted(STAGES, key=lambda s: s.ordinal)]


def stages_for_layer(layer: str) -> list[StageSpec]:
    return [stage for stage in STAGES if stage.layer == layer]


def release_name(*, dataset_name: str, snapshot: str, protocol: str,
                 suffix: str = "dedup-v1") -> str:
    """Return the canonical release directory name.

    Matches the frozen naming used by ``run_stage14_final_release.sbatch``:
    ``PDBClean-<snapshot>-<protocol>-dedup-v1``.
    """

    return f"{dataset_name}-{snapshot}-{protocol}-{suffix}"


def resolve_directory(
    stage: StageSpec,
    *,
    output_root: str,
    snapshot: str,
    protocol: str,
    release_root: str,
    release: str,
) -> str | None:
    """Substitute path variables into a stage directory template."""

    if stage.directory is None:
        return None

    return stage.directory.format(
        output_root=output_root,
        snapshot=snapshot,
        protocol=protocol,
        release_root=release_root,
        release_name=release,
    )


def stage_catalogue() -> list[dict[str, Any]]:
    """Return a JSON-serialisable description of the whole pipeline."""

    return [
        {
            "stage_id": stage.stage_id,
            "ordinal": stage.ordinal,
            "canonical_stage": stage.canonical_stage,
            "title": stage.title,
            "layer": stage.layer,
            "purpose": stage.purpose,
            "persisted": stage.persisted,
            "depends_on": list(stage.depends_on),
            "scientific_parameters": list(stage.scientific_parameters),
            "validation": stage.validation,
            "entry_point": stage.entry_point,
        }
        for stage in sorted(STAGES, key=lambda s: s.ordinal)
    ]


# ---------------------------------------------------------------------------
# Canonical scientific timeline
# ---------------------------------------------------------------------------
#
# STAGES above is the *executable* registry: one entry per producer, ordered by
# execution. This is the *scientific* timeline: one entry per canonical
# identity, in canonical order, which is what every interface displays.
#
# The two differ deliberately:
#
#   * canonical Stage 3 and Stage 4 share one producer (`complete_bri`), because
#     `bri.compute_bri` ends with numpy.around(..., 3) -- the representation is
#     applied at the point of computation. Both identities are still shown.
#   * canonical Stage 8 and Stage 9 share one producer (`complete_bri_nn`), for
#     the same reason: the search emits its distances already represented.
#   * canonical Stages 11-13 are investigation and validation passes with no
#     producer on the release path. They appear with an honest status rather
#     than being omitted, so the timeline never looks like it jumps 10 -> 14.
#   * prerequisites are lettered, never numbered, so they cannot be mistaken
#     for scientific stages.
#
# Never sort this alphabetically. `canonical_timeline()` is the order.

ROLE_PREREQUISITE = "prerequisite"
ROLE_SCIENTIFIC = "scientific"
ROLE_INVESTIGATION = "investigation"
ROLE_VALIDATION = "validation"
ROLE_INPUT = "input"


@dataclass(frozen=True)
class CanonicalStage:
    """One canonical scientific identity in the pipeline timeline."""

    #: Stable key, e.g. ``"prerequisite_a"``, ``"stage_3"``, ``"stage_14a"``.
    key: str

    #: Display label, e.g. ``"Prerequisite A"``, ``"Stage 3"``, ``"Stage 14a"``.
    label: str

    #: Human title, e.g. ``"Complete BRI"``.
    title: str

    #: Position in canonical scientific order (1-based, contiguous).
    position: int

    role: str

    #: The executable registry entry that realises this identity, or ``None``
    #: when nothing on the release path produces it.
    producer: str | None = None

    #: Parent scientific stage for engineering subdivisions (14a/b/c -> 14).
    parent: str | None = None

    #: True when the producer is shared with another canonical identity.
    shared_producer: bool = False

    layer: str = LAYER_GOLD

    purpose: str = ""

    #: Where a non-orchestrated stage's frozen evidence lives, if any.
    frozen_output: str | None = None

    note: str = ""

    @property
    def display(self) -> str:
        """``"Stage 3 - Complete BRI"``: the label used across the whole UI."""

        return f"{self.label} — {self.title}"


CANONICAL_TIMELINE: tuple[CanonicalStage, ...] = (
    CanonicalStage(
        key="prerequisite_a",
        label="Prerequisite A",
        title="Snapshot resolution",
        position=1,
        role=ROLE_PREREQUISITE,
        producer="snapshot",
        layer=LAYER_SNAPSHOT,
        purpose=(
            "Resolve the requested snapshot to a concrete immutable identity "
            "before any processing begins."
        ),
    ),
    CanonicalStage(
        key="prerequisite_b",
        label="Prerequisite B",
        title="Bronze source manifest",
        position=2,
        role=ROLE_PREREQUISITE,
        producer="bronze_source_manifest",
        layer=LAYER_BRONZE,
        purpose=(
            "Immutable inventory of every source object in the snapshot. No "
            "scientific filtering."
        ),
    ),
    CanonicalStage(
        key="prerequisite_c",
        label="Prerequisite C",
        title="Silver parsed representation",
        position=3,
        role=ROLE_PREREQUISITE,
        producer="silver_parse",
        layer=LAYER_SILVER,
        purpose=(
            "Deterministic parsed representation, reconstructed on demand from "
            "the immutable Bronze object identity rather than persisted."
        ),
    ),
    CanonicalStage(
        key="stage_1",
        label="Stage 1",
        title="Structural cleaning (Protocol 3.2)",
        position=4,
        role=ROLE_SCIENTIFIC,
        producer="structural_cleaning",
        purpose=(
            "Apply quality rules Q001-Q006 and the terminal-residue exception."
        ),
    ),
    CanonicalStage(
        key="stage_2",
        label="Stage 2",
        title="Geometric validation",
        position=5,
        role=ROLE_SCIENTIFIC,
        producer="geometric_validation",
        purpose=(
            "Quarantine cleaned chains whose backbone geometry is degenerate, "
            "producing the canonical eligible population."
        ),
    ),
    CanonicalStage(
        key="stage_3",
        label="Stage 3",
        title="Complete BRI",
        position=6,
        role=ROLE_SCIENTIFIC,
        producer="complete_bri",
        shared_producer=True,
        purpose=(
            "Compute the complete backbone rigid invariant for each eligible "
            "chain. Complete BRI is the final geometric representation."
        ),
        note=(
            "Shares one producer with Stage 4: bri.compute_bri applies the "
            "precision grid at the point of computation."
        ),
    ),
    CanonicalStage(
        key="stage_4",
        label="Stage 4",
        title="BRI numerical representation",
        position=7,
        role=ROLE_SCIENTIFIC,
        producer="complete_bri",
        shared_producer=True,
        purpose=(
            "Represent complete BRI on the configured precision grid p. For "
            "the validated default p = 0.001 A this is "
            "BRI_units = round(BRI / p) = round(1000 * BRI), one unit = 1 mA."
        ),
        note=(
            "Shares one producer with Stage 3. Both canonical identities are "
            "shown because they are distinct scientific concepts."
        ),
    ),
    CanonicalStage(
        key="stage_5",
        label="Stage 5",
        title="Brain (9-D average BRI)",
        position=8,
        role=ROLE_SCIENTIFIC,
        producer="brain",
        purpose=(
            "Compute the 9-dimensional average BRI vector. Brain is the "
            "filtering and indexing layer only; it never classifies duplicates."
        ),
    ),
    CanonicalStage(
        key="stage_6",
        label="Stage 6",
        title="Exact chain-length grouping",
        position=9,
        role=ROLE_SCIENTIFIC,
        producer="length_buckets",
        purpose=(
            "Partition chains into exact chain-length buckets. BRI vectors of "
            "different lengths are not comparable."
        ),
    ),
    CanonicalStage(
        key="stage_7",
        label="Stage 7",
        title="Brain candidate filtering",
        position=10,
        role=ROLE_SCIENTIFIC,
        producer="candidate_filtering",
        purpose=(
            "Lossless same-length Brain prefilter (scipy cKDTree plus an exact "
            "integer post-filter). Produces candidate pairs only."
        ),
    ),
    CanonicalStage(
        key="stage_8",
        label="Stage 8",
        title="Complete-BRI nearest-neighbour search",
        position=11,
        role=ROLE_SCIENTIFIC,
        producer="complete_bri_nn",
        shared_producer=True,
        purpose=(
            "Exact complete-BRI L-infinity radius search over the Brain "
            "candidates, using the compressed cover tree."
        ),
        note="Shares one producer with Stage 9.",
    ),
    CanonicalStage(
        key="stage_9",
        label="Stage 9",
        title="Complete-BRI distance representation",
        position=12,
        role=ROLE_SCIENTIFIC,
        producer="complete_bri_nn",
        shared_producer=True,
        purpose=(
            "The authoritative distance representation: exact integer "
            "representation units on the configured precision grid."
        ),
        note=(
            "Shares one producer with Stage 8: the search emits its distances "
            "already represented."
        ),
    ),
    CanonicalStage(
        key="stage_10",
        label="Stage 10",
        title="Duplicate classification",
        position=13,
        role=ROLE_SCIENTIFIC,
        producer="duplicate_classification",
        purpose=(
            "Classify each tested pair from its exact complete-BRI L-infinity "
            "distance: exact when d == 0, near duplicate when d <= tau."
        ),
    ),
    CanonicalStage(
        key="stage_11",
        label="Stage 11",
        title="Acta-style downstream investigation",
        position=14,
        role=ROLE_INVESTIGATION,
        producer=None,
        frozen_output="acta_downstream_investigation_v2",
        purpose=(
            "Review pass over detected duplicates, feeding the Acta-style "
            "analysis."
        ),
        note=(
            "Not on the release path and NOT a deletion relation. Not executed "
            "by the orchestrated pipeline."
        ),
    ),
    CanonicalStage(
        key="stage_12",
        label="Stage 12",
        title="Scientific validation gates",
        position=15,
        role=ROLE_VALIDATION,
        producer=None,
        frozen_output="acta_manual_review_manifest_v2",
        purpose="Validation evidence supporting the review.",
        note=(
            "Not on the release path and NOT a deletion relation. Not executed "
            "by the orchestrated pipeline."
        ),
    ),
    CanonicalStage(
        key="stage_13",
        label="Stage 13",
        title="Detailed investigation / review",
        position=16,
        role=ROLE_INVESTIGATION,
        producer=None,
        frozen_output="acta_detailed_review_v2",
        purpose="Manual detailed review of a selected subset of pairs.",
        note=(
            "Explicitly NOT the global Stage-14 deletion relation, and never "
            "used as one. Not executed by the orchestrated pipeline."
        ),
    ),
    CanonicalStage(
        key="stage_14_input",
        label="Stage 14 input",
        title="Downstream entry metadata",
        position=17,
        role=ROLE_INPUT,
        producer="downstream_metadata",
        parent="Stage 14",
        purpose=(
            "Entry-level metadata used to rank candidate representatives. An "
            "input to Stage 14, not a scientific stage of its own."
        ),
    ),
    CanonicalStage(
        key="stage_14a",
        label="Stage 14a",
        title="Geometric redundancy graph",
        position=18,
        role=ROLE_SCIENTIFIC,
        producer="redundancy_graph",
        parent="Stage 14",
        purpose=(
            "Build the near-duplicate graph over detected pairs. A connected "
            "component is NOT a duplicate equivalence class."
        ),
    ),
    CanonicalStage(
        key="stage_14b",
        label="Stage 14b",
        title="Representative selection",
        position=19,
        role=ROLE_SCIENTIFIC,
        producer="representative_selection",
        parent="Stage 14",
        purpose=(
            "Deterministic quality-ordered greedy direct-edge cover. Every "
            "removed chain must have its own direct edge to its representative."
        ),
    ),
    CanonicalStage(
        key="stage_14c",
        label="Stage 14c",
        title="Final Gold release",
        position=20,
        role=ROLE_SCIENTIFIC,
        producer="gold_release",
        parent="Stage 14",
        purpose=(
            "Publish the retained-chain dataset, the removed-chain audit and "
            "the release manifest."
        ),
    ),
)


#: Producer stage_id -> the canonical identities it realises, in order.
CANONICAL_BY_PRODUCER: dict[str, list[CanonicalStage]] = {}

for _entry in CANONICAL_TIMELINE:
    if _entry.producer is not None:
        CANONICAL_BY_PRODUCER.setdefault(_entry.producer, []).append(_entry)

del _entry


def canonical_timeline() -> list[CanonicalStage]:
    """Return the canonical scientific timeline, in canonical order."""

    return sorted(CANONICAL_TIMELINE, key=lambda entry: entry.position)


def canonical_for_producer(stage_id: str) -> list[CanonicalStage]:
    """Return every canonical identity realised by one producer."""

    return list(CANONICAL_BY_PRODUCER.get(stage_id, []))


def canonical_display(stage_id: str) -> str:
    """Return the display label for a producer, e.g. ``"Stage 3-4 - ..."``."""

    entries = canonical_for_producer(stage_id)

    if not entries:
        return stage_id

    if len(entries) == 1:
        return entries[0].display

    labels = "/".join(entry.label for entry in entries)
    titles = " + ".join(entry.title for entry in entries)

    return f"{labels} — {titles}"


def canonical_catalogue() -> list[dict[str, Any]]:
    """JSON-serialisable canonical timeline for the UI and the CLI."""

    return [
        {
            "key": entry.key,
            "label": entry.label,
            "title": entry.title,
            "display": entry.display,
            "position": entry.position,
            "role": entry.role,
            "producer": entry.producer,
            "parent": entry.parent,
            "shared_producer": entry.shared_producer,
            "layer": entry.layer,
            "purpose": entry.purpose,
            "frozen_output": entry.frozen_output,
            "note": entry.note,
        }
        for entry in canonical_timeline()
    ]


#: Stages whose artefacts encode the BRI representation, and are therefore
#: only valid for the precision grid that produced them.
#:
#: The frozen stage summaries predate configurable precision and do not record
#: it, so a summary-key comparison cannot detect a mismatch.  The planner
#: instead treats an existing artefact as having been produced on the
#: implemented grid, which is true by construction: no other grid has ever
#: been executable.
PRECISION_DEPENDENT_STAGES: frozenset[str] = frozenset(
    ['complete_bri', 'brain', 'length_buckets', 'candidate_filtering', 'complete_bri_nn', 'duplicate_classification', 'redundancy_graph', 'representative_selection', 'gold_release']
)


def producer_canonical_label(stage_id: str) -> str:
    """Return the compound canonical label for one producer.

    A producer realising a single identity returns that label ("Stage 5").
    A producer shared by two adjacent identities returns the compound form
    ("Stage 3-4"), which is what ``StageSpec.canonical_stage`` records and what
    per-stage provenance stores.
    """

    entries = canonical_for_producer(stage_id)

    if not entries:
        return stage_id

    if len(entries) == 1:
        return entries[0].label

    numbers = [entry.label.replace("Stage ", "") for entry in entries]

    return f"Stage {numbers[0]}-{numbers[-1]}"
