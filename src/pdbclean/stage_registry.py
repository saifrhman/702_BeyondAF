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

    # -- detailed scientific description -------------------------------
    #
    # Written from the executable implementation and the frozen provenance.
    # Papers are cited for the scientific rationale only; a COMP702
    # computational choice is never attributed to a paper.

    #: What scientific problem this stage addresses.
    rationale: str = ""

    #: Exactly what operation is performed, in method terms.
    scientific_method: str = ""

    #: What arrives from upstream.
    stage_input: str = ""

    #: What the stage emits.
    stage_output: str = ""

    #: How the next stage consumes the result.
    downstream_role: str = ""

    #: Implementation notes: which code performs it, and which parts are
    #: COMP702 engineering rather than method.
    implementation_note: str = ""

    #: Compact method references, e.g. ``("anosova_match_2025",)``.
    references: tuple[str, ...] = ()

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
            "rationale": entry.rationale,
            "scientific_method": entry.scientific_method,
            "stage_input": entry.stage_input,
            "stage_output": entry.stage_output,
            "downstream_role": entry.downstream_role,
            "implementation_note": entry.implementation_note,
            "references": method_references(entry.references),
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


# ---------------------------------------------------------------------------
# Method references
# ---------------------------------------------------------------------------
#
# Cited for scientific rationale only. Where a step is a COMP702 computational
# choice rather than something a paper prescribes, the stage's
# `implementation_note` says so explicitly.

METHOD_REFERENCES: dict[str, dict[str, str]] = {
    "anosova_match_2025": {
        "key": "anosova_match_2025",
        "authors": "Anosova et al.",
        "title": (
            "A Complete and Bi-Continuous Invariant of Protein Backbones "
            "under Rigid Motion"
        ),
        "venue": "MATCH Communications in Mathematical and in Computer "
                 "Chemistry, 94(1), 97",
        "year": "2025",
        "doi": "10.46793/match.94-1.097A",
        "relevance": (
            "Defines the complete backbone rigid invariant and the average "
            "invariant, and establishes their completeness and continuity."
        ),
    },
    "wlodawer_acta_2025": {
        "key": "wlodawer_acta_2025",
        "authors": "Wlodawer et al.",
        "title": (
            "Duplicate entries in the Protein Data Bank: how to detect and "
            "handle them"
        ),
        "venue": "Acta Crystallographica Section D",
        "year": "2025",
        "doi": "10.1107/S2059798325001883",
        "relevance": (
            "Motivates duplicate detection in the PDB and characterises the "
            "categories of duplicate and near-duplicate deposition."
        ),
    },
}


def method_references(keys: tuple[str, ...]) -> list[dict[str, str]]:
    """Resolve reference keys to their bibliographic records."""

    return [
        METHOD_REFERENCES[key] for key in keys if key in METHOD_REFERENCES
    ]


# ---------------------------------------------------------------------------
# Detailed scientific descriptions
# ---------------------------------------------------------------------------
#
# Keyed by canonical stage key. Applied to CANONICAL_TIMELINE below, so the
# timeline stays the single source the UI, CLI and provenance all read.

_DESCRIPTIONS: dict[str, dict[str, Any]] = {
    "prerequisite_a": {
        "rationale": (
            "The PDB changes daily. A result is only reproducible if it names "
            "the exact archive state it was computed from, so a run must be "
            "bound to one immutable snapshot before any data is read."
        ),
        "scientific_method": (
            "The archive publishes dated snapshots. The configured selection "
            "mode is resolved once: 'latest_complete' discovers the newest "
            "complete snapshot and immediately pins it to its concrete "
            "YYYYMMDD identity; an explicit identity is used as given. The "
            "resolved identity, the mode that produced it, the bucket and a "
            "verified sample coordinate key are written to provenance. A "
            "resumed run reads that pinned identity and never re-resolves "
            "'latest', so it cannot silently drift to a newer archive state."
        ),
        "stage_input": "A selection mode, or an explicit snapshot identity.",
        "stage_output": (
            "A pinned snapshot identity in the resolved run configuration and "
            "in provenance. No data files."
        ),
        "downstream_role": (
            "Every later stage reads only from this snapshot, and the "
            "identity forms part of the run's scientific hash."
        ),
        "implementation_note": (
            "COMP702 orchestration. pdbclean.snapshot_selection."
        ),
    },
    "prerequisite_b": {
        "rationale": (
            "Before any scientific filtering, the run needs an immutable "
            "inventory of what the snapshot actually contained, so that every "
            "later decision can be traced back to a specific source object."
        ),
        "scientific_method": (
            "The snapshot's object listing is enumerated and recorded: PDB "
            "identifier, S3 source key, compressed byte size and ETag. This "
            "is an inventory only. No entry is excluded, no structure is "
            "parsed and no quality judgement is made at this layer."
        ),
        "stage_input": "The pinned snapshot identity.",
        "stage_output": (
            "source_manifest.parquet — one row per source object, with its "
            "verified identity."
        ),
        "downstream_role": (
            "Supplies the object identity that Silver parsing reconstructs "
            "from, and the provenance every Gold chain record carries."
        ),
        "implementation_note": (
            "COMP702 data engineering. scripts/pdbclean/create_manifest.py."
        ),
    },
    "prerequisite_c": {
        "rationale": (
            "Scientific filtering needs a deterministic in-memory "
            "representation of each entry, but storing a second full copy of "
            "the archive would double storage for no scientific gain."
        ),
        "scientific_method": (
            "Each mmCIF object is parsed deterministically by a pinned "
            "parser: models, chains in the canonical label_asym_id namespace, "
            "residues and atom sites, preserving author chain identifiers and "
            "entity identifiers alongside the canonical ones. The result is "
            "not persisted; it is reconstructed on demand from the immutable "
            "Bronze object identity."
        ),
        "stage_input": "Bronze source objects, by verified identity.",
        "stage_output": (
            "Nothing persisted, by design. Verified transitively: every Gold "
            "chain record carries the source key and ETag it was parsed from."
        ),
        "downstream_role": (
            "Supplies the deposited model/chain structures that structural "
            "cleaning applies its rules to."
        ),
        "implementation_note": (
            "COMP702 data engineering. pdbclean.mmcif_parser."
        ),
    },
    "stage_1": {
        "rationale": (
            "A geometric invariant is only meaningful for a backbone that is "
            "actually complete and unambiguous. Entries contain non-protein "
            "polymers, alternate conformations, gaps and missing backbone "
            "atoms, all of which would make a computed invariant describe "
            "something other than a single well-defined backbone."
        ),
        "scientific_method": (
            "Protocol 3.2 applies six rules, each recording an explicit "
            "accept/reject decision per chain:\\n\\n"
            "Q001 protein eligibility — the entry must be a polypeptide, "
            "determined from _entity_poly.type.\\n"
            "Q002 disorder — chains carrying alternate conformations or "
            "duplicated N/CA/C backbone atoms are rejected, because the "
            "backbone would not be uniquely defined.\\n"
            "Q003 residue continuity — the residue numbering (label_seq_id) "
            "must be continuous, since a gap breaks the consecutive-triple "
            "geometry the invariant is built from.\\n"
            "Q004 backbone atoms — every retained residue must carry all "
            "three of N, CA and C.\\n"
            "Q005 backbone distance — consecutive backbone atoms must be at "
            "least the configured minimum apart, rejecting coincident atoms.\\n"
            "Q006 amino acids — residue names must map to the canonical "
            "twenty via the pinned mapping.\\n\\n"
            "Model scope is fixed to a single deposited model, so that one "
            "chain contributes exactly one backbone."
        ),
        "stage_input": "Parsed deposited structures for the pinned snapshot.",
        "stage_output": (
            "quality/merged/accepted.parquet with the retained chains and "
            "their exact retained residue ranges, plus a rejected table "
            "recording which rule excluded each chain."
        ),
        "downstream_role": (
            "Defines the candidate population whose geometry is validated in "
            "Stage 2."
        ),
        "implementation_note": (
            "Rules and thresholds are COMP702 protocol decisions, pinned "
            "against the BRI v1.2.2 reference behaviour by a differential "
            "test. pdbclean.cleaning."
        ),
    },
    "stage_2": {
        "rationale": (
            "A chain can satisfy every bookkeeping rule and still be "
            "geometrically degenerate — collinear or coincident backbone "
            "atoms — which makes the local frames the invariant depends on "
            "numerically unstable or undefined."
        ),
        "scientific_method": (
            "Each accepted chain's backbone is checked for degeneracy. The "
            "N-CA-C triangle must have all angles above the configured "
            "minimum, and consecutive backbone atom distances must exceed the "
            "configured minimum. Chains failing either check are quarantined "
            "with an explicit recorded reason rather than being silently "
            "dropped or silently included."
        ),
        "stage_input": "Accepted chains from Stage 1.",
        "stage_output": (
            "finalized/eligible.parquet — the canonical eligible population — "
            "plus a quarantined table with per-chain geometric reasons."
        ),
        "downstream_role": (
            "The eligible population is the denominator for every later "
            "accounting check and the input to complete BRI."
        ),
        "implementation_note": (
            "Thresholds are COMP702 numerical-stability choices. "
            "pdbclean.geometric_validation."
        ),
    },
    "stage_3": {
        "rationale": (
            "Comparing backbones by coordinates is meaningless without "
            "removing the arbitrary rigid placement of each deposited "
            "structure. A complete invariant under rigid motion allows two "
            "backbones to be compared directly, with equality of the "
            "invariant corresponding to congruence of the backbones."
        ),
        "scientific_method": (
            "For a backbone of m residues, the complete backbone rigid "
            "invariant is an m x 9 matrix. Each row describes one residue's "
            "geometry in a local frame built from its neighbours, so the "
            "representation is invariant under rotation and translation of "
            "the whole structure while remaining complete: it determines the "
            "backbone up to rigid motion. Completeness is what makes it "
            "usable as the authoritative comparison basis, rather than a "
            "lossy descriptor that could conflate distinct backbones."
        ),
        "stage_input": "The canonical eligible chain population.",
        "stage_output": (
            "bri/finalized/bri.parquet — one m x 9 invariant per eligible "
            "chain."
        ),
        "downstream_role": (
            "Complete BRI is the authoritative geometric representation used "
            "for the final duplicate decision in Stages 8-10."
        ),
        "implementation_note": (
            "Definition from Anosova et al. Computed by pdbclean.bri against "
            "the pinned BRI v1.2.2 implementation, with a differential gate."
        ),
        "references": ("anosova_match_2025",),
    },
    "stage_4": {
        "rationale": (
            "Floating-point comparison of geometry is not reproducible across "
            "platforms and cannot support exact equality. Representing the "
            "invariant on a fixed grid makes distances exact integers, so "
            "'identical' and 'within threshold' are decidable rather than "
            "approximate."
        ),
        "scientific_method": (
            "The invariant is represented on a configured precision grid p:\\n\\n"
            "    BRI_units = round(BRI / p)\\n\\n"
            "At the validated default p = 0.001 A this is exactly "
            "round(1000 x BRI), and one representation unit is one "
            "milliangstrom. All later distance arithmetic is exact integer "
            "arithmetic in these units.\\n\\n"
            "Representation precision is NOT a duplicate threshold. p is how "
            "finely geometry is recorded; the near-duplicate threshold is a "
            "separate configured value compared against the resulting "
            "distance."
        ),
        "stage_input": "The complete BRI computed in Stage 3.",
        "stage_output": (
            "The same artefact as Stage 3: the representation is applied at "
            "the point of computation, and the exact integer conversion is "
            "verified lossless on load."
        ),
        "downstream_role": (
            "Every subsequent comparison, in Brain filtering and in the "
            "complete-BRI search, operates on these exact integer units."
        ),
        "implementation_note": (
            "The 0.001 A grid is a COMP702 computational choice matching the "
            "pinned BRI v1.2.2 canonicalisation (numpy.around(..., 3)); it is "
            "not a quantisation step prescribed by the paper. Integer "
            "conversion in pdbclean.full_bri_compare asserts losslessness."
        ),
    },
    "stage_5": {
        "rationale": (
            "Comparing every pair of chains by their full m x 9 invariants "
            "would be quadratic in the population and prohibitively "
            "expensive. A cheap summary that provably lower-bounds the full "
            "distance allows most non-matching pairs to be discarded without "
            "ever computing the expensive comparison."
        ),
        "scientific_method": (
            "Brain is the 9-dimensional average invariant: the column means "
            "of the complete BRI matrix, excluding its first row. It is a "
            "single 9-vector per chain regardless of chain length. Because it "
            "is an average of the same coordinates the full invariant is "
            "built from, a bound on the full distance implies a bound on the "
            "Brain distance, which is what makes it a safe filter.\\n\\n"
            "Brain is defined only for m >= 2. Chains of length 1 have no "
            "Brain vector and are handled separately.\\n\\n"
            "Brain is a filtering and indexing layer only. It never "
            "classifies duplicates."
        ),
        "stage_input": (
            "The represented complete BRI for every eligible chain."
        ),
        "stage_output": (
            "brain/finalized/brain.parquet — one 9-vector per chain with "
            "m >= 2, and an explicit record of the m = 1 chains as undefined."
        ),
        "downstream_role": (
            "Indexed in Stage 7 to generate candidate pairs cheaply."
        ),
        "implementation_note": (
            "Average invariant from Anosova et al. Computed from the "
            "canonical represented BRI by pdbclean.brain, which rejects input "
            "not on the representation grid."
        ),
        "references": ("anosova_match_2025",),
    },
    "stage_6": {
        "rationale": (
            "Two complete BRI matrices can only be compared elementwise if "
            "they have the same number of rows. Chains of different lengths "
            "have no correspondence between their invariants under this "
            "comparison, so comparing them is not merely expensive but "
            "undefined."
        ),
        "scientific_method": (
            "The eligible population is partitioned into buckets by exact "
            "retained residue count m. All later comparison happens strictly "
            "inside one bucket. The partition is verified to cover the "
            "population exactly once, with no chain lost or duplicated.\\n\\n"
            "Chains with m = 1 form their own bucket. They have no Brain "
            "vector and bypass Brain filtering; they are all retained by the "
            "deduplication policy."
        ),
        "stage_input": "The eligible population with its chain lengths.",
        "stage_output": (
            "finalized/bucket_index.parquet — bucket membership per chain, "
            "plus per-bucket population counts."
        ),
        "downstream_role": (
            "Scopes candidate generation: the prefilter is run once per "
            "bucket."
        ),
        "implementation_note": (
            "Correctness requirement, not an optimisation. "
            "pdbclean.length_buckets_cli."
        ),
    },
    "stage_7": {
        "rationale": (
            "Even within one length bucket, exhaustive pairwise comparison of "
            "complete BRI matrices is quadratic and dominated by pairs that "
            "are nowhere near each other. Candidate generation exists to "
            "discard those cheaply while provably keeping every pair that "
            "could still qualify."
        ),
        "scientific_method": (
            "Within each exact-length bucket, chains are indexed by their "
            "9-dimensional Brain vectors and queried under the L-infinity "
            "metric at the configured Brain threshold. Because the bucket "
            "shares a common denominator (m - 1), the threshold is expressed "
            "exactly as tau_units x (m - 1) in integer representation-unit "
            "sums, so the query radius is exact rather than floating point.\\n\\n"
            "The index query is followed by an exact integer post-filter, so "
            "the emitted candidate set contains all and only the pairs "
            "satisfying the Brain criterion.\\n\\n"
            "This stage produces CANDIDATES ONLY. No pair is classified here. "
            "The filter is lossless with respect to the final criterion: a "
            "pair discarded here could not have satisfied the complete-BRI "
            "threshold, which is what makes the pipeline's answer identical "
            "to exhaustive search while being tractable."
        ),
        "stage_input": (
            "Brain vectors and the exact chain-length buckets."
        ),
        "stage_output": (
            "finalized/candidates.parquet — candidate pairs, with per-bucket "
            "counts and the m = 1 bypass recorded."
        ),
        "downstream_role": (
            "The candidate set is the input population for the exact "
            "complete-BRI search in Stage 8."
        ),
        "implementation_note": (
            "SciPy cKDTree (p=inf, eps=0) plus an exact integer post-filter "
            "is a COMP702 engineering implementation of the Brain filtering "
            "step. cKDTree is used ONLY here, and is never the final "
            "complete-BRI search engine. pdbclean.brain_prefilter."
        ),
        "references": ("anosova_match_2025",),
    },
    "stage_8": {
        "rationale": (
            "The final duplicate decision must rest on complete BRI, the "
            "complete invariant, and not on the cheap Brain summary used for "
            "filtering. Complete BRI is the authoritative classification "
            "basis; this stage performs that comparison exactly, over the "
            "reduced candidate population."
        ),
        "scientific_method": (
            "For each candidate pair the complete BRI matrices are compared "
            "under the L-infinity metric: the distance is the maximum "
            "absolute difference over all m x 9 entries. The search is a "
            "radius query at the configured near-duplicate threshold, "
            "executed with a compressed cover tree, which returns exactly the "
            "pairs within the radius — it is an exact structure, not an "
            "approximate nearest-neighbour heuristic.\\n\\n"
            "Because complete BRI is a complete invariant, a distance of "
            "zero corresponds to congruent backbones."
        ),
        "stage_input": (
            "Candidate pairs from Stage 7 and the represented complete BRI."
        ),
        "stage_output": (
            "finalized/candidate_near_duplicates.parquet — qualifying pairs "
            "with their exact distances, plus the m = 1 pair table."
        ),
        "downstream_role": (
            "Supplies the distances that Stage 10 classifies."
        ),
        "implementation_note": (
            "Elkin-Kurlin compressed cover tree; the production search "
            "engine. COMP702 implementation choice for exact fast search. "
            "pdbclean.compressed_cover_tree, pdbclean.full_bri_nn_production."
        ),
        "references": ("anosova_match_2025", "wlodawer_acta_2025"),
    },
    "stage_9": {
        "rationale": (
            "The authoritative distance must be stored in a form that is "
            "exactly reproducible and exactly comparable, so that a "
            "classification can never depend on floating-point rounding."
        ),
        "scientific_method": (
            "Distances are emitted as exact integers in representation units "
            "on the configured precision grid. At the validated p = 0.001 A "
            "one unit is one milliangstrom, so a stored value of 10 is "
            "exactly 0.010 A. The angstrom value is derived from the integer, "
            "never the other way round.\\n\\n"
            "This stored integer distance is the authoritative complete-BRI "
            "L-infinity value for the pair."
        ),
        "stage_input": "The distances produced by the Stage 8 search.",
        "stage_output": (
            "The same artefact as Stage 8: the search emits its distances "
            "already represented."
        ),
        "downstream_role": (
            "Classification and the redundancy graph both compare these "
            "integers directly."
        ),
        "implementation_note": (
            "Exact integer representation is a COMP702 computational choice. "
            "Historical artefacts name the column d_bri_mA, correct for the "
            "0.001 A grid those runs used."
        ),
    },
    "stage_10": {
        "rationale": (
            "A distance alone does not say whether two depositions are "
            "duplicates. This stage applies the explicit, configured "
            "criterion that turns an exact distance into a recorded "
            "scientific classification."
        ),
        "scientific_method": (
            "Each tested pair is classified from its exact complete-BRI "
            "L-infinity distance d, in representation units:\\n\\n"
            "    d == 0            exact duplicate — congruent backbones\\n"
            "    0 < d <= tau      non-zero near duplicate\\n"
            "    d > tau           not a near duplicate\\n\\n"
            "The comparison is inclusive: a pair at exactly tau IS a near "
            "duplicate. Classification is always on complete BRI; the Brain "
            "distance plays no part in it.\\n\\n"
            "Pair counts are counts of PAIRS, not of chains. One chain can "
            "participate in many pairs, so the number of near-duplicate pairs "
            "is not the number of chains any policy would remove."
        ),
        "stage_input": "Complete-BRI distances for every tested pair.",
        "stage_output": (
            "finalized/candidate_classifications.parquet — per-pair flags for "
            "exact, near and non-zero-near, with population accounting."
        ),
        "downstream_role": (
            "Qualifying pairs become the edges of the Stage 14a redundancy "
            "graph. Investigation stages read the same table."
        ),
        "implementation_note": (
            "The inclusive threshold is a COMP702 configured decision. "
            "pdbclean.duplicate_classification."
        ),
        "references": ("wlodawer_acta_2025",),
    },
    "stage_11": {
        "rationale": (
            "Geometric duplication is a mathematical statement about "
            "backbones. Whether a duplicate pair is scientifically "
            "interesting — a re-deposition, a re-refinement, a related "
            "experiment — needs deposition metadata, not geometry."
        ),
        "scientific_method": (
            "Detected pairs are joined against deposition metadata "
            "(experimental method, resolution, deposition relationships) and "
            "reviewed in the style of the Acta duplicate analysis. The output "
            "characterises the detected population; it does not decide which "
            "chains are removed."
        ),
        "stage_input": "Classified pairs and downstream entry metadata.",
        "stage_output": (
            "acta_downstream_investigation_v2 — review tables and summaries."
        ),
        "downstream_role": (
            "Feeds the written analysis. Explicitly NOT an input to the "
            "Stage-14 deletion relation."
        ),
        "implementation_note": (
            "Investigation pass, not orchestrated on the release path."
        ),
        "references": ("wlodawer_acta_2025",),
    },
    "stage_12": {
        "rationale": (
            "A pipeline that reports duplicates must be shown to be correct. "
            "These gates provide the evidence that the fast path agrees with "
            "the exhaustive one and that the safety properties actually hold."
        ),
        "scientific_method": (
            "The validation evidence covers: candidate safety, that the Brain "
            "prefilter discards no pair that could satisfy the complete-BRI "
            "criterion; threshold boundary behaviour, that pairs exactly at "
            "the threshold are included; fast-versus-oracle agreement, that "
            "the cover-tree search returns the same pairs as brute force on "
            "test populations; direct-edge semantics, that every removal is "
            "justified by its own qualifying edge; and that connectedness is "
            "never treated as duplicate equivalence."
        ),
        "stage_input": "Outputs of the search, classification and selection.",
        "stage_output": (
            "acta_manual_review_manifest_v2 — validation manifests and "
            "verdicts."
        ),
        "downstream_role": (
            "Supports publication of the release. Not a deletion relation."
        ),
        "implementation_note": (
            "Validation pass, not orchestrated on the release path."
        ),
    },
    "stage_13": {
        "rationale": (
            "Automatic classification says two backbones are congruent. It "
            "does not say what happened scientifically. A manually reviewed "
            "subset gives the qualitative account behind the counts."
        ),
        "scientific_method": (
            "A selected subset of deposition pairs is reviewed by hand and "
            "categorised — exact congruence, non-zero near duplication, m = 1 "
            "cases, and cases resolved by examining the depositions "
            "themselves. The review is deliberately a subset, chosen for "
            "scientific interest."
        ),
        "stage_input": "Classified pairs and deposition metadata.",
        "stage_output": "acta_detailed_review_v2 — the reviewed subset.",
        "downstream_role": (
            "Feeds the written analysis. This subset is explicitly NOT the "
            "global Stage-14 deletion relation and is never used as one."
        ),
        "implementation_note": (
            "Manual investigation pass, not orchestrated."
        ),
        "references": ("wlodawer_acta_2025",),
    },
    "stage_14_input": {
        "rationale": (
            "Choosing which chain of a duplicate group to keep requires "
            "deposition-level information that geometry does not provide."
        ),
        "scientific_method": (
            "Entry-level metadata — experimental method, resolution and "
            "related deposition attributes — is assembled per entry so that "
            "candidate representatives can be ranked deterministically."
        ),
        "stage_input": "The snapshot's entry metadata.",
        "stage_output": "finalized/entry_metadata.parquet.",
        "downstream_role": (
            "Supplies the quality ordering used by representative selection."
        ),
        "implementation_note": (
            "An input to Stage 14, not a scientific stage of its own."
        ),
    },
    "stage_14a": {
        "rationale": (
            "Redundancy is a relation over many chains, not a property of a "
            "single pair. Representing it as a graph makes the structure of "
            "the redundancy explicit and auditable."
        ),
        "scientific_method": (
            "An undirected graph is built over chains with m >= 2. Each "
            "qualifying complete-BRI near-duplicate pair contributes exactly "
            "one edge; every edge is asserted to satisfy the threshold and "
            "the minimum chain length. Connected components are computed and "
            "characterised, including whether each component is a clique.\\n\\n"
            "A connected component is NOT a duplicate equivalence class. "
            "Near-duplication is not transitive: two chains joined by a path "
            "may be far apart, so component membership alone never justifies "
            "removing a chain."
        ),
        "stage_input": (
            "Qualifying near-duplicate pairs and the m = 1 pair table."
        ),
        "stage_output": (
            "stage14_geometric_graph — component summary, per-node component "
            "assignment, and a global summary recording the threshold used."
        ),
        "downstream_role": (
            "The graph is the structure representative selection walks."
        ),
        "implementation_note": (
            "COMP702 redundancy-resolution design. "
            "scripts/build_stage14_geometric_graph.py."
        ),
    },
    "stage_14b": {
        "rationale": (
            "Given a group of mutually redundant chains, exactly one should "
            "be retained — but only where the geometric evidence actually "
            "justifies removing the others. A naive 'keep one per component' "
            "rule would delete chains that were never shown to be duplicates "
            "of the chain that survived."
        ),
        "scientific_method": (
            "Within each component, candidate representatives are ranked "
            "deterministically by quality using entry metadata, with a fixed "
            "tie-break so the outcome does not depend on input order.\\n\\n"
            "A clique component is fully connected, so one representative "
            "covers every other member directly.\\n\\n"
            "A non-clique component is covered by a greedy DIRECT-EDGE cover: "
            "a chain may only be removed in favour of a representative it "
            "shares an actual qualifying edge with. Where one representative "
            "cannot cover the whole component, additional representatives are "
            "retained. This is why non-clique components can legitimately "
            "keep more than one chain.\\n\\n"
            "Two invariants are enforced and audited: every removed chain has "
            "its own direct edge, within threshold, to the chain it was "
            "assigned to; and there is NO transitive removal. Chains with "
            "m = 1 are retained unconditionally."
        ),
        "stage_input": (
            "The redundancy graph, qualifying edges, the accepted chain "
            "population and entry metadata."
        ),
        "stage_output": (
            "representative_mapping.parquet — the per-chain decision, its "
            "assigned representative and the direct edge distance justifying "
            "it — plus representatives.parquet and a summary asserting the "
            "direct-edge property held for every removal."
        ),
        "downstream_role": (
            "Determines the retained population that the Gold release "
            "publishes."
        ),
        "implementation_note": (
            "COMP702 representative policy v1, a project decision rather than "
            "a published algorithm. scripts/select_stage14_representatives.py."
        ),
        "references": ("wlodawer_acta_2025",),
    },
    "stage_14c": {
        "rationale": (
            "A result is only usable if it is published as an immutable, "
            "fully audited artefact whose every removal can be re-checked "
            "independently."
        ),
        "scientific_method": (
            "The retained chain dataset is written together with the "
            "removed-chain audit, the representative mapping and the edge "
            "tables. The release manifest records the population counts, the "
            "threshold and criterion used, the model scope, the policy "
            "version and hash, and explicit statements of the semantics that "
            "were NOT applied: no transitive removal, connectedness not "
            "treated as equivalence, no m = 1 deduplication, and the Stage-13 "
            "review subset not used as a global edge set.\\n\\n"
            "Publication happens only after every upstream validation gate "
            "has passed. Expectation gates, where configured for a known "
            "snapshot, must match exactly."
        ),
        "stage_input": (
            "The representative mapping and the full chain population."
        ),
        "stage_output": (
            "An immutable release: data/retained_chains.parquet, "
            "audit/removed_chain_audit.parquet, the representative mapping "
            "and edge tables, release_manifest.json and _SUCCESS, all with "
            "recorded SHA256 digests."
        ),
        "downstream_role": (
            "The published deduplicated dataset used by downstream work."
        ),
        "implementation_note": (
            "COMP702 release engineering. "
            "scripts/build_stage14_final_release.py."
        ),
    },
}


def _apply_descriptions() -> tuple[CanonicalStage, ...]:
    """Attach the detailed descriptions to the canonical timeline."""

    import dataclasses

    enriched = []

    for entry in CANONICAL_TIMELINE:
        detail = _DESCRIPTIONS.get(entry.key)

        enriched.append(
            dataclasses.replace(entry, **detail) if detail else entry
        )

    return tuple(enriched)


CANONICAL_TIMELINE = _apply_descriptions()

CANONICAL_BY_PRODUCER = {}

for _entry in CANONICAL_TIMELINE:
    if _entry.producer is not None:
        CANONICAL_BY_PRODUCER.setdefault(_entry.producer, []).append(_entry)

del _entry
