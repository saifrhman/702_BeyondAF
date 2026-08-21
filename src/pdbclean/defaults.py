"""Authoritative built-in scientific defaults for the COMP702 PDBClean pipeline.

This module is the single source of truth for the *validated* default value of
every scientific and dataset-version parameter the pipeline consumes.

Provenance of every value in :data:`VALIDATED_DEFAULTS`
------------------------------------------------------

Each default below reproduces the behaviour that actually produced the frozen
COMP702 release ``PDBClean-20260101-protocol3.2-comp702-v1-dedup-v1``.  Values
were taken from, in the project's authority order:

1. the executable production code
   (``brain_prefilter.py``, ``full_bri_nn_production.py``,
   ``duplicate_classification.py``, ``geometric_validation.py``);
2. the generated production summaries
   (``brain_prefilter/global_summary.json``,
   ``full_bri_nn/finalized/global_summary.json``,
   ``releases/.../provenance/stage8_full_bri_summary.json``,
   ``releases/.../release_manifest.json``);
3. the frozen configuration files
   ``config/pdbclean/protocol_3_2_comp702_v1.yaml`` (sha256 ``50bc5f98...``) and
   ``config/pdbclean/stage14_representative_policy_v1.yaml``
   (sha256 ``d4ca8bcb...``).

Both frozen YAML files are byte-immutable: their SHA256 digests are recorded
inside frozen provenance, so this module never edits them.  It restates their
values so that a run can be resolved without them, and
``tests/pdbclean/test_defaults.py`` asserts the restatement still agrees with
the files and with the frozen production summaries.

Deliberately NOT defaulted here
-------------------------------

``geometric_search`` in the frozen protocol YAML (``query_radius: 1.0``,
``retain_when: less_than 1.0``) describes the *historical PDB707K* design
documented in ``docs/pdbclean/development_status.md`` section 23-24 and
``docs/pdbclean/pipeline_spec.md`` sections 9-10.  The 2026 production pipeline
does not implement it: the executed Brain filter uses an inclusive integer
radius of ``tau_mA * (m - 1)`` and the final classification uses inclusive
``d_bri_mA <= 10``.  Executable code and generated manifests outrank
documentation, so the executed semantics are what this module encodes, under
the key :data:`brain_filter`.  See ``docs/CONFIGURATION.md``.
"""

from __future__ import annotations

import copy
from typing import Any


DEFAULTS_SCHEMA_NAME = "pdbclean_validated_defaults"
DEFAULTS_SCHEMA_VERSION = "1.0"

#: Version of the *default value set*.  Bump only when a default changes, which
#: is a scientific event and must be accompanied by an explicit decision.
DEFAULTS_VERSION = "1.0"

#: The BRI representation precision the executable pipeline implements.
#:
#: p = 0.001 A is currently *structural*, not parametric, in the production
#: code: ``bri.compute_bri`` ends with ``numpy.around(result, 3)`` to reproduce
#: the pinned BRI v1.2.2 canonicalisation, ``full_bri_compare`` scales by
#: ``BRI_MA_PER_ANGSTROM = 1000``, ``brain.compute_brain`` requires canonical
#: 3-decimal input, and ``full_bri_nn`` requires exact int64 milliangstroms.
#:
#: Precision is exposed as configuration so that a precision study is an
#: explicitly identified scientific configuration.  A run that requests any
#: other value is accepted as a *distinct* scientific configuration -- it gets
#: its own scientific hash and can never reuse the frozen release -- but the
#: production stages refuse to execute it, because implementing another grid
#: means changing the BRI canonicalisation, which is a scientific change.
IMPLEMENTED_REPRESENTATION_PRECISION_ANGSTROM = 0.001


#: The frozen COMP702 protocol configuration whose SHA256 is embedded in
#: production provenance.  Never modified; referenced for reproduction.
FROZEN_PROTOCOL_CONFIG_SHA256 = (
    "50bc5f98b26b3da7a2d57054c92235fde10afab7093201c89aeab2a194839bb5"
)

#: The frozen Stage-14 representative policy whose SHA256 is embedded in the
#: release manifest.  Never modified; referenced for reproduction.
FROZEN_REPRESENTATIVE_POLICY_SHA256 = (
    "d4ca8bcb1ed121bba6681a0de5ece6f5d9b8aff7ad6283e2512d7a2e762fd7c5"
)


VALIDATED_DEFAULTS: dict[str, Any] = {
    "schema_name": DEFAULTS_SCHEMA_NAME,
    "schema_version": DEFAULTS_SCHEMA_VERSION,
    "defaults_version": DEFAULTS_VERSION,
    # ------------------------------------------------------------------
    # Release identity
    # ------------------------------------------------------------------
    "release": {
        "dataset_name": "PDBClean",
        "protocol_version": "protocol3.2-comp702-v1",
    },
    # ------------------------------------------------------------------
    # Dataset-version parameters
    # ------------------------------------------------------------------
    "snapshot": {
        # Generic pipeline default.  The frozen COMP702 result is reproduced by
        # explicitly selecting mode=explicit with snapshot_id=20260101.
        "mode": "latest_complete",
        "snapshot_id": None,
        "bucket_url": "https://pdbsnapshots.s3.us-west-2.amazonaws.com",
    },
    # ------------------------------------------------------------------
    # Structural scope
    # ------------------------------------------------------------------
    "selection": {
        "models": {
            "policy": "first_model",
            "model_id": 1,
        },
        "chain_namespace": {
            "canonical": "label_asym_id",
            "preserve_author_chain_id": True,
            "preserve_entity_id": True,
        },
    },
    # ------------------------------------------------------------------
    # Protocol 3.2 structural cleaning rules (Q001 - Q006)
    # ------------------------------------------------------------------
    "quality_rules": {
        "entry_protein": {
            "rule_id": "Q001",
            "enabled": True,
            "scope": "entry",
            "source": "_entity_poly.type",
            "match": "startswith_polypeptide",
        },
        "disorder": {
            "rule_id": "Q002",
            "enabled": True,
            "atom_scope": ["N", "CA", "C"],
            "occupancy_semantics": "bri_v1.2.2_raw_token",
            "duplicate_backbone_atoms_are_disorder": True,
            "reject_alternate_locations_independently": False,
        },
        "residue_continuity": {
            "rule_id": "Q003",
            "enabled": True,
            "identifier": "label_seq_id",
        },
        "backbone_atoms": {
            "rule_id": "Q004",
            "enabled": True,
            "required_atoms": ["N", "CA", "C"],
        },
        "backbone_distance": {
            "rule_id": "Q005",
            "enabled": True,
            "minimum_distance_angstrom": 0.01,
        },
        "amino_acids": {
            "rule_id": "Q006",
            "enabled": True,
            "deposited_label_source": "_atom_site.label_comp_id",
            "mapping_source": "bri_v1.2.2_amino_acid_short",
            "accepted_mapped_codes": "canonical_20_one_letter",
        },
    },
    # ------------------------------------------------------------------
    # Post-cleaning geometric validity
    # ------------------------------------------------------------------
    "post_cleaning_geometric_validation": {
        "enabled": True,
        "minimum_triangle_angle_degrees": 3.0,
    },
    # ------------------------------------------------------------------
    # Complete BRI
    # ------------------------------------------------------------------
    "bri": {
        "enabled": True,
        "implementation": "backbone_rigid_invariant",
        "implementation_version": "1.2.2",
        "vector_partition_key": "chain_length",
        "require_full_accounting": True,
        # Representation precision p: the grid complete BRI is represented on.
        # This is NOT a duplicate threshold -- see `duplicate_search` for that.
        #
        # 0.001 A is the validated COMP702 value and the only value the
        # executable pipeline currently implements, because it is the
        # canonicalisation of the pinned BRI v1.2.2 reference implementation
        # (numpy.around(..., 3)).  It is exposed as configuration so that a
        # precision study is an explicit, separately identified scientific
        # configuration rather than a source edit -- see
        # `IMPLEMENTED_REPRESENTATION_PRECISION_ANGSTROM`.
        "representation_precision_angstrom": 0.001,
        # Generic quantisation rule: BRI_units = round(BRI / p).  For the
        # validated default p = 0.001 A this is exactly the frozen
        # `exact_representation` below, and one unit is one milliangstrom.
        "exact_representation": "BRI_mA = round(1000 * BRI)",
    },
    # ------------------------------------------------------------------
    # Brain (Definition 5.1): the filtering / indexing layer only
    # ------------------------------------------------------------------
    "brain": {
        "dimension": 9,
        "definition": "average_complete_BRI",
        "role": "filtering_and_indexing_only",
    },
    # Stage-7 Brain filter.  Executed semantics: within one exact chain-length
    # bucket the common denominator (m - 1) makes tau = 0.010 A exactly
    # tau_mA * (m - 1) integer milliangstrom-sum units, inclusive.
    "brain_filter": {
        "threshold_angstrom": 0.010,
        "threshold_milliangstrom": 10,
        "operator": "less_than_or_equal",
        "metric": "L_infinity",
        "engine": "scipy_cKDTree",
        "kdtree_p_norm": "inf",
        "kdtree_eps": 0.0,
        "exact_integer_post_filter": True,
        "grouping": "exact_chain_length",
    },
    # ------------------------------------------------------------------
    # Stage-8 complete-BRI nearest-neighbour search and Stage-10 classification
    # ------------------------------------------------------------------
    "duplicate_search": {
        # Key name matches the frozen protocol YAML so an existing config file
        # overrides it without translation.
        "near_duplicate_threshold_angstrom": 0.010,
        "operator": "less_than_or_equal",
        "metric": "L_infinity",
        "representation": "exact_integer_milliangstrom",
        "final_classification_basis": "complete_BRI",
        "exact_duplicate_criterion": "d_bri_mA == 0",
        "search_engine": "elkin_kurlin_compressed_cover_tree",
        "grouping": "exact_chain_length",
        "scientific_sequence": (
            "exact_length_then_brain_then_complete_bri_fast_nn"
        ),
    },
    # ------------------------------------------------------------------
    # Stage-14 near-duplicate graph and representative selection
    # ------------------------------------------------------------------
    "graph": {
        "build_connected_components": True,
        "calculate_edge_density": True,
        "calculate_clique_status": True,
        "require_direct_edge_for_removal": True,
        "connected_component_is_duplicate_equivalence": False,
        "automatic_transitive_removal": False,
    },
    "representative_selection": {
        "policy_name": "comp702_stage14_representative_selection",
        "policy_version": "1.0",
        "policy_config": (
            "config/pdbclean/stage14_representative_policy_v1.yaml"
        ),
        "minimum_deduplicated_chain_length": 2,
        "m1_policy": "retain_all",
        "ranking": [
            "terminal_trimmed_false_preferred",
            "lower_dirty_residue_count_preferred",
            "comparable_method_resolution_preferred",
            "canonical_chain_key_tiebreak",
        ],
        "nonclique_algorithm": (
            "deterministic_quality_ordered_greedy_direct_edge_cover"
        ),
        "use_stage13_review_subset_as_global_edge_set": False,
        "old_snapshot_comparison": False,
    },
    # ------------------------------------------------------------------
    # Infrastructure / execution parameters (not scientific)
    # ------------------------------------------------------------------
    "execution": {
        "batch_size": 500,
        "download_concurrency": 4,
        "quality_array_worker_count": 64,
        "quality_array_concurrency": 4,
        "max_retries": 3,
        "connection_timeout_seconds": 60,
        "atomic_writes": True,
        "write_success_markers": True,
        "executor": "dry-run",
    },
    "storage": {
        "temporary_root": "${TMPDIR}/pdbclean",
        "output_root": "outputs/pdbclean",
        "release_root": "outputs/releases",
        "run_root": "outputs/runs",
        # Durable, content-addressed snapshot preservation. Keeps a completed
        # run reproducible after any temporary cache expires.
        "durable_snapshot_root": "outputs/snapshot_store",
        # Disposable hot materialisation, optimised for parsing throughput.
        # May be deleted and rebuilt from the durable layer at any time.
        "hot_cache_root": "outputs/snapshot_cache",
        "retain_downloaded_mmcif": False,
    },
    "observability": {
        "structured_logging": True,
        "record_slurm_identifiers": True,
        "record_runtime": True,
        "record_peak_memory": True,
    },
    # ------------------------------------------------------------------
    # Dataset-version expectation gates.
    #
    # These are per-dataset acceptance counts, not scientific parameters.  They
    # are null by default: a generic snapshot has no known expected counts.  The
    # frozen 20260101 profile supplies the validated values so that reproducing
    # the frozen run keeps exactly the gates the frozen run enforced.
    # ------------------------------------------------------------------
    "expectations": {
        "canonical_input_chain_count": None,
        "mge2_node_count": None,
        "edge_count": None,
        "m1_edge_count": None,
        "retained_chain_count": None,
        "removed_chain_count": None,
        "total_near_duplicate_pair_count": None,
        "exact_duplicate_pair_count": None,
        "nonzero_near_duplicate_pair_count": None,
    },
}


def validated_defaults() -> dict[str, Any]:
    """Return a deep copy of the built-in validated defaults.

    A copy is returned so that callers can never mutate the authoritative
    default set in place.
    """

    return copy.deepcopy(VALIDATED_DEFAULTS)


def near_duplicate_threshold_milliangstrom(config: dict[str, Any]) -> int:
    """Return the exact integer-milliangstrom near-duplicate threshold.

    Complete BRI is represented exactly in integer milliangstroms, so the
    configured angstrom threshold must land exactly on the 0.001 A grid.  This
    mirrors the validation already performed by
    :func:`pdbclean.config.load_config` and by
    ``full_bri_nn_production._validate_threshold``.
    """

    section = config.get("duplicate_search")

    if not isinstance(section, dict):
        raise ValueError("duplicate_search section is missing")

    value = section.get("near_duplicate_threshold_angstrom")

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(
            "duplicate_search.near_duplicate_threshold_angstrom "
            "must be numeric"
        )

    value = float(value)

    if not value > 0.0 or value == float("inf"):
        raise ValueError(
            "duplicate_search.near_duplicate_threshold_angstrom "
            "must be finite and greater than zero"
        )

    scaled = value * 1000.0

    if abs(scaled - round(scaled)) > 1.0e-9:
        raise ValueError(
            "duplicate_search.near_duplicate_threshold_angstrom "
            "must be an exact integer milliangstrom value"
        )

    return int(round(scaled))


def brain_filter_threshold_milliangstrom(config: dict[str, Any]) -> int:
    """Return the exact integer-milliangstrom Brain filtering threshold.

    Brain coordinates are means of canonical 3-decimal-place BRI coordinates,
    so within one exact chain-length bucket the threshold is exact in
    milliangstrom-sum units.  The threshold must therefore also lie on the
    0.001 A grid.
    """

    section = config.get("brain_filter")

    if not isinstance(section, dict):
        raise ValueError("brain_filter section is missing")

    value = section.get("threshold_angstrom")

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(
            "brain_filter.threshold_angstrom must be numeric"
        )

    value = float(value)

    if not value > 0.0 or value == float("inf"):
        raise ValueError(
            "brain_filter.threshold_angstrom must be finite and "
            "greater than zero"
        )

    scaled = value * 1000.0

    if abs(scaled - round(scaled)) > 1.0e-9:
        raise ValueError(
            "brain_filter.threshold_angstrom must be an exact "
            "integer milliangstrom value"
        )

    return int(round(scaled))


# ---------------------------------------------------------------------------
# Representation precision (p) and precision-grid arithmetic
# ---------------------------------------------------------------------------
#
# The frozen pipeline works in exact integer milliangstroms because its
# validated precision is p = 0.001 A.  These helpers state the same arithmetic
# generically, in "representation units", so that configuration, validation and
# display do not assume one unit is always one milliangstrom.
#
# For the validated default p = 0.001 A:
#     1 representation unit == 1 mA
#     tau = 0.010 A         == 10 units


#: Tolerance for deciding that a ratio is an exact integer.  Matches the
#: tolerance the frozen threshold validators already use.
_GRID_TOLERANCE = 1.0e-9


def representation_precision_angstrom(config: dict[str, Any]) -> float:
    """Return the configured BRI representation precision p, in angstroms."""

    section = config.get("bri")

    if not isinstance(section, dict):
        raise ValueError("bri section is missing")

    value = section.get("representation_precision_angstrom")

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(
            "bri.representation_precision_angstrom must be numeric"
        )

    value = float(value)

    if not value > 0.0 or value == float("inf"):
        raise ValueError(
            "bri.representation_precision_angstrom must be finite and "
            "greater than zero"
        )

    return value


def quantise_angstrom_to_units(
    value_angstrom: float,
    precision_angstrom: float,
    *,
    name: str = "value",
) -> int:
    """Return ``value / p`` as an exact integer number of representation units.

    Raises when the value does not land exactly on the configured grid, so a
    threshold can never be silently rounded into a different threshold.
    """

    if isinstance(value_angstrom, bool) or not isinstance(
        value_angstrom, (int, float)
    ):
        raise ValueError(f"{name} must be numeric")

    value = float(value_angstrom)

    if not value > 0.0 or value == float("inf"):
        raise ValueError(f"{name} must be finite and greater than zero")

    ratio = value / float(precision_angstrom)
    units = round(ratio)

    if abs(ratio - units) > _GRID_TOLERANCE:
        raise ValueError(
            f"{name} ({value:g} A) is not exactly representable on the "
            f"configured {precision_angstrom:g} A precision grid: "
            f"{value:g} / {precision_angstrom:g} = {ratio:.6f}, which is not "
            "a whole number of representation units. Choose a threshold that "
            "is an exact multiple of the precision, or change the precision."
        )

    if units <= 0:
        raise ValueError(
            f"{name} must be at least one representation unit"
        )

    return int(units)


def near_duplicate_threshold_units(config: dict[str, Any]) -> int:
    """Return tau as an exact integer number of representation units."""

    section = config.get("duplicate_search")

    if not isinstance(section, dict):
        raise ValueError("duplicate_search section is missing")

    return quantise_angstrom_to_units(
        section.get("near_duplicate_threshold_angstrom"),
        representation_precision_angstrom(config),
        name="duplicate_search.near_duplicate_threshold_angstrom",
    )


def brain_filter_threshold_units(config: dict[str, Any]) -> int:
    """Return the Brain filter threshold in representation units."""

    section = config.get("brain_filter")

    if not isinstance(section, dict):
        raise ValueError("brain_filter section is missing")

    return quantise_angstrom_to_units(
        section.get("threshold_angstrom"),
        representation_precision_angstrom(config),
        name="brain_filter.threshold_angstrom",
    )


def representation_unit_label(precision_angstrom: float) -> str:
    """Return a human label for one representation unit.

    Only the validated precision has an established physical name, so any
    other grid is described generically rather than being given one.
    """

    if abs(float(precision_angstrom) - 0.001) <= _GRID_TOLERANCE:
        return "mA"

    return "unit"


def precision_is_implemented(config: dict[str, Any]) -> bool:
    """Whether the executable pipeline implements the configured precision."""

    return abs(
        representation_precision_angstrom(config)
        - IMPLEMENTED_REPRESENTATION_PRECISION_ANGSTROM
    ) <= _GRID_TOLERANCE


class PrecisionNotImplementedError(RuntimeError):
    """Raised when a stage is asked to run on an unimplemented precision grid."""


def require_implemented_precision(config: dict, *, stage: str) -> float:
    """Refuse to execute a stage under an unimplemented precision grid.

    ``bri.representation_precision_angstrom`` is exposed as configuration so a
    precision study is an explicitly identified scientific configuration.  The
    executable pipeline, however, implements only the validated grid: the BRI
    canonicalisation is ``numpy.around(..., 3)`` (pinned BRI v1.2.2) and every
    downstream representation assumes exact integer milliangstroms.

    Rather than silently producing v1.2.2-canonical output labelled with some
    other precision, a stage refuses to run.  Implementing another grid is a
    scientific change and needs an explicit decision.
    """

    precision = representation_precision_angstrom(config)

    if not precision_is_implemented(config):
        raise PrecisionNotImplementedError(
            f"{stage}: bri.representation_precision_angstrom is set to "
            f"{precision:g} A, but the executable pipeline implements only "
            f"{IMPLEMENTED_REPRESENTATION_PRECISION_ANGSTROM:g} A. That value "
            "is the pinned BRI v1.2.2 canonicalisation "
            "(numpy.around(..., 3)); representing complete BRI on a different "
            "grid would change the scientific method and is not implemented. "
            "The configuration is still valid and has its own scientific "
            "identity, but no stage will execute under it."
        )

    return precision
