"""Post-cleaning geometric validation for BRI-compatible protein chains.

This module does not perform Protocol 3.2 cleaning.  It validates the
geometry of an already-retained cleaned chain before BRI computation.

Scientific thresholds are explicit configuration values so they can be
changed without modifying the validation algorithm.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math

import numpy as np

from pdbclean.mmcif_parser import AtomObservation, ChainObservation
from pdbclean.quality import (
    protocol32_backbone_distances,
    protocol32_backbone_projection,
)


BACKBONE_ATOMS = ("N", "CA", "C")


def reconstruct_retained_backbone_chain(
    chain: ChainObservation,
    retained_label_seq_ids: list[int] | tuple[int, ...],
) -> ChainObservation:
    """Reconstruct exactly the Gold-retained Protocol 3.2 backbone.

    The source chain is projected with the same ATOM/N/CA/C semantics
    used during Protocol 3.2 cleaning. Only the exact residue IDs
    recorded in Gold are retained.

    Missing, duplicate, or reordered Gold lineage is rejected rather
    than silently repaired.
    """

    requested = tuple(retained_label_seq_ids)

    if not requested:
        raise ValueError(
            "retained_label_seq_ids must contain at least one residue"
        )

    if any(
        not isinstance(residue_id, int)
        or isinstance(residue_id, bool)
        for residue_id in requested
    ):
        raise ValueError(
            "retained_label_seq_ids must contain only integers"
        )

    if len(set(requested)) != len(requested):
        raise ValueError(
            "retained_label_seq_ids must not contain duplicates"
        )

    projected = protocol32_backbone_projection(
        chain,
        backbone_atoms=BACKBONE_ATOMS,
    )

    requested_set = set(requested)

    retained = replace(
        projected,
        atoms=[
            atom
            for atom in projected.atoms
            if atom.label_seq_id in requested_set
        ],
    )

    observed: list[int] = []
    seen: set[int] = set()

    for atom in retained.atoms:
        residue_id = atom.label_seq_id

        if residue_id is None:
            continue

        if residue_id not in seen:
            seen.add(residue_id)
            observed.append(residue_id)

    if tuple(observed) != requested:
        raise ValueError(
            "Reconstructed retained residue IDs do not match Gold lineage: "
            f"expected={requested!r}, observed={tuple(observed)!r}"
        )

    return retained


@dataclass(frozen=True)
class GeometricValidationConfig:
    """Scientific parameters for post-cleaning geometric validation."""

    minimum_backbone_distance_angstrom: float = 0.010
    minimum_triangle_angle_degrees: float = 3.0

    def __post_init__(self) -> None:
        if (
            not math.isfinite(self.minimum_backbone_distance_angstrom)
            or self.minimum_backbone_distance_angstrom < 0
        ):
            raise ValueError(
                "minimum_backbone_distance_angstrom must be finite and >= 0"
            )

        if (
            not math.isfinite(self.minimum_triangle_angle_degrees)
            or not 0 <= self.minimum_triangle_angle_degrees <= 180
        ):
            raise ValueError(
                "minimum_triangle_angle_degrees must be finite and in [0, 180]"
            )


@dataclass(frozen=True)
class GeometricViolation:
    """One explicitly recorded post-cleaning geometric violation."""

    violation_type: str
    residue_id: int | None
    details: str


@dataclass(frozen=True)
class GeometricValidationResult:
    """Validation result for one cleaned chain."""

    passed: bool
    minimum_observed_backbone_distance_angstrom: float | None
    minimum_observed_triangle_angle_degrees: float | None
    minimum_observed_basis_h_norm_angstrom: float | None
    violations: tuple[GeometricViolation, ...]


def _xyz(atom: AtomObservation) -> np.ndarray:
    """Return one atom coordinate as a float64 vector."""

    return np.asarray(
        [atom.x, atom.y, atom.z],
        dtype=np.float64,
    )


def _angle_degrees(
    point1: np.ndarray,
    vertex: np.ndarray,
    point2: np.ndarray,
) -> float:
    """Return the internal angle point1-vertex-point2 in degrees."""

    vector1 = point1 - vertex
    vector2 = point2 - vertex

    norm1 = float(np.linalg.norm(vector1))
    norm2 = float(np.linalg.norm(vector2))

    if (
        not math.isfinite(norm1)
        or not math.isfinite(norm2)
        or norm1 <= 0.0
        or norm2 <= 0.0
    ):
        return math.nan

    cosine = float(
        np.dot(vector1, vector2) / (norm1 * norm2)
    )

    if not math.isfinite(cosine):
        return math.nan

    # Protect arccos from tiny floating-point excursions outside [-1, 1].
    cosine = float(np.clip(cosine, -1.0, 1.0))

    return float(np.degrees(np.arccos(cosine)))


def _index_retained_backbone(
    chain: ChainObservation,
) -> dict[int, dict[str, AtomObservation]]:
    """Index exactly one N, CA and C observation for every retained residue."""

    backbone = protocol32_backbone_projection(
        chain,
        backbone_atoms=BACKBONE_ATOMS,
    )

    indexed: dict[
        int,
        dict[str, list[AtomObservation]],
    ] = {}

    for atom in backbone.atoms:
        if atom.label_seq_id is None:
            raise ValueError(
                "Post-cleaning geometric validation requires label_seq_id"
            )

        indexed.setdefault(
            atom.label_seq_id,
            {name: [] for name in BACKBONE_ATOMS},
        )

        if atom.atom_name in BACKBONE_ATOMS:
            indexed[atom.label_seq_id][atom.atom_name].append(atom)

    if not indexed:
        raise ValueError(
            "Post-cleaning geometric validation requires at least one residue"
        )

    result: dict[int, dict[str, AtomObservation]] = {}

    for residue_id in sorted(indexed):
        result[residue_id] = {}

        for atom_name in BACKBONE_ATOMS:
            observations = indexed[residue_id][atom_name]

            if len(observations) != 1:
                raise ValueError(
                    "Expected exactly one "
                    f"{atom_name} atom for residue {residue_id}; "
                    f"found {len(observations)}"
                )

            result[residue_id][atom_name] = observations[0]

    return result


def validate_post_cleaning_geometry(
    chain: ChainObservation,
    *,
    config: GeometricValidationConfig | None = None,
) -> GeometricValidationResult:
    """Validate one already-cleaned chain before BRI computation.

    Checks:

    1. retained residue IDs are consecutive;
    2. all N/CA/C coordinates are finite;
    3. the six pinned Q005 distances are at least the configured threshold;
    4. all three internal N-CA-C triangle angles satisfy the configured
       minimum angle;
    5. the Definition 3.4 residue basis is mathematically defined:
       |CA->N| > 0 and |h| > 0;
    6. no derived geometric quantity is non-finite.

    No residues are removed here.  Violations are recorded explicitly.
    """

    if config is None:
        config = GeometricValidationConfig()

    violations: list[GeometricViolation] = []

    try:
        indexed = _index_retained_backbone(chain)
    except ValueError as exc:
        return GeometricValidationResult(
            passed=False,
            minimum_observed_backbone_distance_angstrom=None,
            minimum_observed_triangle_angle_degrees=None,
            minimum_observed_basis_h_norm_angstrom=None,
            violations=(
                GeometricViolation(
                    violation_type="backbone_precondition_failed",
                    residue_id=None,
                    details=str(exc),
                ),
            ),
        )

    residue_ids = sorted(indexed)

    for left, right in zip(residue_ids, residue_ids[1:]):
        if right != left + 1:
            violations.append(
                GeometricViolation(
                    violation_type="nonconsecutive_residue_ids",
                    residue_id=left,
                    details=(
                        f"residue {left} is followed by residue {right}"
                    ),
                )
            )

    for residue_id, atoms in indexed.items():
        for atom_name, atom in atoms.items():
            coordinate = _xyz(atom)

            if not np.isfinite(coordinate).all():
                violations.append(
                    GeometricViolation(
                        violation_type="nonfinite_coordinate",
                        residue_id=residue_id,
                        details=f"atom={atom_name}",
                    )
                )

    observed_distances: list[float] = []

    try:
        q005_distances = protocol32_backbone_distances(
            chain,
            required_atoms=BACKBONE_ATOMS,
        )
    except ValueError as exc:
        violations.append(
            GeometricViolation(
                violation_type="q005_distance_precondition_failed",
                residue_id=None,
                details=str(exc),
            )
        )
        q005_distances = ()

    for item in q005_distances:
        distance = float(item.distance_angstrom)

        if math.isfinite(distance):
            observed_distances.append(distance)
        else:
            violations.append(
                GeometricViolation(
                    violation_type="nonfinite_backbone_distance",
                    residue_id=item.residue_id,
                    details=f"{item.atom1}-{item.atom2}",
                )
            )
            continue

        if distance < config.minimum_backbone_distance_angstrom:
            violations.append(
                GeometricViolation(
                    violation_type="backbone_distance_below_minimum",
                    residue_id=item.residue_id,
                    details=(
                        f"{item.atom1}-{item.atom2}: "
                        f"{distance:.12g} < "
                        f"{config.minimum_backbone_distance_angstrom:.12g}"
                    ),
                )
            )

    observed_angles: list[float] = []
    observed_h_norms: list[float] = []

    for residue_id in residue_ids:
        atoms = indexed[residue_id]

        n = _xyz(atoms["N"])
        ca = _xyz(atoms["CA"])
        c = _xyz(atoms["C"])

        triangle_angles = {
            "N": _angle_degrees(ca, n, c),
            "CA": _angle_degrees(n, ca, c),
            "C": _angle_degrees(n, c, ca),
        }

        for vertex, angle in triangle_angles.items():
            if not math.isfinite(angle):
                violations.append(
                    GeometricViolation(
                        violation_type="nonfinite_triangle_angle",
                        residue_id=residue_id,
                        details=f"vertex={vertex}",
                    )
                )
                continue

            observed_angles.append(angle)

            if angle < config.minimum_triangle_angle_degrees:
                violations.append(
                    GeometricViolation(
                        violation_type="triangle_angle_below_minimum",
                        residue_id=residue_id,
                        details=(
                            f"vertex={vertex}: "
                            f"{angle:.12g} < "
                            f"{config.minimum_triangle_angle_degrees:.12g}"
                        ),
                    )
                )

        ca_to_n = n - ca
        ca_to_c = c - ca

        ca_to_n_norm = float(np.linalg.norm(ca_to_n))

        if (
            not math.isfinite(ca_to_n_norm)
            or ca_to_n_norm <= 0.0
        ):
            violations.append(
                GeometricViolation(
                    violation_type="definition_3_4_undefined_ca_to_n",
                    residue_id=residue_id,
                    details=f"norm={ca_to_n_norm!r}",
                )
            )
            continue

        u = ca_to_n / ca_to_n_norm

        h = ca_to_c - np.dot(ca_to_c, u) * u
        h_norm = float(np.linalg.norm(h))

        if math.isfinite(h_norm):
            observed_h_norms.append(h_norm)

        if not math.isfinite(h_norm) or h_norm <= 0.0:
            violations.append(
                GeometricViolation(
                    violation_type="definition_3_4_undefined_h",
                    residue_id=residue_id,
                    details=f"norm={h_norm!r}",
                )
            )

    minimum_distance = (
        min(observed_distances)
        if observed_distances
        else None
    )

    minimum_angle = (
        min(observed_angles)
        if observed_angles
        else None
    )

    minimum_h_norm = (
        min(observed_h_norms)
        if observed_h_norms
        else None
    )

    return GeometricValidationResult(
        passed=not violations,
        minimum_observed_backbone_distance_angstrom=minimum_distance,
        minimum_observed_triangle_angle_degrees=minimum_angle,
        minimum_observed_basis_h_norm_angstrom=minimum_h_norm,
        violations=tuple(violations),
    )
