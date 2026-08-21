"""Mol* scene construction, extracted from the validated preparation script.

Every function below is lifted **verbatim in behaviour** from
``reports/molstar_exact_duplicate_examples/prepare_scenes.py``, which produced
the frozen example scenes. That script executes at import (it writes scenes and
prints a report), so the reusable pure core is extracted here instead of being
imported from it. ``tests/pdbclean/test_molstar.py`` asserts that this module
still reproduces the frozen ``.mvsj`` scenes, so the two cannot drift.

No scientific decision is made here. Mol* is **visual inspection only**: a
scene never determines whether two chains are duplicates. Duplicate
classification comes from the complete-BRI L-infinity calculation alone, and
the scene is built *from* that recorded result.

The superposition here (Kabsch, on paired backbone atoms) exists solely to put
two structures in a comparable orientation on screen. It is not part of the
duplicate criterion and is never fed back into it.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import gemmi
import numpy as np


ATOM_ORDER = {
    "N": 0,
    "CA": 1,
    "C": 2,
}


# ----------------------------------------------------------------------
# mmCIF helpers
# ----------------------------------------------------------------------

def col(block, tag):
    values = block.find_values(tag)
    if not values:
        raise RuntimeError(f"Missing mmCIF column: {tag}")
    return list(values)


def atom_site_columns(path):
    doc = gemmi.cif.read(str(path))
    block = doc.sole_block()

    group = col(block, "_atom_site.group_PDB")
    label_chain = col(block, "_atom_site.label_asym_id")
    label_seq = col(block, "_atom_site.label_seq_id")
    label_atom = col(block, "_atom_site.label_atom_id")
    xs = col(block, "_atom_site.Cartn_x")
    ys = col(block, "_atom_site.Cartn_y")
    zs = col(block, "_atom_site.Cartn_z")

    model_col = block.find_values("_atom_site.pdbx_PDB_model_num")
    models = list(model_col) if model_col else ["1"] * len(group)

    return {
        "group": group,
        "chain": label_chain,
        "seq": label_seq,
        "atom": label_atom,
        "x": xs,
        "y": ys,
        "z": zs,
        "model": models,
    }


def backbone_coordinates(path, chain_id):
    """
    Extract model-1 ATOM N/CA/C coordinates using:
      label_asym_id
      label_seq_id
      label_atom_id

    This is deliberately the same label-chain namespace used by
    the Stage-10 duplicate-classification table.
    """
    a = atom_site_columns(path)

    rows = {}

    for i in range(len(a["group"])):
        if a["group"][i] != "ATOM":
            continue

        if a["model"][i] != "1":
            continue

        if a["chain"][i] != chain_id:
            continue

        atom = a["atom"][i]

        if atom not in ATOM_ORDER:
            continue

        seq = a["seq"][i]

        if seq in (".", "?"):
            continue

        key = (int(seq), atom)

        if key in rows:
            raise RuntimeError(
                f"Duplicate backbone atom in {path.name}: "
                f"chain={chain_id}, key={key}"
            )

        rows[key] = np.array(
            [
                float(a["x"][i]),
                float(a["y"][i]),
                float(a["z"][i]),
            ],
            dtype=np.float64,
        )

    keys = sorted(
        rows,
        key=lambda k: (k[0], ATOM_ORDER[k[1]])
    )

    if not keys:
        raise RuntimeError(
            f"No model-1 N/CA/C atoms found for "
            f"{path.name}, label_asym_id={chain_id}"
        )

    xyz = np.stack([rows[k] for k in keys])

    return keys, xyz


def selected_atom_coordinates(path, chain_ids):
    """
    Extract all model-1 protein ATOM coordinates for selected label chains.

    Used only to calculate a clean side-by-side presentation layout.
    It does NOT alter the molecular structures.
    """
    chain_ids = set(chain_ids)

    a = atom_site_columns(path)

    xyz = []

    for i in range(len(a["group"])):
        if a["group"][i] != "ATOM":
            continue

        if a["model"][i] != "1":
            continue

        if a["chain"][i] not in chain_ids:
            continue

        xyz.append(
            [
                float(a["x"][i]),
                float(a["y"][i]),
                float(a["z"][i]),
            ]
        )

    if not xyz:
        raise RuntimeError(
            f"No ATOM coordinates found in {path.name} "
            f"for chains {sorted(chain_ids)}"
        )

    return np.asarray(xyz, dtype=np.float64)


# ----------------------------------------------------------------------
# Exact chain pairing
# ----------------------------------------------------------------------

def paired_backbones(
    ref_path,
    ref_chain,
    mov_path,
    mov_chain,
    expected_residues,
):
    ref_keys, ref_xyz = backbone_coordinates(
        ref_path,
        ref_chain,
    )

    mov_keys, mov_xyz = backbone_coordinates(
        mov_path,
        mov_chain,
    )

    ref = dict(zip(ref_keys, ref_xyz))
    mov = dict(zip(mov_keys, mov_xyz))

    common_keys = sorted(
        set(ref) & set(mov),
        key=lambda k: (k[0], ATOM_ORDER[k[1]])
    )

    expected_atoms = expected_residues * 3

    if len(common_keys) != expected_atoms:
        raise RuntimeError(
            f"{ref_path.stem}:{ref_chain} vs "
            f"{mov_path.stem}:{mov_chain}: "
            f"expected {expected_atoms} common N/CA/C atoms "
            f"for {expected_residues} residues, "
            f"but found {len(common_keys)}"
        )

    A = np.stack([ref[k] for k in common_keys])
    B = np.stack([mov[k] for k in common_keys])

    return common_keys, A, B


# ----------------------------------------------------------------------
# Rigid alignment
# ----------------------------------------------------------------------

def kabsch_reference_from_moving(reference, moving):
    """
    Compute R,t such that:

        aligned_moving = (R @ moving.T).T + t

    optimally superposes `moving` onto `reference`.
    """
    ref_centroid = reference.mean(axis=0)
    mov_centroid = moving.mean(axis=0)

    Q = reference - ref_centroid
    P = moving - mov_centroid

    H = P.T @ Q

    U, _, Vt = np.linalg.svd(H)

    R = Vt.T @ U.T

    # Require a proper rotation, not reflection.
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    t = ref_centroid - R @ mov_centroid

    return R, t


def apply_transform(xyz, R, t):
    return (R @ xyz.T).T + t


def rmsd(A, B):
    d = A - B
    return float(
        np.sqrt(
            np.mean(
                np.sum(d * d, axis=1)
            )
        )
    )


def max_atom_distance(A, B):
    return float(
        np.max(
            np.linalg.norm(A - B, axis=1)
        )
    )


def mvs_matrix(R, t):
    """
    MolViewSpec 4x4 transform matrix.

    MVS requires column-major flattening.
    """
    M = np.eye(4, dtype=np.float64)

    M[:3, :3] = R
    M[:3, 3] = t

    return M.flatten(order="F").tolist()


# ----------------------------------------------------------------------
# Automatic side-by-side layout
# ----------------------------------------------------------------------

def side_by_side_translations(left_xyz, right_xyz, gap_A=30.0):
    """
    Place two structures symmetrically along world X.

    Both structures keep their original deposited orientation.
    Only translations are applied.

    Separation is derived from their X extents, so this does not
    depend on arbitrary fixed shift values.
    """
    left_centroid = left_xyz.mean(axis=0)
    right_centroid = right_xyz.mean(axis=0)

    left_width = float(
        left_xyz[:, 0].max() - left_xyz[:, 0].min()
    )

    right_width = float(
        right_xyz[:, 0].max() - right_xyz[:, 0].min()
    )

    centre_separation = (
        0.5 * left_width
        + 0.5 * right_width
        + gap_A
    )

    left_target = np.array(
        [-centre_separation / 2.0, 0.0, 0.0],
        dtype=np.float64,
    )

    right_target = np.array(
        [centre_separation / 2.0, 0.0, 0.0],
        dtype=np.float64,
    )

    left_translation = left_target - left_centroid
    right_translation = right_target - right_centroid

    return (
        left_translation.tolist(),
        right_translation.tolist(),
        centre_separation,
    )


# ----------------------------------------------------------------------
# MolViewSpec builders
# ----------------------------------------------------------------------

def representation(
    rep_type,
    color,
    opacity=1.0,
):
    children = [
        {
            "kind": "color",
            "params": {
                "color": color,
            },
        }
    ]

    if opacity < 1.0:
        children.append(
            {
                "kind": "opacity",
                "params": {
                    "opacity": opacity,
                },
            }
        )

    return {
        "kind": "representation",
        "params": {
            "type": rep_type,
        },
        "children": children,
    }


def structure_branch(
    url,
    selector,
    color,
    opacity=1.0,
    rep_type="cartoon",
    matrix=None,
    translation=None,
):
    """
    Build a valid MolViewSpec structure branch.

    IMPORTANT:
    transform and component are BOTH direct children of structure.

    Valid tree:

        structure
        ├── transform
        └── component
            └── representation
                ├── color
                └── opacity

    The previous version incorrectly nested component below transform,
    which violates the MolViewSpec node grammar.
    """

    component = {
        "kind": "component",
        "params": {
            "selector": selector,
        },
        "children": [
            representation(
                rep_type,
                color,
                opacity,
            )
        ],
    }

    if matrix is not None and translation is not None:
        raise ValueError(
            "Use either matrix or translation, not both."
        )

    structure_children = []

    if matrix is not None:
        structure_children.append(
            {
                "kind": "transform",
                "params": {
                    "matrix": matrix,
                },
            }
        )

    elif translation is not None:
        structure_children.append(
            {
                "kind": "transform",
                "params": {
                    "translation": translation,
                },
            }
        )

    structure_children.append(component)

    return {
        "kind": "download",
        "params": {
            "url": url,
        },
        "children": [
            {
                "kind": "parse",
                "params": {
                    "format": "mmcif",
                },
                "children": [
                    {
                        "kind": "structure",
                        "params": {
                            "type": "model",
                            "model_index": 0,
                        },
                        "children": structure_children,
                    }
                ],
            }
        ],
    }

def water_branch(
    url,
    color="#38bdf8",
    opacity=0.85,
):
    return structure_branch(
        url=url,
        selector="water",
        color=color,
        opacity=opacity,
        rep_type="spacefill",
    )


def scene(title, branches):
    return {
        "metadata": {
            "title": title,
            "version": "1",
            "timestamp": datetime.now(
                timezone.utc
            ).isoformat(),
        },
        "root": {
            "kind": "root",
            "children": branches + [
                {
                    "kind": "canvas",
                    "params": {
                        "background_color": "#ffffff",
                    },
                }
            ],
        },
    }


