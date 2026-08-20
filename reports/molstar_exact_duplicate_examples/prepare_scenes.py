#!/usr/bin/env python3

import json
from datetime import datetime, timezone
from pathlib import Path

import gemmi
import numpy as np


HERE = Path(__file__).resolve().parent

FILES = {
    "1a0t": HERE / "1a0t.cif",
    "1oh2": HERE / "1oh2.cif",
    "9eur": HERE / "9eur.cif",
    "9ewn": HERE / "9ewn.cif",
}

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


def save(name, data):
    path = HERE / name

    path.write_text(
        json.dumps(
            data,
            indent=2,
        )
        + "\n"
    )

    print(f"WROTE {path.name}")


# ----------------------------------------------------------------------
# Verify inputs exist
# ----------------------------------------------------------------------

for pdb_id, path in FILES.items():
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing visualization input: {path}"
        )


# ----------------------------------------------------------------------
# Scientifically relevant pair calculations
# ----------------------------------------------------------------------

# Pair 1: chain B is exact BRI but different raw coordinates.
keys_B, A_B, B_B = paired_backbones(
    FILES["1a0t"],
    "B",
    FILES["1oh2"],
    "B",
    expected_residues=413,
)

R_B, t_B = kabsch_reference_from_moving(
    A_B,
    B_B,
)

B_B_aligned = apply_transform(
    B_B,
    R_B,
    t_B,
)

pair1_raw_rmsd = rmsd(
    A_B,
    B_B,
)

pair1_raw_max = max_atom_distance(
    A_B,
    B_B,
)

pair1_aligned_rmsd = rmsd(
    A_B,
    B_B_aligned,
)

pair1_aligned_max = max_atom_distance(
    A_B,
    B_B_aligned,
)

matrix_B = mvs_matrix(
    R_B,
    t_B,
)


# Pair 2: raw-coordinate-identical control.
keys_9, A_9, B_9 = paired_backbones(
    FILES["9eur"],
    "A",
    FILES["9ewn"],
    "A",
    expected_residues=300,
)

pair2_raw_rmsd = rmsd(
    A_9,
    B_9,
)

pair2_raw_max = max_atom_distance(
    A_9,
    B_9,
)


# ----------------------------------------------------------------------
# Automatic presentation layout
# ----------------------------------------------------------------------

pair1_left_xyz = selected_atom_coordinates(
    FILES["1a0t"],
    ["A", "B", "C"],
)

pair1_right_xyz = selected_atom_coordinates(
    FILES["1oh2"],
    ["A", "B", "C"],
)

(
    pair1_left_translation,
    pair1_right_translation,
    pair1_side_separation,
) = side_by_side_translations(
    pair1_left_xyz,
    pair1_right_xyz,
    gap_A=35.0,
)


pair2_left_xyz = selected_atom_coordinates(
    FILES["9eur"],
    ["A"],
)

pair2_right_xyz = selected_atom_coordinates(
    FILES["9ewn"],
    ["A"],
)

(
    pair2_left_translation,
    pair2_right_translation,
    pair2_side_separation,
) = side_by_side_translations(
    pair2_left_xyz,
    pair2_right_xyz,
    gap_A=30.0,
)


# ----------------------------------------------------------------------
# Save metrics/provenance
# ----------------------------------------------------------------------

metrics = {
    "pair1": {
        "reference": "1a0t:B",
        "moving": "1oh2:B",
        "retained_residue_count": 413,
        "backbone_atom_count": len(keys_B),
        "stage10_d_brain": 0.0,
        "stage10_d_bri_mA": 0,
        "stage10_d_bri_A": 0.0,
        "stage10_is_zero_duplicate": True,
        "raw_backbone_rmsd_A": pair1_raw_rmsd,
        "raw_backbone_max_distance_A": pair1_raw_max,
        "aligned_backbone_rmsd_A": pair1_aligned_rmsd,
        "aligned_backbone_max_distance_A": pair1_aligned_max,
        "transform_1oh2B_to_1a0tB_column_major": matrix_B,
        "side_by_side_centre_separation_A": pair1_side_separation,
        "side_by_side_left_translation": pair1_left_translation,
        "side_by_side_right_translation": pair1_right_translation,
    },

    "pair2": {
        "reference": "9eur:A",
        "moving": "9ewn:A",
        "retained_residue_count": 300,
        "backbone_atom_count": len(keys_9),
        "stage10_d_brain": 0.0,
        "stage10_d_bri_mA": 0,
        "stage10_d_bri_A": 0.0,
        "stage10_is_zero_duplicate": True,
        "raw_backbone_rmsd_A": pair2_raw_rmsd,
        "raw_backbone_max_distance_A": pair2_raw_max,
        "side_by_side_centre_separation_A": pair2_side_separation,
        "side_by_side_left_translation": pair2_left_translation,
        "side_by_side_right_translation": pair2_right_translation,
    },
}

save(
    "metrics.json",
    metrics,
)


# ----------------------------------------------------------------------
# Visual definitions
# ----------------------------------------------------------------------

BLUE = "#2563eb"
ORANGE = "#f97316"

ABC = [
    {
        "label_asym_id": "A",
    },
    {
        "label_asym_id": "B",
    },
    {
        "label_asym_id": "C",
    },
]

B_BRI = [
    {
        "label_asym_id": "B",
        "label_atom_id": "N",
    },
    {
        "label_asym_id": "B",
        "label_atom_id": "CA",
    },
    {
        "label_asym_id": "B",
        "label_atom_id": "C",
    },
]

A_CHAIN = {
    "label_asym_id": "A",
}


# ======================================================================
# PAIR 1
# ======================================================================

# ----------------------------------------------------------------------
# Scene 1: side-by-side presentation view
#
# No alignment is performed.
# Only display translations are applied.
# ----------------------------------------------------------------------

save(
    "pair1_side_by_side.mvsj",
    scene(
        "1a0t vs 1oh2 — side by side",
        [
            structure_branch(
                "./1a0t.cif",
                ABC,
                BLUE,
                1.0,
                "cartoon",
                translation=pair1_left_translation,
            ),

            structure_branch(
                "./1oh2.cif",
                ABC,
                ORANGE,
                1.0,
                "cartoon",
                translation=pair1_right_translation,
            ),
        ],
    ),
)


# ----------------------------------------------------------------------
# Scene 2: deposited coordinates
#
# Absolutely no transformation.
# A/C raw coordinates coincide.
# B raw coordinates differ.
# ----------------------------------------------------------------------

save(
    "pair1_deposited.mvsj",
    scene(
        "1a0t vs 1oh2 — deposited coordinates",
        [
            structure_branch(
                "./1a0t.cif",
                ABC,
                BLUE,
                1.0,
                "cartoon",
            ),

            structure_branch(
                "./1oh2.cif",
                ABC,
                ORANGE,
                0.58,
                "cartoon",
            ),
        ],
    ),
)


# ----------------------------------------------------------------------
# Scene 3: exact BRI chain B rigid superposition
#
# Only the N/CA/C atoms used by BRI are shown.
# 1oh2:B is rigidly transformed onto 1a0t:B.
# ----------------------------------------------------------------------

save(
    "pair1_chainB_superposed.mvsj",
    scene(
        "1a0t:B vs 1oh2:B — N/CA/C rigid superposition",
        [
            structure_branch(
                "./1a0t.cif",
                B_BRI,
                BLUE,
                1.0,
                "ball_and_stick",
            ),

            structure_branch(
                "./1oh2.cif",
                B_BRI,
                ORANGE,
                0.58,
                "ball_and_stick",
                matrix=matrix_B,
            ),
        ],
    ),
)


# ----------------------------------------------------------------------
# Scene 4: water discrepancy
#
# Protein polymers remain in deposited coordinates.
# 1a0t waters are shown.
# The finalized evidence records:
#   1a0t = 330 water atoms
#   1oh2 = 0 water atoms
# ----------------------------------------------------------------------

save(
    "pair1_waters.mvsj",
    scene(
        "1a0t vs 1oh2 — water discrepancy",
        [
            structure_branch(
                "./1a0t.cif",
                ABC,
                BLUE,
                0.48,
                "cartoon",
            ),

            structure_branch(
                "./1oh2.cif",
                ABC,
                ORANGE,
                0.38,
                "cartoon",
            ),

            water_branch(
                "./1a0t.cif",
            ),
        ],
    ),
)


# ======================================================================
# PAIR 2
# ======================================================================

# ----------------------------------------------------------------------
# Scene 5: 9eur / 9ewn side-by-side presentation view
# ----------------------------------------------------------------------

save(
    "pair2_side_by_side.mvsj",
    scene(
        "9eur:A vs 9ewn:A — side by side",
        [
            structure_branch(
                "./9eur.cif",
                A_CHAIN,
                BLUE,
                1.0,
                "cartoon",
                translation=pair2_left_translation,
            ),

            structure_branch(
                "./9ewn.cif",
                A_CHAIN,
                ORANGE,
                1.0,
                "cartoon",
                translation=pair2_right_translation,
            ),
        ],
    ),
)


# ----------------------------------------------------------------------
# Scene 6: 9eur / 9ewn raw-coordinate overlay
#
# No alignment at all.
# Their model-1 chain-A backbone coordinates are already identical.
# ----------------------------------------------------------------------

save(
    "pair2_exact_overlay.mvsj",
    scene(
        "9eur:A vs 9ewn:A — deposited-coordinate overlay",
        [
            structure_branch(
                "./9eur.cif",
                A_CHAIN,
                BLUE,
                1.0,
                "cartoon",
            ),

            structure_branch(
                "./9ewn.cif",
                A_CHAIN,
                ORANGE,
                0.55,
                "cartoon",
            ),
        ],
    ),
)


# ----------------------------------------------------------------------
# Console summary
# ----------------------------------------------------------------------

print()
print("============================================================")
print("NUMERICAL DUPLICATE CHECK")
print("============================================================")

print()
print("PAIR 1 — 1a0t:B <-> 1oh2:B")
print("----------------------------------------")
print("Stage-10 d_bri_mA           : 0")
print(f"Backbone atoms              : {len(keys_B)}")
print(f"Raw RMSD (A)                : {pair1_raw_rmsd:.9f}")
print(f"Raw maximum distance (A)    : {pair1_raw_max:.9f}")
print(f"Aligned RMSD (A)            : {pair1_aligned_rmsd:.9f}")
print(f"Aligned maximum distance (A): {pair1_aligned_max:.9f}")
print(
    f"Side-by-side separation (A) : "
    f"{pair1_side_separation:.3f}"
)

print()
print("PAIR 2 — 9eur:A <-> 9ewn:A")
print("----------------------------------------")
print("Stage-10 d_bri_mA           : 0")
print(f"Backbone atoms              : {len(keys_9)}")
print(f"Raw RMSD (A)                : {pair2_raw_rmsd:.9f}")
print(f"Raw maximum distance (A)    : {pair2_raw_max:.9f}")
print(
    f"Side-by-side separation (A) : "
    f"{pair2_side_separation:.3f}"
)

print()
print("============================================================")
print("GENERATED MOL* SCENES")
print("============================================================")
print("pair1_side_by_side.mvsj")
print("pair1_deposited.mvsj")
print("pair1_chainB_superposed.mvsj")
print("pair1_waters.mvsj")
print("pair2_side_by_side.mvsj")
print("pair2_exact_overlay.mvsj")
print()

# ======================================================================
# AUTO_FIX_MVS_TREE_V1
#
# MolViewSpec schema:
#   structure
#   ├── transform
#   └── component
#
# A component must NOT be nested below transform.
#
# This post-processing step makes every generated MVSJ conform to that
# parent-child requirement, even if an earlier builder emitted:
#
#   structure
#   └── transform
#       └── component
# ======================================================================

def _normalise_mvs_node(node):
    children = node.get("children", [])

    # First recursively normalize descendants.
    for child in children:
        _normalise_mvs_node(child)

    if node.get("kind") != "structure":
        return

    fixed_children = []

    for child in children:

        # Invalid historical form:
        #
        # structure
        #   transform
        #     component
        #
        # Convert to:
        #
        # structure
        #   transform
        #   component
        #
        if (
            child.get("kind") == "transform"
            and child.get("children")
        ):
            grandchildren = child.pop("children")

            fixed_children.append(child)
            fixed_children.extend(grandchildren)

        else:
            fixed_children.append(child)

    node["children"] = fixed_children


def _normalise_mvs_file(path):
    data = json.loads(path.read_text())

    _normalise_mvs_node(
        data["root"]
    )

    path.write_text(
        json.dumps(
            data,
            indent=2,
        )
        + "\n"
    )


_MVS_FILES_TO_FIX = [
    "pair1_side_by_side.mvsj",
    "pair1_deposited.mvsj",
    "pair1_chainB_superposed.mvsj",
    "pair1_waters.mvsj",
    "pair2_side_by_side.mvsj",
    "pair2_exact_overlay.mvsj",
]

print()
print("============================================================")
print("NORMALISING MOLVIEWSPEC TREE STRUCTURE")
print("============================================================")

for _mvs_name in _MVS_FILES_TO_FIX:
    _mvs_path = HERE / _mvs_name

    if not _mvs_path.exists():
        raise FileNotFoundError(
            f"Expected generated MVS file missing: {_mvs_path}"
        )

    _normalise_mvs_file(
        _mvs_path
    )

    print(
        "NORMALISED:",
        _mvs_name,
    )

print()
