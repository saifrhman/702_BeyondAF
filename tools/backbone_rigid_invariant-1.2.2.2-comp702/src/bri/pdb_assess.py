import multiprocessing as mp
import logging
from pathlib import Path
from functools import partial

import click
import matplotlib.pyplot as plt
import natsort
import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree

# import plotly.graph_objects as go
# from plotly.subplots import make_subplots

# TODO: tausion angle -> [-90,270]
from bri.invariant_compare import (
    INVARIANT_VALUE_COLS,
    extract_chain_info,
    generate_index_table,
    convert_index_table,
    coordinate_value_reshape,
)
from bri.invariant import get_invariant_BTP, get_invariant
from bri.base.math_base import FloatArray
from bri.filter import residue_continuity_check
from bri.pdbx2df import MiniChain
from bri.pdb707k_extract import extract_vectors

# from bri.filter import integrated_chainwise_filter, residue_continuity_check

BTP_ATTR = [
    "|C-N-A|",
    "|N-A-C|",
    "|A-C-N|",
    "TP(NA)_x",
    "TP(NA)_y",
    "TP(AC)_x",
    "TP(AC)_y",
    "TP(CN)_x",
    "TP(CN)_y",
]

LAI = [
    "length(N)",
    "length(A)",
    "length(C)",
    "angle(N)",
    "angle(A)",
    "angle(C)",
    "tau(NA)",
    "tau(AC)",
    "tau(CN)",
]


def get_RMS(x: FloatArray, y: FloatArray, n_null: int = 0) -> float:
    """Computes RMSD between 2 point sets."""
    N = np.size(x) - n_null
    return np.sqrt(np.sum((x - y) ** 2) / N)


def group_invariant_compare(
    x: pd.DataFrame,
    y: pd.DataFrame | None = None,
    seq_compare: bool = False,
    BRI_flag: bool = False,
    metric: str = "RMS",
) -> pd.DataFrame | None:
    """Compare backbone rigid invariants using either RMS or Chebyshev (BRI) metrics."""

    same_col = ["pdb_id", "model_id", "chain_id"]
    seq_compare_col = ["pdb_id", "model_id", "chain_id", "seq"]
    reduce_matrix = False

    if y is None:
        y = x
        reduce_matrix = True

    chains_x = extract_chain_info(x)
    chains_y = extract_chain_info(y)

    if len(chains_x) < 1 or len(chains_y) < 1:
        return None

    chain_lengths = list(chains_x["chain_length"].unique())
    if len(chain_lengths) > 1:
        logging.warning("Chain lengths are not equal.")
        return None
    chain_length = chain_lengths[0]

    attr_list = INVARIANT_VALUE_COLS if BRI_flag else BTP_ATTR
    base_value = coordinate_value_reshape(x, chain_length, attr_list)
    search_value = coordinate_value_reshape(y, chain_length, attr_list)

    if metric == "RMS":
        RMS = partial(get_RMS, n_null=6)
        tree = BallTree(base_value, metric=RMS)
    else:
        tree = BallTree(base_value, metric="chebyshev")

    nn_res = tree.query(search_value, k=len(base_value), return_distance=True)
    idx_table = generate_index_table(nn_res, reduce_matrix)
    distance_res = convert_index_table(idx_table, chains_x, chains_x, same_col)

    if not seq_compare:
        return distance_res

    # Sequence distance comparison
    base_seq_value = coordinate_value_reshape(x, chain_length, ["residue_label_ascii"])
    targ_seq_value = coordinate_value_reshape(y, chain_length, ["residue_label_ascii"])
    seq_tree = BallTree(base_seq_value, metric="hamming")
    nn_res_seq = seq_tree.query(targ_seq_value, k=len(chains_x), return_distance=True)

    idx_table_seq = generate_index_table(nn_res_seq, reduce_matrix)
    result_seq = convert_index_table(idx_table_seq, chains_x, chains_x, seq_compare_col)

    result_seq = result_seq.rename(columns={"distance": "seq_diff"})
    result_seq["seq_diff"] = (
        (result_seq["seq_diff"] * chain_length).round().astype("int")
    )

    merge_keys = [f"{col}1" for col in same_col] + [f"{col}2" for col in same_col]
    return distance_res.merge(result_seq, on=merge_keys, how="inner")


def get_data(folder: Path) -> pd.DataFrame:
    """Load and concatenate all CSVs in a target directory."""
    files = folder.glob("*.csv")
    data = []
    for f in files:
        df = pd.read_csv(f)
        df["pdb_id"] = "AF2".join(f.stem.split("alphafold2_ptm_model"))[2:-4]
        df["chain_length"] = len(df)
        df = df.fillna(0)
        data.append(df)

    if not data:
        raise ValueError(f"No CSV files found in {folder}")
    return pd.concat(data)


def compute_distance_matrix(input: Path, output: Path):
    """Calculate and save RMS and Linf distance matrices and a single comprehensive CSV table."""

    process_distance_matrices(input, output)
    process_distance_tables(input, output)


def distance_table_to_matrix(
    distance_table: pd.DataFrame,
    id_col1: str = "prediction_id1",
    id_col2: str = "prediction_id2",
    distance_col: str = "distance",
) -> pd.DataFrame:
    """Convert a distance table to a symmetric distance matrix."""
    all_ids = natsort.natsorted(
        set(distance_table[id_col1].unique()) | set(distance_table[id_col2].unique())
    )
    matrix = pd.DataFrame(index=all_ids, columns=all_ids, dtype=float)

    for _, row in distance_table.iterrows():
        id1, id2, dist = row[id_col1], row[id_col2], row[distance_col]
        matrix.loc[id1, id2] = dist
        matrix.loc[id2, id1] = dist

    np.fill_diagonal(matrix.values, 0.0)
    return matrix


# --- Workflow Execution Functions ---
def process_distance_matrices(input_dir: Path, output_dir: Path):
    """Generates RMS and Linf distance matrices."""
    data = get_data(input_dir)

    configurations = [
        ("RMS", True, "distance_matrix_BRI_RMS.csv"),
        ("Chebyshev", True, "distance_matrix_BRI_Linf.csv"),
        ("RMS", False, "distance_matrix_BTI_RMS.csv"),
        ("Chebyshev", False, "distance_matrix_BTI_Linf.csv"),
    ]

    for metric, is_bri, filename in configurations:
        dist = group_invariant_compare(data, BRI_flag=is_bri, metric=metric)
        if dist is not None:
            matrix = distance_table_to_matrix(dist, "pdb_id1", "pdb_id2")
            matrix.to_csv(output_dir / filename)
            logging.info(f"Saved matrix: {filename}")


def process_distance_tables(input_dir: Path, output_dir: Path):
    """Compiles distance tables into a single comprehensive CSV."""
    data = get_data(input_dir)

    dist_bri_linf = group_invariant_compare(
        data, BRI_flag=True, metric="Chebyshev"
    ).rename(columns={"distance": "BRI_Linf_dist, Angstrom"})
    dist_bri_rms = group_invariant_compare(data, BRI_flag=True, metric="RMS").rename(
        columns={"distance": "BRI_RMS_dist, Angstrom"}
    )
    dist_lai_linf = group_invariant_compare(data, metric="Chebyshev").rename(
        columns={"distance": "BTI_Linf_dist, Angstrom"}
    )
    dist_lai_rms = group_invariant_compare(data, metric="RMS").rename(
        columns={"distance": "BTI_RMS_dist, Angstrom"}
    )

    bri_merged = pd.merge(dist_bri_linf, dist_bri_rms, on=["pdb_id1", "pdb_id2"])
    lai_merged = pd.merge(dist_lai_linf, dist_lai_rms, on=["pdb_id1", "pdb_id2"])
    final_table = pd.merge(bri_merged, lai_merged, on=["pdb_id1", "pdb_id2"])

    final_table = final_table[
        [
            "BRI_Linf_dist, Angstrom",
            "BRI_RMS_dist, Angstrom",
            "BTI_Linf_dist, Angstrom",
            "BTI_RMS_dist, Angstrom",
            "pdb_id1",
            "pdb_id2",
        ]
    ]
    final_table = final_table.rename(
        columns={"pdb_id1": "prediction_id1", "pdb_id2": "prediction_id2"}
    )

    output_path = output_dir / "distance_table.csv"
    _ = final_table.to_csv(output_path, index=False)
    logging.info(f"Saved distance table: {output_path}")


def refine_invariant_data(target_dir: Path):
    """Refines invariants into a clean target directory."""
    source_dir = target_dir / "compute_invariants"
    out_dir = target_dir / "invariants"
    out_dir.mkdir(exist_ok=True)

    for f in source_dir.glob("*.csv"):
        df = pd.read_csv(f).rename(columns={"pdb_id": "prediction_id"}).fillna(0)
        file_name = (
            f"AF2{''.join(f.stem.split('alphafold2_ptm_model')[1:])}_BRI_LAI_BTI.csv"
        )
        df.to_csv(out_dir / file_name, index=False)
    logging.info(f"Refined files saved to {out_dir}")


# --- Visualization ---
def scatter_projection(input: Path, output: Path):

    data = get_data(input)

    bri_stat = (
        data.groupby(["pdb_id"])[INVARIANT_VALUE_COLS]
        .agg(["mean", "std"])
        .reset_index()
    )
    bri_stat.columns = ["_".join(i) for i in bri_stat.columns]
    _ = bri_stat.to_csv(output / f"{input.name}_BRI_proj.csv", index=False)
    btp_stat = data.groupby(["pdb_id"])[BTP_ATTR].agg(["mean"]).reset_index()
    btp_stat.columns = ["_".join(i) for i in btp_stat.columns]
    _ = btp_stat.to_csv(output / f"{input.name}_BTI_proj.csv", index=False)

    # creating subplots
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    idx = 0
    for i in range(3):
        for j in range(3):
            ax = axes[i, j]
            ax.scatter(bri_stat.iloc[:, idx + 1], bri_stat.iloc[:, idx + 2], alpha=0.6)
            title = bri_stat.columns[idx + 1][:-5]
            ax.set_title(title)
            ax.set_xlabel(bri_stat.columns[idx + 1])
            ax.set_ylabel(bri_stat.columns[idx + 2])
            idx += 2

    fig.tight_layout()

    save_name = f"{input.name}_BRI_proj.png"
    save_path = output / save_name
    plt.savefig(save_path, dpi=300)

    # creating subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    idx = 3
    for i in range(3):
        ax = axes[i]
        ax.scatter(btp_stat.iloc[:, idx + 1], btp_stat.iloc[:, idx + 2], alpha=0.6)
        title = btp_stat.columns[idx + 1][:-7]
        ax.set_title(title)
        ax.set_xlabel(btp_stat.columns[idx + 1])
        ax.set_ylabel(btp_stat.columns[idx + 2])
        idx += 2

    fig.tight_layout()

    save_name = f"{input.name}_BTI_proj.png"
    save_path = output / save_name
    plt.savefig(save_path, dpi=300)


def plot_comparison_curve(
    pdb_inv_root: Path, columns: list[str], ref_path: str = "", off_set: int = 0
):
    """Plot features from BRI files against ColabFold data."""

    pdb_data = {col: [] for col in columns}

    # Load ploting data
    tick_r = 0
    for file in pdb_inv_root.glob("*.csv"):
        try:
            df = pd.read_csv(file)
            tick_r = len(df)
            for col in columns:
                if col in df.columns:
                    values: FloatArray = pd.to_numeric(df[col], errors="coerce")
                    pdb_data[col].append(values)
        except Exception as e:
            logging.warning(f"Failed to read {file}: {e}")

    fig, axes = plt.subplots(len(columns), 1, figsize=(15, 8), sharex=True)

    # Plot loaded data
    for i, col in enumerate(columns):
        ax = axes[i] if len(columns) > 1 else axes
        for arr in pdb_data[col]:
            ax.plot(range(1, 1 + len(arr)), arr, color="lightgray", alpha=0.4)

        ax.set_ylabel(col)
        _apply_y_limits(ax, col)
        ax.set_xticks(list(range(1, tick_r, 50)) + [tick_r])
        ax.grid(alpha=0.3)

    if ref_path:
        ref_info = Path(ref_path)
        if ref_info.is_file():
            logging.warning("File reference not supported.")

        try:
            if len(ref_path.split("-")) > 1:
                p_id, m_id, c_id = ref_path.split("-")
                chain = MiniChain(p_id, int(m_id), c_id, 1, 750)
            else:
                chain = MiniChain(ref_path, 1, "A", 1, 750)
        except Exception:
            raise ValueError(
                f"Reference {ref_path} cannot be found. Please check the format or the accessibility of PDB."
            )

        try:
            missing_residue = residue_continuity_check(chain.get_feature())
            if missing_residue is not None:
                missing_id_list = [
                    int(i) for i in missing_residue.at[0, "missed_residues"].split(",")
                ]
                seq_start = find_sequence_boundaries(
                    missing_id_list, chain._start_residue
                )
                pdb_invariant = plot_broken_chain(chain, seq_start, columns)

            else:
                pdb_invariant = [get_invariant(chain.get_feature(), ext=True)]
                if not set(columns).issubset(set(pdb_invariant[0].columns)):
                    pdb_invariant = [get_invariant_BTP(pdb_invariant[0])]
            # pdb_invariant.fillna(0, inplace=True)
        except Exception as e:
            logging.error(f"Error loading invariants from {ref_path} due to: {e}")
            raise RuntimeError(f"Error loading invariants from {ref_path} due to: {e}")

        for i, col in enumerate(columns):
            ax = axes[i] if len(columns) > 1 else axes
            for invariant in pdb_invariant:
                pdb_res_id_col = invariant["residue_id"]
                if off_set:
                    pdb_res_id_col = pdb_res_id_col - off_set
                pdb_res_id = pdb_res_id_col.to_list()
                ax.plot(
                    pdb_res_id,
                    invariant[col],
                    color="darkred",
                    linewidth=1.0,
                    marker="o",
                    markersize=1.5,
                    label="PDB",
                )

    plt.xlabel("Residue Index")
    fig.tight_layout()
    return fig


def plot_broken_chain(chain: MiniChain, seq_start: list[int], cols: list[str]):
    path = chain.path
    m_id = chain._model_id
    c_id = chain._chain_id

    length = np.array(seq_start[1:]) - np.array(seq_start[:-1])
    segments = list(zip(seq_start[:-1], length))[::2]
    if set(cols).issubset(set(BTP_ATTR)):
        std_inv = [
            MiniChain(path, m_id, c_id, *seg).get_chain_invariant_BTP()
            for seg in segments
        ]
    else:
        std_inv = [
            MiniChain(path, m_id, c_id, *seg).get_chain_invariant(True)
            for seg in segments
        ]

    return std_inv


def find_sequence_boundaries(arr: list[int], domain_start: int = 1):
    """Finds the starting numbers of continuous sequences (both present and missing)."""
    if not arr:
        return [domain_start]

    boundaries = []

    # Beginning of the domain
    if arr[0] > domain_start:
        boundaries.append(domain_start)
    boundaries.append(arr[0])

    # 2. Iterate through the array to find gaps (state transitions)
    for i in range(1, len(arr)):
        difference = arr[i] - arr[i - 1]

        if difference > 1:
            # A gap is detected
            missing_start = arr[i - 1] + 1
            present_start = arr[i]

            boundaries.append(missing_start)
            boundaries.append(present_start)

    return boundaries


def _apply_y_limits(ax, col: str):
    """Helper to consistently apply plotting limits based on data type."""
    if col in INVARIANT_VALUE_COLS:
        ax.set_ylim(-2.3, 2.3)
        ax.set_yticks([-2, -1, 0, 1, 2])
    elif col.startswith("tau"):
        ax.set_ylim(-200, 200)
    elif col.startswith("angle"):
        ax.set_ylim(0, 180)
        ax.set_yticks([0, 45, 90, 135, 180])
    elif col.startswith("length"):
        ax.set_ylim(0, 2.5)
        ax.set_yticks([0, 1, 2])
    elif col.startswith("|"):
        ax.set_ylim(0, 4.5)
        ax.set_yticks([0, 1, 2, 3, 4])
    elif col in BTP_ATTR:
        ax.set_ylim(-2.3, 2.3)
        ax.set_yticks([-2, 0, 2])


def compute_dir_invariants(input: Path, output: Path, n_process: int):
    "Compute and save invariants from given .pdb files."

    pdb_files = list(input.glob("*.pdb"))
    res_folder = Path(output / "compute_invariants")
    if not res_folder.exists():
        res_folder.mkdir()

    pool = mp.Pool(n_process)
    func = partial(compute_and_save, res_folder=res_folder)
    result = pool.map(func, pdb_files)
    pool.close()
    result = sum(list(result))


def compute_and_save(file: Path, res_folder: Path):
    chain = MiniChain.from_pdb(file)
    bri = chain.get_chain_invariant(True)
    bri.drop(columns=["chain_length"], inplace=True)
    btp = chain.get_chain_invariant_BTP()
    res = bri.merge(btp)
    res.to_csv(res_folder / f"{file.stem}_inv.csv", index=False)
    return 0


@click.group()
def cli():
    """Protein Backbone Invariant Analysis CLI Tool.

    A suite of tools for processing, analyzing, and visualizing rigid
    backbone invariants and torsion parameters of protein structures.
    """
    pass


@cli.command("inv")
@click.argument(
    "input_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.argument(
    "output_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option(
    "-n",
    "--n_process",
    default=int(mp.cpu_count() // 2),
    show_default=True,
    type=int,
    help="Number of CPU cores to utilize for parallel processing.",
)
def compute_invariants(input_dir: Path, output_dir: Path, n_process: int):
    """Calculates invariants from raw PDB files and saves them as CSVs, including Backbone Rigid Invariants (BRI), Length Angle Invariants (LAI) and Backbone Torsion Invariants (BTI).

    \b
    INPUT_DIR: Directory containing your source .pdb files.
    OUTPUT_DIR: Directory where the resulting invariant CSVs will be saved."""

    compute_dir_invariants(input_dir, output_dir, n_process)


@cli.command("compare")
@click.argument(
    "input_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.argument(
    "output_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
def run_matrix(input_dir: Path, output_dir: Path):
    """Calculate and save RMS and Chebyshev (L-infinity) distance matrices.

    Measures distances between all processed protein chain invariants in the
    input directory and compiles them into symmetric distance matrices and a
    single comprehensive CSV table.

    \b
    INPUT_DIR: Directory containing computed invariant CSVs (e.g., from 'inv').
    OUTPUT_DIR: Directory where the distance matrices and table will be saved.
    """
    compute_distance_matrix(input_dir, output_dir)


@cli.command("proj")
@click.argument(
    "input_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.argument(
    "output_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
def run_projection(input_dir: Path, output_dir: Path):
    """Calculate statistical summaries and generate scatter projections.

    Calculates the mean and standard deviation for the invariants and generates
    2D scatter plots to visualize these distributions across the dataset.

    \b
    INPUT_DIR: Directory containing computed invariant CSVs.
    OUTPUT_DIR: Directory where the statistical CSVs and plots will be saved.
    """

    scatter_projection(input_dir, output_dir)


@cli.command("plot")
@click.argument(
    "input_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.argument(
    "output_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option(
    "-s",
    "--structure",
    type=str,
    help="Experimental structure for reference to compare against (e.g., 1HHO-1-A will specify chain `A` of model `1` in entry `1HHO`, or simply 1HHO for short).",
)
@click.option(
    "-ofs",
    "--offset",
    type=int,
    help="Use the existing gap between invariants and the reference to help alignment.",
)
def run_plot(input_dir: Path, structure: str, output_dir: Path, offset: int):
    """Generate comparison plots within the dataset.

    Creates overlay plots comparing the invariant sequences
    of computed invariants, and optionally against a known experimental reference structure.

    \b
    INPUT_DIR: Directory containing computed invariant CSVs.
    OUTPUT_DIR: Directory where the comparison PNG plots will be saved.
    """

    inv_dir = input_dir
    if not inv_dir.exists():
        logging.error(f"Invariants directory not found at {inv_dir}.")
        return

    fig = plot_comparison_curve(inv_dir, INVARIANT_VALUE_COLS, structure, offset)
    save_path = output_dir / "BRI_comparison.png"
    fig.savefig(save_path, dpi=300)
    logging.info(f"Comparison figure saved: {save_path}")

    fig = plot_comparison_curve(inv_dir, BTP_ATTR, structure, offset)
    save_path = output_dir / "BTI_comparison.png"
    fig.savefig(save_path, dpi=300)
    logging.info(f"Comparison figure saved: {save_path}")

    fig = plot_comparison_curve(inv_dir, LAI, structure, offset)
    save_path = output_dir / "LAI_comparison.png"
    fig.savefig(save_path, dpi=300)
    logging.info(f"Comparison figure saved: {save_path}")


@cli.command("pipe")
@click.argument(
    "input_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.argument(
    "output_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option(
    "-n",
    "--n_process",
    default=int(mp.cpu_count() // 2),
    show_default=True,
    type=int,
    help="Number of CPU processors used to accelerate computation",
)
def pipeline(input_dir: Path, output_dir: Path, n_process: int):
    """Execute the full pipeline: invariants, matrices, and projections.

    Sequentially executes the 'inv', 'compare', and 'proj' commands.
    Ideal for processing a raw batch of PDBs through to final statistical
    projections in a single execution.

    \b
    INPUT_DIR: Directory containing raw source .pdb files.
    OUTPUT_DIR: Directory where all resulting CSVs, matrices, and plots will be saved.
    """

    inv_dir = Path(output_dir / "compute_invariants")
    compute_dir_invariants(input_dir, output_dir, n_process)
    compute_distance_matrix(inv_dir, output_dir)
    scatter_projection(inv_dir, output_dir)


@cli.command("pdb707k-extract")
@click.argument(
    "manifest",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.argument(
    "mmcif_root",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.argument(
    "output_dir",
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--task-index",
    required=True,
    type=click.IntRange(min=0),
    help="Zero-based batch index.",
)
@click.option(
    "--batch-size",
    required=True,
    type=click.IntRange(min=1),
    help="Number of unique PDB entries processed in this batch.",
)
def run_pdb707k_extract(
    manifest: Path,
    mmcif_root: Path,
    output_dir: Path,
    task_index: int,
    batch_size: int,
):
    """Extract compact BRI vectors for one PDB707K batch.

    Uses the same MiniChain invariant representation and coordinate reshaping
    mathematics as the installed BRI 1.2.2.2 package.
    """
    extract_vectors(
        manifest_path=manifest,
        mmcif_root_path=mmcif_root,
        output_dir_path=output_dir,
        task_index=task_index,
        batch_size=batch_size,
    )


if __name__ == "__main__":
    cli()
