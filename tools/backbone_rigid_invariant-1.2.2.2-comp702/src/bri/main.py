import argparse
import datetime
import functools
import json
import logging
import multiprocessing as mp
from pathlib import Path
import sys
import time
from typing import Optional, Union, List, Tuple, Any

import pandas as pd
import tqdm

from bri import Entry, MiniEntry, Chain, MiniChain, group_invariant_compare
from bri.filter import (
    entry_integrated_cleaning,
    minientry_integrated_cleaning,
    MINI_ENTRY_CLEAN_COL,
)


_TASK_TYPE = ["clean", "bri", "duplicate"]


# Set up basic logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger("BRI")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments for the BRI evaluation tool.

    :return: Parsed command line arguments
    """
    parser = argparse.ArgumentParser(description="BRI evaluation tool")

    parser.add_argument(
        "--task",
        type=str,
        choices=_TASK_TYPE,
        required=True,
        help="Type of task to run",
    )

    parser.add_argument(
        "--input-data",
        type=str,
        default=None,
        help="Dataset path",
    )

    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output results path",
    )

    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Maximum number of input entries (default: all)",
    )

    parser.add_argument(
        "--process",
        type=bool,
        default=False,
        action=argparse.BooleanOptionalAction,
        help="Show the real-time process of task",
    )

    parser.add_argument(
        "--extend",
        type=bool,
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Whether to use extended option for cleaning or bri computation pipline",
    )

    return parser.parse_args()


def input_process(
    input_path: Union[str, Path], start: int = 0, end: Optional[int] = None
) -> Union[List[str], List[Any], pd.DataFrame]:
    """Process input files and return data in appropriate format.

    :param input_path: Path to input file (.txt, .json, or .csv)
    :param start: Starting index for data slice
    :param end: Ending index for data slice
    :return: Processed input data
    :raises ValueError: If file type is not supported
    """
    end = end + start if end else None
    ext = Path(input_path).suffix
    ext = ext.lower()

    if ext == ".txt":
        with open(input_path, "r") as file:
            text = file.read().strip("\n")
            return text.split(",")[start:end]

    elif ext == ".json":
        with open(input_path, "r") as file:
            return json.load(file)[start:end]

    elif ext == ".csv":
        return pd.read_csv(input_path, na_values="", keep_default_na=False).iloc[
            start:end, :
        ]

    elif ext == "":
        all_paths = list(Path(input_path).glob("*"))
        return all_paths[start:end]

    else:
        raise ValueError(f"Unsupported file type: {ext}")


def pre_process_clean(input_data: Any) -> Any:
    """Pre-process input data for cleaning task.

    :param input_data: Input data to be processed
    :return: Processed data ready for cleaning
    """
    logger.info("Pre-processing data for cleaning task")
    # cleaning-specific pre-processing
    return input_data


def pre_process_bri(input_data: Any) -> Any:
    """Pre-process input data for BRI computation task.

    :param input_data: Input data to be processed
    :return: Processed data ready for BRI computation
    """
    logger.info("Pre-processing data for BRI computation task")
    # BRI-specific pre-processing
    if isinstance(input_data, pd.DataFrame):
        if set(MINI_ENTRY_CLEAN_COL).issubset(set(input_data.columns)):
            return input_data[MINI_ENTRY_CLEAN_COL].to_dict("records")
        else:
            logger.error(f"Cannot find required data for BRI computation")
            sys.exit(1)

    elif isinstance(input_data, list):
        return input_data

    # TODO pre-process for json(dict) input
    else:
        logger.error(f"Unsupported input data type for BRI computation")
        sys.exit(1)


def pre_process_duplicate(input_data: pd.DataFrame):
    """Pre-process input data for duplicate detection task.

    :param input_data: Input data to be processed
    :return: Processed data ready for duplicate detection
    """
    logger.info("Pre-processing data for duplicate detection task")

    cleaned_input = input_data.drop(columns=["seq"])

    return cleaned_input


def post_process_clean(result: Any, output_path: Optional[str] = None) -> None:
    """Post-process results after cleaning task.

    :param result: Result data from cleaning task
    :param output_path: Optional path to save results
    """
    logger.info("Post-processing cleaning results")
    # remove None type
    clean_result = [i for i in result if i is not None]
    # separate cleaned chains and noisy chains
    clean = [i[0] for i in clean_result]
    dirty = [i[1] for i in clean_result]
    chain_num = sum([i[2] for i in clean_result])
    clean_result = pd.concat(clean, ignore_index=True)
    dirty_result = pd.concat(dirty, ignore_index=True)

    output_path = output_path or "."
    clean_result.to_csv(Path(output_path) / "chains_cleaned.csv", index=False)
    dirty_result.to_csv(Path(output_path) / "chains_filtered.csv", index=False)

    # Post-process results
    logger.info(
        f"Task completed. Cleaned {len(clean_result)} chains out of {chain_num}."
    )


def post_process_bri(
    result: List[pd.DataFrame], output_path: Optional[str] = None
) -> None:
    """Post-process results after BRI computation task.

    :param result: Result data from BRI computation
    :param output_path: Optional path to save results
    """
    logger.info("Post-processing BRI computation results")
    # Filter out None results
    valid_results = [r for r in result if r is not None]
    if valid_results:
        output_path = output_path or "bri_results.csv"  # set save location
        save_loc = Path(output_path)
        if save_loc.suffix:
            bri_results = pd.concat(valid_results)
            if save_loc.suffix == ".parquet":
                bri_results.to_parquet(save_loc, index=False)
            else:
                bri_results.to_csv(save_loc, index=False)
            logger.info(f"Saved {len(valid_results)} BRIs to {save_loc}")

        else:  # Results already saved
            logger.info(f"Saved {len(valid_results)} BRIs to {save_loc}")

    else:
        logger.warning("No valid BRI results to save")


def post_process_duplicate(
    result: pd.DataFrame, output_path: Optional[str] = None
) -> None:
    """Post-process results after duplicate detection task.

    :param result: Result data from duplicate detection
    :param output_path: Optional path to save results
    """
    logger.info("Post-processing duplicate detection results")
    timestr = datetime.datetime.now().strftime("%d%m%y")

    # threshold
    result1 = result[result["distance"] < 1]
    result1.to_csv(f"L_inf_lt1_{timestr}.csv", index=False)
    result001_seq_e = result[(result["distance"] < 0.01) & (result["seq_diff"] < 1)]
    result001_seq_e.to_csv(f"L_inf_lt001_eq_seq{timestr}.csv", index=False)
    logger.info(
        f"Concise BRI distance results saved to L_inf_lt1_{timestr}.csv, L_inf_lt001_eq_seq{timestr}.csv"
    )

    # save all results by chain length only if output path specified
    if output_path:
        logger.info(f"Saving all BRI distance results to {output_path}")
        columns = "pdb_id1,model_id1,chain_id1,pdb_id2,model_id2,chain_id2,L_inf_invariant,seq_diff,seq1,seq2".split(
            ","
        )
        results_group = result.groupby("chain_length")
        for length, results in results_group:
            file_path = Path(output_path) / f"BRI_L_inf_{length}.csv"
            results.to_csv(file_path, index=False, columns=columns)


def clean_task(
    path: Union[str, Path], extend: bool = False
) -> Optional[Tuple[pd.DataFrame, pd.DataFrame]]:
    """Execute the cleaning task on a given PDB/mmCIF file.

    :param path: Path to the PDB/mmCIF file
    :param extend: Whether to use extended cleaning process
    :return: Tuple of cleaned and dirty entries, or None if cleaning fails
    """
    try:
        if extend:
            result = entry_integrated_cleaning(path)
        else:
            result = minientry_integrated_cleaning(path)
    except Exception as e:
        logger.warning(f"Failed to clean {path}. Exception: {e}")
        return None

    return result


def bri_task(
    input_id: Union[dict, str, Path], paras: argparse.Namespace
) -> Optional[pd.DataFrame]:
    """Compute Backbone Rigid Invariants for a given chain data.

    :param chain_dict: Dictionary containing chain information
    :param paras: Program arguments
    :return: Computed invariant as pandas DataFrame, or None if computation fails
    """
    # If to save separately
    save_path = None
    if paras.output is not None and Path(paras.output).is_dir():
        save_path = Path(paras.output)

    if isinstance(input_id, dict):
        try:
            chain = MiniChain(**input_id)
            invariant = chain.get_chain_invariant(angles=paras.extend)
            if save_path is not None:
                invariant.to_csv(save_path / f"{chain.__repr__()}.csv", index=False)
            return invariant
        except Exception as e:
            logger.warning(
                f"Failed to compute invariant for {input_id.get('pdb_id', 'unknown')}-{input_id.get('model_id', 'unknown')}-{input_id.get('chain_id', 'unknown')}: {e}"
            )
            return None

    elif isinstance(input_id, (str, Path)):
        try:
            entry = MiniEntry(input_id)
            invariant = entry.get_entry_invariant(angles=paras.extend)
            if save_path is not None:
                entry.save_invariants(paras.output, entry.chains, paras.extend)
            return invariant
        except Exception as e:
            logger.warning(f"Failed to compute invariant for {input_id}: {e}")
            return None

    return None


def process_invariant_with_pdb_id(chain_id):
    # get invariants of chains having same length
    invariant = MiniChain(**chain_id).invariant
    invariant["pdb_id"] = chain_id["pdb_id"]
    return invariant


def duplicate_task(dataset: pd.DataFrame, pool) -> Optional[Any]:
    """Detect duplicates in a given structure.

    :param path: Path to the PDB/mmCIF file
    :return: Duplicate detection results, or None if detection fails
    """
    compare_results = []

    # get invariants by length
    chain_groups = dataset.groupby("chain_length")
    for length, chain_group in chain_groups:
        chain_ids = chain_group.to_dict("records")

        # apply multiprocessing
        try:
            invariants_same_length = pool.map(
                process_invariant_with_pdb_id, chain_ids
            )  # get invariants in group
            invariants_same_length = list(invariants_same_length)

            fmt_bri_group = pd.concat(invariants_same_length, ignore_index=True)
            fmt_bri_group = fmt_bri_group[
                [
                    "pdb_id",
                    "model_id",
                    "chain_id",
                    "residue_id",
                    "residue_label",
                    "x(N)",
                    "y(N)",
                    "z(N)",
                    "x(A)",
                    "y(A)",
                    "z(A)",
                    "x(C)",
                    "y(C)",
                    "z(C)",
                    "chain_length",
                ]
            ]  # choose columns

            # compare invariants
            compare_result = group_invariant_compare(fmt_bri_group, seq_compare=True)
            compare_result["chain_length"] = length
            compare_results.append(compare_result)
        except Exception as e:
            logger.warning(
                f"Failed duplicate detection on {len(chain_ids)} chains having {length} residues"
            )
            logger.warning(f"{e}")

    compare_results = pd.concat(compare_results, ignore_index=True)
    return compare_results


def main():
    args = parse_args()

    # Set up file logging only after parsing args (not on --help)
    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"eval_{timestamp}.log"
    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)

    logger.info(f"Start at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("=" * 40)
    logger.info(f"Available CPUs: {mp.cpu_count()}")

    # Parse input data
    logger.info(f"Processing {args.input_data} as input")
    dataset = input_process(args.input_data, 0, args.max_samples)

    # Pre-process based on task type
    pre_process_func = globals().get(f"pre_process_{args.task}")
    if pre_process_func:
        dataset = pre_process_func(dataset)
    else:
        logger.error(f"Pre-process not defined for task '{args.task}'")

    # Assign task
    task = globals().get(f"{args.task}_task")
    if task is None:
        logger.error(f"Task function '{args.task}_task' not found")
        sys.exit(1)
    # For clean task, apply extend parameter
    if args.task == "clean":
        task = functools.partial(task, extend=args.extend)
    elif args.task == "bri":
        task = functools.partial(task, paras=args)

    # Task start
    logger.info(f"Start running task: {args.task}")
    start_time = time.time()
    if args.task != "duplicate":
        pool = mp.Pool(mp.cpu_count())
        if args.process:
            result = tqdm.tqdm(pool.imap_unordered(task, dataset), total=len(dataset))
        else:
            result = pool.map(task, dataset)

        result = list(result)
        pool.close()
    else:
        pool = mp.Pool(mp.cpu_count())
        result = task(dataset, pool)
        pool.close()

    # Post-process based on task type
    post_process_func = globals().get(f"post_process_{args.task}")
    if post_process_func:
        post_process_func(result, args.output)

    logger.info(f"Completed with total time: {(time.time() - start_time):.2f}s")
