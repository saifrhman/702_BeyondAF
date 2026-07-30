#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Build the model-1 PDB707K manifest for BRI geometric analysis."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--missing-pdb-ids", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, keep_default_na=False)
    df = df[df["model_id"] == 1].copy()

    with open(args.missing_pdb_ids) as f:
        missing = {line.strip().upper() for line in f if line.strip()}

    df = df[~df["pdb_id"].str.upper().isin(missing)].copy()

    # 8RX0 contains three duplicated label-chain keys in the source CSV.
    # These retained records were verified against the current RCSB mmCIF:
    # B -> start 8, length 148, auth chain B
    # F -> start 1, length 140, auth chain G
    # I -> start 1, length 75, auth chain U
    keep_8rx0 = {
        "B": (8, 148),
        "F": (1, 140),
        "I": (1, 75),
    }

    remove = pd.Series(False, index=df.index)

    for chain_id, (start, length) in keep_8rx0.items():
        key = (
            (df["pdb_id"].str.upper() == "8RX0")
            & (df["chain_id"] == chain_id)
        )
        correct = (
            (df["start_residue"] == start)
            & (df["chain_length"] == length)
        )
        remove |= key & ~correct

    df = df[~remove].copy()

    key_cols = ["pdb_id", "model_id", "chain_id"]

    if df.duplicated(key_cols).any():
        raise RuntimeError("Duplicate chain keys remain after manifest preparation.")

    cols = [
        "pdb_id",
        "model_id",
        "chain_id",
        "start_residue",
        "chain_length",
        "seq",
    ]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    df[cols].to_csv(output, index=False)

    print(f"Production chains: {len(df)}")
    print(f"Unique chain keys: {df[key_cols].drop_duplicates().shape[0]}")
    print(f"Unique PDB entries: {df['pdb_id'].str.upper().nunique()}")
    print(f"Saved: {output}")


if __name__ == "__main__":
    main()
