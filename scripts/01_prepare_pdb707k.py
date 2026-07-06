from pathlib import Path
import pandas as pd

INPUT_CSV = Path.home() / "COMP702_BeyondAF/data/checked_PDB707K_cleaned_chains_sequences_19Feb2025.csv"
OUTPUT_CSV = Path.home() / "COMP702_BeyondAF/outputs/pdb707k_model1_chain_keys.csv"

def normalise_pdb_chain(pdb_id, chain_id):
    return f"{str(pdb_id).lower()}_{str(chain_id).upper()}"

def main():
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(INPUT_CSV)

    print("Loaded rows:", len(df))
    print("Columns:", list(df.columns))

    if "model_id" in df.columns:
        df = df[df["model_id"] == 1].copy()

    df["chain_key"] = [
        normalise_pdb_chain(pdb_id, chain_id)
        for pdb_id, chain_id in zip(df["pdb_id"], df["chain_id"])
    ]

    keep_cols = [
        "chain_key",
        "pdb_id",
        "chain_id",
        "model_id",
        "chain_length",
        "seq",
    ]

    keep_cols = [c for c in keep_cols if c in df.columns]
    df = df[keep_cols].drop_duplicates("chain_key")

    df.to_csv(OUTPUT_CSV, index=False)

    print("Model-1 unique chains:", len(df))
    print("Saved to:", OUTPUT_CSV)

if __name__ == "__main__":
    main()
