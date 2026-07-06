from pathlib import Path
import os
import re
import pandas as pd

# ===== 1. Centralised input paths =====
# These paths can be overridden using environment variables.
# By default, they assume the standard COMP702 Barkla2 layout:
# ~/COMP702_BeyondAF/code/COMP390_code/Result/...

COMP702_ROOT = Path(os.environ.get("COMP702_ROOT", Path.home() / "COMP702_BeyondAF"))
COMP390_ROOT = Path(os.environ.get("COMP390_ROOT", COMP702_ROOT / "code" / "COMP390_code"))
RESULT_ROOT = Path(os.environ.get("COMP702_RESULT_ROOT", COMP390_ROOT / "Result"))

# Root folder containing repeated OpenFold inference outputs.
ROOT = Path(os.environ.get("OPENFOLD_RUNS5_ROOT", RESULT_ROOT / "openfold_inference_runs5times"))

# ===== 2. Centralised output path =====
# Output folder for the extracted pLDDT summary CSV.
OUT_DIR = Path(os.environ.get("PLDDT_OUT_DIR", RESULT_ROOT / "plddt_analysis_results"))
OUT_DIR.mkdir(exist_ok=True, parents=True)

OUT_CSV = Path(os.environ.get("UNRELAXED_PLDDT_CSV", OUT_DIR / "unrelaxed_plddt_summary.csv"))

# Match only unrelaxed prediction PDB paths
PATH_PAT = re.compile(
    r"epoch_(\d+)/your_model_epoch_\d+/repeat_(\d+)/predictions/([^/]+?)_model_1_unrelaxed\.pdb$"
)

def extract_plddt_from_pdb(pdb_path: Path):
    """
    Extract whole-protein pLDDT from a PDB file.

    OpenFold stores pLDDT values in the B-factor column.
    This function uses the CA atom B-factor for each residue,
    then returns summary statistics across the chain.
    """
    vals = []

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                if atom_name == "CA":
                    try:
                        bfactor = float(line[60:66].strip())
                        vals.append(bfactor)
                    except ValueError:
                        continue

    if not vals:
        return None

    return {
        "mean_plddt": sum(vals) / len(vals),
        "min_plddt": min(vals),
        "max_plddt": max(vals),
        "num_residues": len(vals),
    }

def main():
    rows = []

    # Scan only unrelaxed PDB files
    pdb_files = sorted(ROOT.rglob("*_unrelaxed.pdb"))

    print(f"Found {len(pdb_files)} unrelaxed pdb files")

    for pdb_path in pdb_files:
        rel_str = str(pdb_path)
        m = PATH_PAT.search(rel_str)
        if not m:
            continue

        epoch, repeat, chain = m.groups()
        plddt_stats = extract_plddt_from_pdb(pdb_path)

        if plddt_stats is None:
            print(f"WARNING: no CA atoms found in {pdb_path}")
            continue

        row = {
            "epoch": int(epoch),
            "repeat": int(repeat),
            "chain": chain,
            "state": "unrelaxed",
            "pdb_path": str(pdb_path),
            **plddt_stats,
        }
        rows.append(row)

    df = pd.DataFrame(rows).sort_values(
        ["epoch", "repeat", "chain"]
    ).reset_index(drop=True)

    df.to_csv(OUT_CSV, index=False)

    print(f"Saved summary to: {OUT_CSV}")
    print()
    print("Rows:", len(df))
    print()
    print("Example rows:")
    print(df.head())

if __name__ == "__main__":
    main()