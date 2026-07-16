from pathlib import Path
import os
import re

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def extract_repeat(path: Path):
    m = re.search(r"repeat_(\d+)", str(path))
    return int(m.group(1)) if m else None


def extract_state(path: Path):
    name = path.name.lower()
    if "unrelaxed" in name:
        return "unrelaxed"
    if "relaxed" in name:
        return "relaxed"
    return "unknown"


def extract_plddt_from_pdb(pdb_path: Path) -> pd.DataFrame:
    rows = []
    seen = set()

    with open(pdb_path, "r", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue

            chain = line[21].strip() or "_"
            resseq = line[22:26].strip()
            icode = line[26].strip()
            resname = line[17:20].strip()

            key = (chain, resseq, icode)
            if key in seen:
                continue
            seen.add(key)

            try:
                plddt = float(line[60:66].strip())
            except ValueError:
                continue

            rows.append(
                {
                    "pdb_file": str(pdb_path),
                    "repeat": extract_repeat(pdb_path),
                    "state": extract_state(pdb_path),
                    "chain": chain,
                    "residue_index": int(resseq),
                    "insertion_code": icode,
                    "residue_name": resname,
                    "plddt": plddt,
                }
            )

    return pd.DataFrame(rows)


def confidence_band(plddt: float) -> str:
    if plddt >= 90:
        return "very_high"
    if plddt >= 70:
        return "confident"
    if plddt >= 50:
        return "low"
    return "very_low"


def main():
    root = Path(os.environ["TWOLO_PRED_ROOT"]).expanduser().resolve()
    out_dir = Path(os.environ["TWOLO_ANALYSIS_OUT"]).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb_files = sorted(root.rglob("*.pdb"))

    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found under: {root}")

    all_rows = []

    for pdb in pdb_files:
        df = extract_plddt_from_pdb(pdb)
        if not df.empty:
            all_rows.append(df)

    if not all_rows:
        raise RuntimeError("No CA pLDDT values extracted from PDB files.")

    per_residue = pd.concat(all_rows, ignore_index=True)
    per_residue["confidence_band"] = per_residue["plddt"].apply(confidence_band)

    per_file = (
        per_residue.groupby(["pdb_file", "repeat", "state"])
        .agg(
            residues=("plddt", "count"),
            mean_plddt=("plddt", "mean"),
            median_plddt=("plddt", "median"),
            min_plddt=("plddt", "min"),
            max_plddt=("plddt", "max"),
            std_plddt=("plddt", "std"),
            very_high=("confidence_band", lambda x: (x == "very_high").sum()),
            confident=("confidence_band", lambda x: (x == "confident").sum()),
            low=("confidence_band", lambda x: (x == "low").sum()),
            very_low=("confidence_band", lambda x: (x == "very_low").sum()),
        )
        .reset_index()
        .sort_values(["state", "repeat"])
    )

    per_repeat = (
        per_file.groupby(["repeat", "state"])
        .agg(
            mean_plddt=("mean_plddt", "mean"),
            min_plddt=("min_plddt", "min"),
            max_plddt=("max_plddt", "max"),
            residues=("residues", "first"),
        )
        .reset_index()
        .sort_values(["state", "repeat"])
    )

    per_residue_summary = (
        per_residue.groupby(["state", "residue_index", "residue_name"])
        .agg(
            mean_plddt=("plddt", "mean"),
            std_plddt=("plddt", "std"),
            min_plddt=("plddt", "min"),
            max_plddt=("plddt", "max"),
        )
        .reset_index()
        .sort_values(["state", "residue_index"])
    )

    per_residue.to_csv(out_dir / "2olo_per_residue_plddt.csv", index=False)
    per_file.to_csv(out_dir / "2olo_per_file_plddt_summary.csv", index=False)
    per_repeat.to_csv(out_dir / "2olo_per_repeat_summary.csv", index=False)
    per_residue_summary.to_csv(out_dir / "2olo_per_residue_summary_across_repeats.csv", index=False)

    # Plot 1: mean pLDDT per repeat
    plt.figure(figsize=(10, 4))
    for state, sub in per_repeat.groupby("state"):
        plt.plot(sub["repeat"], sub["mean_plddt"], marker="o", label=state)
    plt.xlabel("Repeat")
    plt.ylabel("Mean pLDDT")
    plt.title("2OLO pretrained model: mean pLDDT across repeats")
    plt.ylim(0, 100)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "2olo_mean_plddt_by_repeat.png", dpi=200)
    plt.close()

    # Plot 2: residue-level profile, unrelaxed preferred
    profile_state = "unrelaxed" if "unrelaxed" in set(per_residue_summary["state"]) else per_residue_summary["state"].iloc[0]
    profile = per_residue_summary[per_residue_summary["state"] == profile_state]

    plt.figure(figsize=(10, 4))
    plt.plot(profile["residue_index"], profile["mean_plddt"])
    plt.fill_between(
        profile["residue_index"],
        profile["mean_plddt"] - profile["std_plddt"].fillna(0),
        profile["mean_plddt"] + profile["std_plddt"].fillna(0),
        alpha=0.25,
    )
    plt.xlabel("Residue index")
    plt.ylabel("Mean pLDDT")
    plt.title(f"2OLO pretrained model: residue pLDDT profile ({profile_state})")
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(out_dir / "2olo_residue_plddt_profile.png", dpi=200)
    plt.close()

    # Plot 3: distribution
    plt.figure(figsize=(7, 4))
    plt.hist(per_residue["plddt"], bins=20)
    plt.xlabel("pLDDT")
    plt.ylabel("Residue count")
    plt.title("2OLO pretrained model: pLDDT distribution")
    plt.tight_layout()
    plt.savefig(out_dir / "2olo_plddt_distribution.png", dpi=200)
    plt.close()

    print("Input root:", root)
    print("Output directory:", out_dir)
    print()
    print("Number of PDB files analysed:", len(pdb_files))
    print("Number of residue records:", len(per_residue))
    print()
    print("Per-state summary:")
    print(
        per_file.groupby("state")
        .agg(
            files=("pdb_file", "count"),
            mean_plddt_mean=("mean_plddt", "mean"),
            mean_plddt_std=("mean_plddt", "std"),
            min_plddt=("min_plddt", "min"),
            max_plddt=("max_plddt", "max"),
        )
        .reset_index()
        .to_string(index=False)
    )
    print()
    print("Lowest-confidence residues by mean pLDDT:")
    print(profile.sort_values("mean_plddt").head(15).to_string(index=False))
    print()
    print("Saved outputs:")
    for p in sorted(out_dir.iterdir()):
        print(p)


if __name__ == "__main__":
    main()
