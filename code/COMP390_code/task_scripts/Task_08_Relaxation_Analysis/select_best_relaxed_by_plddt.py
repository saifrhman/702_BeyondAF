from pathlib import Path
import os

import pandas as pd


def env_path(name: str, default: Path) -> Path:
    return Path(os.environ.get(name, str(default))).expanduser()


HOME = Path.home()
COMP702_ROOT = env_path("COMP702_ROOT", HOME / "COMP702_BeyondAF")
COMP390_ROOT = env_path("COMP390_ROOT", COMP702_ROOT / "code" / "COMP390_code")
RESULT_ROOT = env_path("COMP702_RESULT_ROOT", COMP390_ROOT / "Result")

RELAX_SAMPLE_CSV = env_path(
    "RELAX_SAMPLE_CSV",
    RESULT_ROOT / "relax_analysis_results_by_files" / "relax_sample_table.csv",
)
PLDDT_CSV = env_path(
    "UNRELAXED_PLDDT_CSV",
    RESULT_ROOT / "plddt_analysis_results" / "unrelaxed_plddt_summary.csv",
)

OUT_DIR = env_path("BEST_DISPLAY_SELECTION_DIR", RESULT_ROOT / "best_display_selection_by_plddt")
OUT_DIR.mkdir(exist_ok=True, parents=True)

OUT_BEST = OUT_DIR / "best_relaxed_file_per_epoch_chain.csv"
OUT_TOP_EPOCHS = OUT_DIR / "recommended_epochs_for_display.csv"


def main():
    if not RELAX_SAMPLE_CSV.exists():
        raise FileNotFoundError(f"Relax sample CSV not found: {RELAX_SAMPLE_CSV}")

    if not PLDDT_CSV.exists():
        raise FileNotFoundError(f"pLDDT CSV not found: {PLDDT_CSV}")

    relax_df = pd.read_csv(RELAX_SAMPLE_CSV)
    plddt_df = pd.read_csv(PLDDT_CSV)

    if "state" in plddt_df.columns:
        plddt_df = plddt_df[plddt_df["state"] == "unrelaxed"].copy()

    relax_keep = [
        "epoch", "repeat", "chain",
        "relax_success", "relax_error",
        "unrelaxed_path", "relaxed_path", "relaxed_exists",
    ]
    plddt_keep = [
        "epoch", "repeat", "chain",
        "mean_plddt", "min_plddt", "max_plddt", "num_residues", "pdb_path",
    ]

    relax_df = relax_df[[c for c in relax_keep if c in relax_df.columns]].copy()
    plddt_df = plddt_df[[c for c in plddt_keep if c in plddt_df.columns]].copy()

    required_relax = {"epoch", "repeat", "chain", "relax_success"}
    required_plddt = {"epoch", "repeat", "chain", "mean_plddt"}

    missing_relax = required_relax - set(relax_df.columns)
    missing_plddt = required_plddt - set(plddt_df.columns)

    if missing_relax:
        raise ValueError(f"Relax CSV missing required columns: {sorted(missing_relax)}")

    if missing_plddt:
        raise ValueError(f"pLDDT CSV missing required columns: {sorted(missing_plddt)}")

    merged = relax_df.merge(
        plddt_df,
        on=["epoch", "repeat", "chain"],
        how="left",
    )

    success_df = merged[(merged["relax_success"] == 1) & merged["mean_plddt"].notna()].copy()

    if success_df.empty:
        empty = pd.DataFrame()
        empty.to_csv(OUT_BEST, index=False)
        empty.to_csv(OUT_TOP_EPOCHS, index=False)
        print("No successful relaxed samples with pLDDT values were found.")
        print("Wrote empty outputs to:", OUT_DIR)
        return

    best_idx = success_df.groupby(["epoch", "chain"])["mean_plddt"].idxmax()
    best_df = success_df.loc[best_idx].copy()
    best_df = best_df.sort_values(["chain", "epoch"]).reset_index(drop=True)
    best_df.to_csv(OUT_BEST, index=False)

    rec_rows = []

    for chain, sub in best_df.groupby("chain"):
        sub = sub.sort_values("epoch").reset_index(drop=True)
        n = len(sub)

        chosen = [sub.iloc[0], sub.iloc[n // 2], sub.iloc[-1]]
        chosen_df = pd.DataFrame(chosen).drop_duplicates(subset=["epoch", "chain"])

        for _, row in chosen_df.iterrows():
            rec_rows.append(
                {
                    "chain": row["chain"],
                    "epoch": row["epoch"],
                    "repeat": row["repeat"],
                    "mean_plddt": row["mean_plddt"],
                    "relaxed_path": row.get("relaxed_path", ""),
                    "unrelaxed_path": row.get("unrelaxed_path", ""),
                }
            )

    rec_df = pd.DataFrame(rec_rows).sort_values(["chain", "epoch"]).reset_index(drop=True)
    rec_df.to_csv(OUT_TOP_EPOCHS, index=False)

    print("Saved best file table to:", OUT_BEST)
    print("Saved recommended epochs table to:", OUT_TOP_EPOCHS)
    print()
    print("Example best selections:")
    cols = [c for c in ["epoch", "chain", "repeat", "mean_plddt", "relaxed_path"] if c in best_df.columns]
    print(best_df[cols].head(20))


if __name__ == "__main__":
    main()
