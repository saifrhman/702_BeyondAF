from pathlib import Path
import os

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def env_path(name: str, default: Path) -> Path:
    return Path(os.environ.get(name, str(default))).expanduser()


HOME = Path.home()
COMP702_ROOT = env_path("COMP702_ROOT", HOME / "COMP702_BeyondAF")
COMP390_ROOT = env_path("COMP390_ROOT", COMP702_ROOT / "code" / "COMP390_code")
RESULT_ROOT = env_path("COMP702_RESULT_ROOT", COMP390_ROOT / "Result")

INPUT_CSV = env_path(
    "RELAX_SAMPLE_CSV",
    RESULT_ROOT / "relax_analysis_results_by_files" / "relax_sample_table.csv",
)

OUT_DIR = env_path(
    "DISPLAY_SELECTION_OUT_DIR",
    RESULT_ROOT / "display_selection_analysis",
)
OUT_DIR.mkdir(exist_ok=True, parents=True)


def main():
    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Input CSV not found: {INPUT_CSV}")

    df = pd.read_csv(INPUT_CSV)

    required_cols = {"epoch", "repeat", "chain", "relax_success", "relax_error"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    summary = (
        df.groupby(["epoch", "chain"])
        .agg(
            total=("relax_success", "size"),
            success=("relax_success", "sum"),
            error=("relax_error", "sum"),
        )
        .reset_index()
    )

    summary["success_rate"] = summary["success"] / summary["total"]
    summary["all5_success"] = (summary["success"] == 5).astype(int)
    summary["atleast4_success"] = (summary["success"] >= 4).astype(int)
    summary["atleast3_success"] = (summary["success"] >= 3).astype(int)

    summary = summary.sort_values(["chain", "epoch"]).reset_index(drop=True)
    summary.to_csv(OUT_DIR / "epoch_chain_relax_counts.csv", index=False)

    all5_df = summary[summary["success"] == 5].copy().sort_values(["chain", "epoch"])
    all5_df.to_csv(OUT_DIR / "epoch_chain_all5_success.csv", index=False)

    atleast4_df = summary[summary["success"] >= 4].copy().sort_values(["chain", "epoch"])
    atleast4_df.to_csv(OUT_DIR / "epoch_chain_atleast4_success.csv", index=False)

    atleast3_df = summary[summary["success"] >= 3].copy().sort_values(["chain", "epoch"])
    atleast3_df.to_csv(OUT_DIR / "epoch_chain_atleast3_success.csv", index=False)

    matrix_success = summary.pivot(index="chain", columns="epoch", values="success")
    matrix_success.to_csv(OUT_DIR / "chain_epoch_success_matrix.csv")

    matrix_all5 = summary.pivot(index="chain", columns="epoch", values="all5_success")
    matrix_all5.to_csv(OUT_DIR / "chain_epoch_all5_matrix.csv")

    by_chain = (
        summary.groupby("chain")
        .agg(
            total_epochs=("epoch", "nunique"),
            total_success_runs=("success", "sum"),
            total_error_runs=("error", "sum"),
            num_all5_epochs=("all5_success", "sum"),
            num_atleast4_epochs=("atleast4_success", "sum"),
            num_atleast3_epochs=("atleast3_success", "sum"),
            mean_success_per_epoch=("success", "mean"),
        )
        .reset_index()
        .sort_values(["num_all5_epochs", "mean_success_per_epoch"], ascending=False)
    )
    by_chain.to_csv(OUT_DIR / "chain_display_priority_summary.csv", index=False)

    top_good = by_chain.head(10)
    top_good.to_csv(OUT_DIR / "recommended_chains_for_display.csv", index=False)

    if all5_df.empty:
        first_all5 = pd.DataFrame(
            columns=["chain", "first_epoch_all5", "last_epoch_all5", "num_all5_epochs"]
        )
    else:
        first_all5 = (
            all5_df.groupby("chain")
            .agg(
                first_epoch_all5=("epoch", "min"),
                last_epoch_all5=("epoch", "max"),
                num_all5_epochs=("epoch", "count"),
            )
            .reset_index()
            .sort_values("first_epoch_all5")
        )

    first_all5.to_csv(OUT_DIR / "chain_first_epoch_with_all5_success.csv", index=False)

    if not by_chain.empty:
        plt.figure(figsize=(9, 5))
        plt.bar(by_chain["chain"], by_chain["total_success_runs"])
        plt.xticks(rotation=45)
        plt.ylabel("Total successful relaxed runs")
        plt.title("Total Relax Success Count by Chain")
        plt.tight_layout()
        plt.savefig(OUT_DIR / "success_count_by_chain.png", dpi=200)
        plt.close()

        plt.figure(figsize=(9, 5))
        plt.bar(by_chain["chain"], by_chain["num_all5_epochs"])
        plt.xticks(rotation=45)
        plt.ylabel("Number of epochs with 5/5 relaxed results")
        plt.title("All-5-Success Epoch Count by Chain")
        plt.tight_layout()
        plt.savefig(OUT_DIR / "all5_count_by_chain.png", dpi=200)
        plt.close()

    print("Saved outputs to:", OUT_DIR)
    print()
    print("Top chains for display:")
    print(top_good)
    print()
    print("First epoch with 5/5 success per chain:")
    print(first_all5)
    print()
    print("Example of all 5/5 epoch-chain combinations:")
    print(all5_df.head(20))


if __name__ == "__main__":
    main()
