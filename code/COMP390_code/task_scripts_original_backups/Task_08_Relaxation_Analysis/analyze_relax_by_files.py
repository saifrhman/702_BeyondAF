from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt

# ===== 1. input =====
ROOT = Path("/users/sgmwu14/scratch/openfold_inference_runs5times")

# ===== 2. output =====
OUT_DIR = Path("/users/sgmwu14/scratch/relax_analysis_results_by_files")
OUT_DIR.mkdir(exist_ok=True, parents=True)

#  match epoch / repeat / chain
UNRELAX_PAT = re.compile(
    r"epoch_(\d+)/your_model_epoch_\d+/repeat_(\d+)/predictions/([^/]+?)_model_1_unrelaxed\.pdb$"
)


def build_sample_table(root: Path) -> pd.DataFrame:
    rows = []

    unrelaxed_files = sorted(root.rglob("*_unrelaxed.pdb"))
    print(f"Found {len(unrelaxed_files)} unrelaxed pdb files")

    for unrelaxed_path in unrelaxed_files:
        path_str = str(unrelaxed_path)
        m = UNRELAX_PAT.search(path_str)
        if not m:
            continue

        epoch, repeat, chain = m.groups()

        # match relaxed filename
        relaxed_path = Path(str(unrelaxed_path).replace("_unrelaxed.pdb", "_relaxed.pdb"))

        relaxed_exists = relaxed_path.exists()

        status = "success" if relaxed_exists else "error"

        rows.append(
            {
                "epoch": int(epoch),
                "repeat": int(repeat),
                "chain": chain,
                "status": status,
                "relax_success": int(relaxed_exists),
                "relax_error": int(not relaxed_exists),
                "unrelaxed_path": str(unrelaxed_path),
                "relaxed_path": str(relaxed_path),
                "relaxed_exists": relaxed_exists,
            }
        )

    df = pd.DataFrame(rows).sort_values(["epoch", "repeat", "chain"]).reset_index(drop=True)
    return df


def summarize(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    out = (
        df.groupby(group_col)
        .agg(
            total=("status", "size"),
            success=("relax_success", "sum"),
            errors=("relax_error", "sum"),
        )
        .reset_index()
    )
    out["error_rate"] = out["errors"] / out["total"]
    out["success_rate"] = out["success"] / out["total"]
    return out


def plot_bar(summary_df: pd.DataFrame, x_col: str, y_col: str, title: str, out_path: Path):
    plt.figure(figsize=(10, 5))
    plt.bar(summary_df[x_col].astype(str), summary_df[y_col])
    plt.title(title)
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def main():
    print("Scanning root:", ROOT)
    print("Saving results to:", OUT_DIR)

    df = build_sample_table(ROOT)

    # 保存样本级总表
    sample_csv = OUT_DIR / "relax_sample_table.csv"
    df.to_csv(sample_csv, index=False)

    # 汇总表
    by_epoch = summarize(df, "epoch")
    by_chain = summarize(df, "chain")
    by_repeat = summarize(df, "repeat")

    by_epoch.to_csv(OUT_DIR / "summary_by_epoch.csv", index=False)
    by_chain.to_csv(OUT_DIR / "summary_by_chain.csv", index=False)
    by_repeat.to_csv(OUT_DIR / "summary_by_repeat.csv", index=False)

    # 图
    plot_bar(
        by_epoch.sort_values("epoch"),
        "epoch",
        "error_rate",
        "Relax Error Rate by Epoch",
        OUT_DIR / "error_rate_by_epoch.png",
    )

    plot_bar(
        by_chain.sort_values("error_rate", ascending=False),
        "chain",
        "error_rate",
        "Relax Error Rate by Protein Chain",
        OUT_DIR / "error_rate_by_chain.png",
    )

    plot_bar(
        by_repeat.sort_values("repeat"),
        "repeat",
        "error_rate",
        "Relax Error Rate by Repeat",
        OUT_DIR / "error_rate_by_repeat.png",
    )

    # 打印结果
    print()
    print("Saved sample-level table to:", sample_csv)
    print()
    print("Overall status counts:")
    print(df["status"].value_counts(dropna=False))
    print()
    print("Summary by epoch:")
    print(by_epoch.head(10))
    print()
    print("Summary by chain:")
    print(by_chain.sort_values("error_rate", ascending=False))
    print()
    print("Summary by repeat:")
    print(by_repeat.sort_values("repeat"))

    # 检查理论总数
    expected_total = 33 * 7 * 5
    print()
    print(f"Expected total samples: {expected_total}")
    print(f"Observed unrelaxed samples: {len(df)}")
    if len(df) != expected_total:
        print("WARNING: observed sample count does not match expected total!")


if __name__ == "__main__":
    main()