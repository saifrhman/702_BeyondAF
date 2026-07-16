from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# 1. 输入文件
# =========================
INPUT_CSV = Path("/users/sgmwu14/scratch/relax_analysis_results_by_files/relax_sample_table.csv")

# =========================
# 2. 输出目录
# =========================
OUT_DIR = Path("/users/sgmwu14/scratch/display_selection_analysis")
OUT_DIR.mkdir(exist_ok=True, parents=True)

def main():
    df = pd.read_csv(INPUT_CSV)

    # 保险检查
    required_cols = {"epoch", "repeat", "chain", "relax_success", "relax_error"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # =========================
    # 3. 按 epoch + chain 汇总
    # =========================
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

    # =========================
    # 4. 直接给你可选展示对象
    # =========================
    all5_df = summary[summary["success"] == 5].copy().sort_values(["chain", "epoch"])
    all5_df.to_csv(OUT_DIR / "epoch_chain_all5_success.csv", index=False)

    atleast4_df = summary[summary["success"] >= 4].copy().sort_values(["chain", "epoch"])
    atleast4_df.to_csv(OUT_DIR / "epoch_chain_atleast4_success.csv", index=False)

    atleast3_df = summary[summary["success"] >= 3].copy().sort_values(["chain", "epoch"])
    atleast3_df.to_csv(OUT_DIR / "epoch_chain_atleast3_success.csv", index=False)

    # =========================
    # 5. 做矩阵，方便看趋势
    # =========================
    matrix_success = summary.pivot(index="chain", columns="epoch", values="success")
    matrix_success.to_csv(OUT_DIR / "chain_epoch_success_matrix.csv")

    matrix_all5 = summary.pivot(index="chain", columns="epoch", values="all5_success")
    matrix_all5.to_csv(OUT_DIR / "chain_epoch_all5_matrix.csv")

    # =========================
    # 6. 统计每条 chain 有多少个 epoch 达到 5/5
    # =========================
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

    # =========================
    # 7. 推荐你展示的对象
    # =========================
    # (A) 最适合展示 improvement 的 chain：all5 epoch 多、平均 success 高
    top_good = by_chain.head(10)
    top_good.to_csv(OUT_DIR / "recommended_chains_for_display.csv", index=False)

    # (B) 每条 chain 最早达到 5/5 的 epoch
    first_all5 = (
        all5_df.groupby("chain")
        .agg(first_epoch_all5=("epoch", "min"),
             last_epoch_all5=("epoch", "max"),
             num_all5_epochs=("epoch", "count"))
        .reset_index()
        .sort_values("first_epoch_all5")
    )
    first_all5.to_csv(OUT_DIR / "chain_first_epoch_with_all5_success.csv", index=False)

    # =========================
    # 8. 画两个简单图
    # =========================
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

    # =========================
    # 9. 打印结果
    # =========================
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