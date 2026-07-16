from pathlib import Path
import pandas as pd

# =========================
# 1. 输入文件
# =========================
RELAX_SAMPLE_CSV = Path("/users/sgmwu14/scratch/relax_analysis_results_by_files/relax_sample_table.csv")
PLDDT_CSV = Path("/users/sgmwu14/scratch/plddt_analysis_results/unrelaxed_plddt_summary.csv")

# =========================
# 2. 输出目录
# =========================
OUT_DIR = Path("/users/sgmwu14/scratch/best_display_selection_by_plddt")
OUT_DIR.mkdir(exist_ok=True, parents=True)

OUT_BEST = OUT_DIR / "best_relaxed_file_per_epoch_chain.csv"
OUT_TOP_EPOCHS = OUT_DIR / "recommended_epochs_for_display.csv"


def main():
    relax_df = pd.read_csv(RELAX_SAMPLE_CSV)
    plddt_df = pd.read_csv(PLDDT_CSV)

    # 只保留 unrelaxed pLDDT
    if "state" in plddt_df.columns:
        plddt_df = plddt_df[plddt_df["state"] == "unrelaxed"].copy()

    # 只保留需要列
    relax_keep = [
        "epoch", "repeat", "chain",
        "relax_success", "relax_error",
        "unrelaxed_path", "relaxed_path", "relaxed_exists"
    ]
    relax_keep = [c for c in relax_keep if c in relax_df.columns]
    relax_df = relax_df[relax_keep].copy()

    plddt_keep = [
        "epoch", "repeat", "chain",
        "mean_plddt", "min_plddt", "max_plddt", "num_residues", "pdb_path"
    ]
    plddt_keep = [c for c in plddt_keep if c in plddt_df.columns]
    plddt_df = plddt_df[plddt_keep].copy()

    # merge
    merged = relax_df.merge(
        plddt_df,
        on=["epoch", "repeat", "chain"],
        how="left",
        validate="one_to_one"
    )

    # 只考虑成功 relax 的样本
    success_df = merged[merged["relax_success"] == 1].copy()

    # 对每个 epoch + chain，选 mean_plddt 最高的那个 repeat
    best_idx = success_df.groupby(["epoch", "chain"])["mean_plddt"].idxmax()
    best_df = success_df.loc[best_idx].copy()

    best_df = best_df.sort_values(["chain", "epoch"]).reset_index(drop=True)

    # 保存完整结果
    best_df.to_csv(OUT_BEST, index=False)

    # 再做一个适合挑图的汇总：
    # 每条 chain 推荐哪几个 epoch 用来展示
    # 这里简单选：
    # - 最早成功 epoch
    # - 中间 epoch
    # - 最晚 epoch
    rec_rows = []

    for chain, sub in best_df.groupby("chain"):
        sub = sub.sort_values("epoch").reset_index(drop=True)
        n = len(sub)
        if n == 0:
            continue

        chosen = []
        chosen.append(sub.iloc[0])          # earliest
        chosen.append(sub.iloc[n // 2])     # middle
        chosen.append(sub.iloc[-1])         # latest

        # 去重（防止样本太少时 earliest/middle/latest 重复）
        chosen_df = pd.DataFrame(chosen).drop_duplicates(subset=["epoch", "chain"])

        for _, row in chosen_df.iterrows():
            rec_rows.append({
                "chain": row["chain"],
                "epoch": row["epoch"],
                "repeat": row["repeat"],
                "mean_plddt": row["mean_plddt"],
                "relaxed_path": row["relaxed_path"],
                "unrelaxed_path": row["unrelaxed_path"],
            })

    rec_df = pd.DataFrame(rec_rows).sort_values(["chain", "epoch"]).reset_index(drop=True)
    rec_df.to_csv(OUT_TOP_EPOCHS, index=False)

    print("Saved best file table to:", OUT_BEST)
    print("Saved recommended epochs table to:", OUT_TOP_EPOCHS)
    print()
    print("Example best selections:")
    print(best_df[["epoch", "chain", "repeat", "mean_plddt", "relaxed_path"]].head(20))


if __name__ == "__main__":
    main()