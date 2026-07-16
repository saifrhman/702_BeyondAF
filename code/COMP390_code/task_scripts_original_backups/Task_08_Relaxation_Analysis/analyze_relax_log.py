from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt

# ===== 1. 修改成你的日志路径 =====
LOG_PATH = Path("/users/sgmwu14/fastscratch/of_relax_all_2789454.log")

# ===== 2. 输出目录 =====
OUT_DIR = Path("/users/sgmwu14/scratch/relax_analysis_results")
OUT_DIR.mkdir(exist_ok=True, parents=True)

# 匹配路径中的 epoch / repeat / chain
RUN_PAT = re.compile(
    r"epoch_(\d+)/your_model_epoch_\d+/repeat_(\d+)/predictions/([^/]+?)_model_1_(unrelaxed|relaxed)\.pdb"
)

def parse_relax_log(log_path: Path) -> pd.DataFrame:
    text = log_path.read_text(errors="ignore")
    lines = text.splitlines()

    records = {}
    current_key = None

    for line in lines:
        s = line.strip()

        # 1) 已存在 relaxed 文件
        if "Skip (already relaxed):" in s:
            path = s.split("Skip (already relaxed):", 1)[1].strip()
            m = RUN_PAT.search(path)
            if m:
                epoch, repeat, chain, kind = m.groups()
                key = (int(epoch), int(repeat), chain)
                rec = records.setdefault(
                    key,
                    {
                        "epoch": int(epoch),
                        "repeat": int(repeat),
                        "chain": chain,
                        "status": None,
                        "error_message": None,
                        "violations": None,
                        "source_path": path,
                    },
                )
                if rec["status"] is None:
                    rec["status"] = "skipped_existing_relaxed"

        # 2) 开始 relax
        elif "Relaxing:" in s:
            path = s.split("Relaxing:", 1)[1].strip()
            m = RUN_PAT.search(path)
            if m:
                epoch, repeat, chain, kind = m.groups()
                key = (int(epoch), int(repeat), chain)
                rec = records.setdefault(
                    key,
                    {
                        "epoch": int(epoch),
                        "repeat": int(repeat),
                        "chain": chain,
                        "status": None,
                        "error_message": None,
                        "violations": None,
                        "source_path": path,
                    },
                )
                rec["status"] = "running"
                rec["source_path"] = path
                current_key = key

        # 3) relax 成功保存
        elif s.startswith("Saved:"):
            path = s.split("Saved:", 1)[1].strip()
            m = RUN_PAT.search(path)
            if m:
                epoch, repeat, chain, kind = m.groups()
                key = (int(epoch), int(repeat), chain)
                rec = records.setdefault(
                    key,
                    {
                        "epoch": int(epoch),
                        "repeat": int(repeat),
                        "chain": chain,
                        "status": None,
                        "error_message": None,
                        "violations": None,
                        "source_path": path,
                    },
                )
                rec["status"] = "success"
                rec["source_path"] = path
                current_key = key

        # 4) violations
        elif s.startswith("Violations:"):
            val = s.split(":", 1)[1].strip()
            try:
                val = float(val)
            except ValueError:
                pass
            if current_key in records:
                records[current_key]["violations"] = val

        # 5) relax 失败
        elif s.startswith("ERROR processing:"):
            path = s.split("ERROR processing:", 1)[1].strip()
            m = RUN_PAT.search(path)
            if m:
                epoch, repeat, chain, kind = m.groups()
                key = (int(epoch), int(repeat), chain)
                rec = records.setdefault(
                    key,
                    {
                        "epoch": int(epoch),
                        "repeat": int(repeat),
                        "chain": chain,
                        "status": None,
                        "error_message": None,
                        "violations": None,
                        "source_path": path,
                    },
                )
                rec["status"] = "error"
                rec["source_path"] = path
                current_key = key

        # 6) error 的下一行通常是错误信息
        elif (
            current_key is not None
            and current_key in records
            and records[current_key]["status"] == "error"
            and records[current_key]["error_message"] is None
            and s
            and not s.startswith("[")
        ):
            records[current_key]["error_message"] = s

    df = pd.DataFrame(records.values()).sort_values(["epoch", "repeat", "chain"]).reset_index(drop=True)

    # 构造二分类标签
    df["relax_error"] = (df["status"] == "error").astype(int)
    df["relax_success"] = (df["status"] == "success").astype(int)

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
    df = parse_relax_log(LOG_PATH)

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

    # 基本打印
    print("Saved sample-level table to:", sample_csv)
    print()
    print("Overall status counts:")
    print(df["status"].value_counts(dropna=False))
    print()
    print("Top error messages:")
    print(df["error_message"].value_counts(dropna=False).head(10))
    print()
    print("Summary by epoch:")
    print(by_epoch.head())
    print()
    print("Summary by chain:")
    print(by_chain.sort_values("error_rate", ascending=False))

if __name__ == "__main__":
    main()