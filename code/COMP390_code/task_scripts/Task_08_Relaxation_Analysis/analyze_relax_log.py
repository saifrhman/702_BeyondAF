from pathlib import Path
import os
import re

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

LOG_PATH = env_path("RELAX_LOG_PATH", COMP702_ROOT / "logs" / "slurm" / "of_relax_all.log")
OUT_DIR = env_path("RELAX_LOG_OUT_DIR", RESULT_ROOT / "relax_analysis_results")
OUT_DIR.mkdir(exist_ok=True, parents=True)

RUN_PAT = re.compile(
    r"epoch_(\d+).*?repeat_(\d+).*?/predictions/([^/]+?)_model_1_(unrelaxed|relaxed)\.pdb"
)


def parse_relax_log(log_path: Path, count_skipped_as_success: bool = False) -> pd.DataFrame:
    text = log_path.read_text(errors="ignore")
    records = {}
    current_key = None

    for line in text.splitlines():
        s = line.strip()

        if "Skip" in s and "relaxed" in s:
            path = s.split(":", 1)[1].strip() if ":" in s else s
            m = RUN_PAT.search(path)
            if m:
                epoch, repeat, chain, _ = m.groups()
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
                    rec["source_path"] = path

        elif "Relaxing:" in s:
            path = s.split("Relaxing:", 1)[1].strip()
            m = RUN_PAT.search(path)
            if m:
                epoch, repeat, chain, _ = m.groups()
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

        elif s.startswith("Saved:"):
            path = s.split("Saved:", 1)[1].strip()
            m = RUN_PAT.search(path)
            if m:
                epoch, repeat, chain, _ = m.groups()
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

        elif s.startswith("Violations:"):
            val = s.split(":", 1)[1].strip()
            try:
                val = float(val)
            except ValueError:
                pass
            if current_key in records:
                records[current_key]["violations"] = val

        elif s.startswith("ERROR processing:"):
            path = s.split("ERROR processing:", 1)[1].strip()
            m = RUN_PAT.search(path)
            if m:
                epoch, repeat, chain, _ = m.groups()
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

        elif (
            current_key is not None
            and current_key in records
            and records[current_key]["status"] == "error"
            and records[current_key]["error_message"] is None
            and s
            and not s.startswith("[")
        ):
            records[current_key]["error_message"] = s

    if not records:
        return pd.DataFrame(
            columns=[
                "epoch", "repeat", "chain", "status", "error_message",
                "violations", "source_path", "relax_error", "relax_success",
            ]
        )

    df = pd.DataFrame(records.values()).sort_values(["epoch", "repeat", "chain"]).reset_index(drop=True)
    df["relax_error"] = (df["status"] == "error").astype(int)

    if count_skipped_as_success:
        df["relax_success"] = df["status"].isin(["success", "skipped_existing_relaxed"]).astype(int)
    else:
        df["relax_success"] = (df["status"] == "success").astype(int)

    return df


def summarize(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=[group_col, "total", "success", "errors", "error_rate", "success_rate"])

    out = (
        df.groupby(group_col)
        .agg(total=("status", "size"), success=("relax_success", "sum"), errors=("relax_error", "sum"))
        .reset_index()
    )
    out["error_rate"] = out["errors"] / out["total"]
    out["success_rate"] = out["success"] / out["total"]
    return out


def plot_bar(summary_df: pd.DataFrame, x_col: str, y_col: str, title: str, out_path: Path):
    if summary_df.empty:
        return

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
    print("Reading log from:", LOG_PATH)
    print("Saving results to:", OUT_DIR)

    if not LOG_PATH.exists():
        raise FileNotFoundError(f"Relaxation log not found: {LOG_PATH}")

    df = parse_relax_log(LOG_PATH, count_skipped_as_success=False)

    sample_csv = OUT_DIR / "relax_sample_table.csv"
    df.to_csv(sample_csv, index=False)

    by_epoch = summarize(df, "epoch")
    by_chain = summarize(df, "chain")
    by_repeat = summarize(df, "repeat")

    by_epoch.to_csv(OUT_DIR / "summary_by_epoch.csv", index=False)
    by_chain.to_csv(OUT_DIR / "summary_by_chain.csv", index=False)
    by_repeat.to_csv(OUT_DIR / "summary_by_repeat.csv", index=False)

    plot_bar(by_epoch.sort_values("epoch"), "epoch", "error_rate", "Relax Error Rate by Epoch", OUT_DIR / "error_rate_by_epoch.png")
    plot_bar(by_chain.sort_values("error_rate", ascending=False), "chain", "error_rate", "Relax Error Rate by Protein Chain", OUT_DIR / "error_rate_by_chain.png")
    plot_bar(by_repeat.sort_values("repeat"), "repeat", "error_rate", "Relax Error Rate by Repeat", OUT_DIR / "error_rate_by_repeat.png")

    print("Saved sample-level table to:", sample_csv)
    print("Overall status counts:")
    print(df["status"].value_counts(dropna=False) if not df.empty else "No records found")


if __name__ == "__main__":
    main()
