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

ROOT = env_path("OPENFOLD_RELAX_ANALYSIS_ROOT", RESULT_ROOT / "openfold_inference_runs5times")
OUT_DIR = env_path("RELAX_ANALYSIS_BY_FILES_DIR", RESULT_ROOT / "relax_analysis_results_by_files")
OUT_DIR.mkdir(exist_ok=True, parents=True)

EPOCH_PAT = re.compile(r"epoch_(\d+)")
REPEAT_PAT = re.compile(r"repeat_(\d+)")


def parse_unrelaxed_path(path: Path):
    path_str = str(path)
    epoch_match = EPOCH_PAT.search(path_str)
    repeat_match = REPEAT_PAT.search(path_str)

    if not epoch_match or not repeat_match:
        return None

    name = path.name
    if not name.endswith("_unrelaxed.pdb"):
        return None

    chain = name.replace("_model_1_unrelaxed.pdb", "").replace("_unrelaxed.pdb", "")

    return {
        "epoch": int(epoch_match.group(1)),
        "repeat": int(repeat_match.group(1)),
        "chain": chain,
    }


def build_sample_table(root: Path) -> pd.DataFrame:
    rows = []
    unrelaxed_files = sorted(root.rglob("*_unrelaxed.pdb"))
    print(f"Found {len(unrelaxed_files)} unrelaxed PDB files")

    for unrelaxed_path in unrelaxed_files:
        parsed = parse_unrelaxed_path(unrelaxed_path)
        if parsed is None:
            continue

        relaxed_path = Path(str(unrelaxed_path).replace("_unrelaxed.pdb", "_relaxed.pdb"))
        relaxed_exists = relaxed_path.exists()

        rows.append(
            {
                **parsed,
                "status": "success" if relaxed_exists else "error",
                "relax_success": int(relaxed_exists),
                "relax_error": int(not relaxed_exists),
                "unrelaxed_path": str(unrelaxed_path),
                "relaxed_path": str(relaxed_path),
                "relaxed_exists": relaxed_exists,
            }
        )

    if not rows:
        return pd.DataFrame(
            columns=[
                "epoch", "repeat", "chain", "status", "relax_success", "relax_error",
                "unrelaxed_path", "relaxed_path", "relaxed_exists",
            ]
        )

    return pd.DataFrame(rows).sort_values(["epoch", "repeat", "chain"]).reset_index(drop=True)


def summarize(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=[group_col, "total", "success", "errors", "error_rate", "success_rate"])

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
    print("Scanning root:", ROOT)
    print("Saving results to:", OUT_DIR)

    if not ROOT.exists():
        raise FileNotFoundError(f"Relaxation root not found: {ROOT}")

    df = build_sample_table(ROOT)

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
    print("Observed unrelaxed samples:", len(df))
    print("Overall status counts:")
    print(df["status"].value_counts(dropna=False) if not df.empty else "No samples found")


if __name__ == "__main__":
    main()
