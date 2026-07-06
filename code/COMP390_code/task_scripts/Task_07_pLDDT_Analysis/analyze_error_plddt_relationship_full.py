from pathlib import Path
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Optional statistical modelling
import statsmodels.formula.api as smf

# =========================
# 1. Centralised input files
# =========================
# These paths can be overridden using environment variables.
# By default, they assume the standard COMP702 Barkla2 layout:
# ~/COMP702_BeyondAF/code/COMP390_code/Result/...

COMP702_ROOT = Path(os.environ.get("COMP702_ROOT", Path.home() / "COMP702_BeyondAF"))
COMP390_ROOT = Path(os.environ.get("COMP390_ROOT", COMP702_ROOT / "code" / "COMP390_code"))
RESULT_ROOT = Path(os.environ.get("COMP702_RESULT_ROOT", COMP390_ROOT / "Result"))

PLDDT_OUT_DIR = Path(os.environ.get("PLDDT_OUT_DIR", RESULT_ROOT / "plddt_analysis_results"))
RELAX_ANALYSIS_BY_FILES_DIR = Path(
    os.environ.get("RELAX_ANALYSIS_BY_FILES_DIR", RESULT_ROOT / "relax_analysis_results_by_files")
)

PLDDT_CSV = Path(os.environ.get("UNRELAXED_PLDDT_CSV", PLDDT_OUT_DIR / "unrelaxed_plddt_summary.csv"))

RELAX_SAMPLE_CSV = Path(os.environ.get("RELAX_SAMPLE_CSV", RELAX_ANALYSIS_BY_FILES_DIR / "relax_sample_table.csv"))
RELAX_CHAIN_CSV = Path(os.environ.get("RELAX_CHAIN_CSV", RELAX_ANALYSIS_BY_FILES_DIR / "summary_by_chain.csv"))
RELAX_EPOCH_CSV = Path(os.environ.get("RELAX_EPOCH_CSV", RELAX_ANALYSIS_BY_FILES_DIR / "summary_by_epoch.csv"))
RELAX_REPEAT_CSV = Path(os.environ.get("RELAX_REPEAT_CSV", RELAX_ANALYSIS_BY_FILES_DIR / "summary_by_repeat.csv"))

# =========================
# 2. Centralised output directory
# =========================
OUT_DIR = Path(os.environ.get("ERROR_PLDDT_OUT_DIR", RESULT_ROOT / "error_plddt_analysis_results"))
OUT_DIR.mkdir(exist_ok=True, parents=True)

MERGED_CSV = OUT_DIR / "merged_error_plddt_table.csv"


def load_inputs():
    plddt_df = pd.read_csv(PLDDT_CSV)
    relax_sample_df = pd.read_csv(RELAX_SAMPLE_CSV)
    relax_chain_df = pd.read_csv(RELAX_CHAIN_CSV)
    relax_epoch_df = pd.read_csv(RELAX_EPOCH_CSV)
    relax_repeat_df = pd.read_csv(RELAX_REPEAT_CSV)

    return plddt_df, relax_sample_df, relax_chain_df, relax_epoch_df, relax_repeat_df


def standardize_and_merge(plddt_df, relax_sample_df):
    # Keep only the required columns
    plddt_keep = [
        "epoch", "repeat", "chain",
        "mean_plddt", "min_plddt", "max_plddt", "num_residues",
        "state", "pdb_path"
    ]
    plddt_keep = [c for c in plddt_keep if c in plddt_df.columns]
    plddt_df = plddt_df[plddt_keep].copy()

    # If a state column exists, keep only unrelaxed predictions
    if "state" in plddt_df.columns:
        plddt_df = plddt_df[plddt_df["state"] == "unrelaxed"].copy()

    relax_keep = [
        "epoch", "repeat", "chain",
        "status", "relax_success", "relax_error",
        "unrelaxed_path", "relaxed_path", "relaxed_exists"
    ]
    relax_keep = [c for c in relax_keep if c in relax_sample_df.columns]
    relax_sample_df = relax_sample_df[relax_keep].copy()

    merged = relax_sample_df.merge(
        plddt_df,
        on=["epoch", "repeat", "chain"],
        how="left",
        validate="one_to_one"
    )

    return merged


def summarize_with_plddt(df, group_col: str):
    out = (
        df.groupby(group_col)
        .agg(
            total=("relax_error", "size"),
            errors=("relax_error", "sum"),
            success=("relax_success", "sum"),
            mean_plddt=("mean_plddt", "mean"),
            median_plddt=("mean_plddt", "median"),
            std_plddt=("mean_plddt", "std"),
            min_plddt=("mean_plddt", "min"),
            max_plddt=("mean_plddt", "max"),
        )
        .reset_index()
    )
    out["error_rate"] = out["errors"] / out["total"]
    out["success_rate"] = out["success"] / out["total"]
    return out


def save_summaries(merged_df):
    by_epoch = summarize_with_plddt(merged_df, "epoch")
    by_chain = summarize_with_plddt(merged_df, "chain")
    by_repeat = summarize_with_plddt(merged_df, "repeat")

    by_epoch.to_csv(OUT_DIR / "summary_error_plddt_by_epoch.csv", index=False)
    by_chain.to_csv(OUT_DIR / "summary_error_plddt_by_chain.csv", index=False)
    by_repeat.to_csv(OUT_DIR / "summary_error_plddt_by_repeat.csv", index=False)

    return by_epoch, by_chain, by_repeat


def save_error_vs_plddt_summary(merged_df):
    grouped = (
        merged_df.groupby("relax_error")["mean_plddt"]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
    )
    grouped.to_csv(OUT_DIR / "summary_plddt_by_error.csv", index=False)
    return grouped


def plot_box_success_vs_error(df):
    plot_df = df.dropna(subset=["mean_plddt"]).copy()

    success_vals = plot_df.loc[plot_df["relax_error"] == 0, "mean_plddt"].values
    error_vals = plot_df.loc[plot_df["relax_error"] == 1, "mean_plddt"].values

    plt.figure(figsize=(7, 5))
    plt.boxplot([success_vals, error_vals], tick_labels=["success", "error"])
    plt.ylabel("Mean pLDDT")
    plt.title("Mean pLDDT by Relaxation Outcome")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "boxplot_plddt_success_vs_error.png", dpi=200)
    plt.close()


def plot_epoch_vs_mean_plddt(df):
    plot_df = (
        df.groupby("epoch")
        .agg(mean_plddt=("mean_plddt", "mean"))
        .reset_index()
        .sort_values("epoch")
    )

    plt.figure(figsize=(9, 5))
    plt.plot(plot_df["epoch"], plot_df["mean_plddt"], marker="o")
    plt.xlabel("Epoch")
    plt.ylabel("Mean pLDDT")
    plt.title("Mean pLDDT by Epoch")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "mean_plddt_by_epoch.png", dpi=200)
    plt.close()


def plot_chain_vs_mean_plddt(df):
    plot_df = (
        df.groupby("chain")
        .agg(mean_plddt=("mean_plddt", "mean"))
        .reset_index()
        .sort_values("mean_plddt", ascending=False)
    )

    plt.figure(figsize=(9, 5))
    plt.bar(plot_df["chain"], plot_df["mean_plddt"])
    plt.xlabel("Protein Chain")
    plt.ylabel("Mean pLDDT")
    plt.title("Mean pLDDT by Protein Chain")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "mean_plddt_by_chain.png", dpi=200)
    plt.close()


def plot_repeat_vs_mean_plddt(df):
    plot_df = (
        df.groupby("repeat")
        .agg(mean_plddt=("mean_plddt", "mean"))
        .reset_index()
        .sort_values("repeat")
    )

    plt.figure(figsize=(8, 5))
    plt.bar(plot_df["repeat"].astype(str), plot_df["mean_plddt"])
    plt.xlabel("Repeat")
    plt.ylabel("Mean pLDDT")
    plt.title("Mean pLDDT by Repeat")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "mean_plddt_by_repeat.png", dpi=200)
    plt.close()


def plot_epoch_vs_error_rate(df):
    plot_df = (
        df.groupby("epoch")
        .agg(error_rate=("relax_error", "mean"))
        .reset_index()
        .sort_values("epoch")
    )

    plt.figure(figsize=(9, 5))
    plt.plot(plot_df["epoch"], plot_df["error_rate"], marker="o")
    plt.xlabel("Epoch")
    plt.ylabel("Error Rate")
    plt.title("Relax Error Rate by Epoch")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "error_rate_by_epoch.png", dpi=200)
    plt.close()


def plot_chain_vs_error_rate(df):
    plot_df = (
        df.groupby("chain")
        .agg(error_rate=("relax_error", "mean"))
        .reset_index()
        .sort_values("error_rate", ascending=False)
    )

    plt.figure(figsize=(9, 5))
    plt.bar(plot_df["chain"], plot_df["error_rate"])
    plt.xlabel("Protein Chain")
    plt.ylabel("Error Rate")
    plt.title("Relax Error Rate by Protein Chain")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "error_rate_by_chain.png", dpi=200)
    plt.close()


def plot_repeat_vs_error_rate(df):
    plot_df = (
        df.groupby("repeat")
        .agg(error_rate=("relax_error", "mean"))
        .reset_index()
        .sort_values("repeat")
    )

    plt.figure(figsize=(8, 5))
    plt.bar(plot_df["repeat"].astype(str), plot_df["error_rate"])
    plt.xlabel("Repeat")
    plt.ylabel("Error Rate")
    plt.title("Relax Error Rate by Repeat")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "error_rate_by_repeat.png", dpi=200)
    plt.close()


def plot_scatter_plddt_vs_error(df):
    plot_df = df.dropna(subset=["mean_plddt"]).copy()

    rng = np.random.default_rng(42)
    jitter = rng.normal(0, 0.03, size=len(plot_df))
    y = plot_df["relax_error"].values + jitter

    plt.figure(figsize=(8, 5))
    plt.scatter(plot_df["mean_plddt"], y, alpha=0.5, s=18)
    plt.xlabel("Mean pLDDT")
    plt.ylabel("Relax Error (jittered)")
    plt.title("Mean pLDDT vs Relax Error")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "scatter_plddt_vs_error.png", dpi=200)
    plt.close()


def compute_correlations(df):
    corr_df = df.dropna(subset=["mean_plddt"]).copy()

    pearson_epoch_plddt = corr_df["epoch"].corr(corr_df["mean_plddt"], method="pearson")
    spearman_epoch_plddt = corr_df["epoch"].corr(corr_df["mean_plddt"], method="spearman")

    pearson_plddt_error = corr_df["mean_plddt"].corr(corr_df["relax_error"], method="pearson")
    spearman_plddt_error = corr_df["mean_plddt"].corr(corr_df["relax_error"], method="spearman")

    corr_out = pd.DataFrame({
        "metric": [
            "pearson(epoch, mean_plddt)",
            "spearman(epoch, mean_plddt)",
            "pearson(mean_plddt, relax_error)",
            "spearman(mean_plddt, relax_error)",
        ],
        "value": [
            pearson_epoch_plddt,
            spearman_epoch_plddt,
            pearson_plddt_error,
            spearman_plddt_error,
        ]
    })

    corr_out.to_csv(OUT_DIR / "correlation_summary.csv", index=False)
    return corr_out


def run_logistic_regression(df):
    model_df = df.dropna(subset=["mean_plddt"]).copy()

    model = smf.logit(
        formula="relax_error ~ mean_plddt + epoch + C(chain) + C(repeat)",
        data=model_df
    ).fit(disp=False)

    with open(OUT_DIR / "logistic_regression_summary.txt", "w") as f:
        f.write(model.summary().as_text())

    coef_df = pd.DataFrame({
        "term": model.params.index,
        "coef": model.params.values,
        "pvalue": model.pvalues.values,
    })
    coef_df.to_csv(OUT_DIR / "logistic_regression_coefficients.csv", index=False)

    return model


def compare_with_existing_summaries(relax_chain_df, relax_epoch_df, relax_repeat_df):
    # Save copies of the original summaries for comparison
    relax_chain_df.to_csv(OUT_DIR / "original_summary_by_chain.csv", index=False)
    relax_epoch_df.to_csv(OUT_DIR / "original_summary_by_epoch.csv", index=False)
    relax_repeat_df.to_csv(OUT_DIR / "original_summary_by_repeat.csv", index=False)


def main():
    print("Reading inputs:")
    print("  pLDDT file       :", PLDDT_CSV)
    print("  relax sample     :", RELAX_SAMPLE_CSV)
    print("  relax by chain   :", RELAX_CHAIN_CSV)
    print("  relax by epoch   :", RELAX_EPOCH_CSV)
    print("  relax by repeat  :", RELAX_REPEAT_CSV)
    print("Saving outputs to:")
    print("  ", OUT_DIR)
    print()

    plddt_df, relax_sample_df, relax_chain_df, relax_epoch_df, relax_repeat_df = load_inputs()

    compare_with_existing_summaries(relax_chain_df, relax_epoch_df, relax_repeat_df)

    merged_df = standardize_and_merge(plddt_df, relax_sample_df)
    merged_df.to_csv(MERGED_CSV, index=False)

    print("Merged table saved to:", MERGED_CSV)
    print("Merged rows:", len(merged_df))
    print("Missing mean_plddt rows:", merged_df["mean_plddt"].isna().sum())
    print()

    by_epoch, by_chain, by_repeat = save_summaries(merged_df)
    by_error = save_error_vs_plddt_summary(merged_df)
    corr_df = compute_correlations(merged_df)

    print("Summary of mean pLDDT by relax_error:")
    print(by_error)
    print()

    print("Summary by epoch:")
    print(by_epoch.head(10))
    print()

    print("Summary by chain:")
    print(by_chain.sort_values("error_rate", ascending=False))
    print()

    print("Summary by repeat:")
    print(by_repeat.sort_values("repeat"))
    print()

    print("Correlations:")
    print(corr_df)
    print()

    # Generate plots
    plot_box_success_vs_error(merged_df)
    plot_epoch_vs_mean_plddt(merged_df)
    plot_chain_vs_mean_plddt(merged_df)
    plot_repeat_vs_mean_plddt(merged_df)
    plot_epoch_vs_error_rate(merged_df)
    plot_chain_vs_error_rate(merged_df)
    plot_repeat_vs_error_rate(merged_df)
    plot_scatter_plddt_vs_error(merged_df)

    # Run regression
    try:
        model = run_logistic_regression(merged_df)
        print("Logistic regression completed.")
        print(model.summary())
    except Exception as e:
        print("WARNING: logistic regression failed:", e)

    print()
    print("All outputs saved in:", OUT_DIR)


if __name__ == "__main__":
    main()