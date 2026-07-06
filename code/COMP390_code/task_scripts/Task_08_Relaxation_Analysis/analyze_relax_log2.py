from analyze_relax_log import LOG_PATH, OUT_DIR, parse_relax_log, summarize, plot_bar


def main():
    print("Reading log from:", LOG_PATH)
    print("Saving results to:", OUT_DIR)

    if not LOG_PATH.exists():
        raise FileNotFoundError(f"Relaxation log not found: {LOG_PATH}")

    df = parse_relax_log(LOG_PATH, count_skipped_as_success=True)

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
    print("Skipped existing relaxed files are counted as successful.")
    print("Overall status counts:")
    print(df["status"].value_counts(dropna=False) if not df.empty else "No records found")


if __name__ == "__main__":
    main()
