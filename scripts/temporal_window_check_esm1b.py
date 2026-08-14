import pandas as pd
from sklearn.metrics import roc_auc_score

# Follow-up to disorder_split_esm1b.py: outputs/auc_over_time_esm1b_disorder.png
# shows ordered and IDR AUC converging noisily from roughly 2012 through 2021 --
# well before ESM-1b's Aug 2021 release. That's the same shape that made part
# of AlphaMissense's IDR jump (session 2, outputs/auc_over_time_plddt.png) a
# pre-existing convergence artifact rather than a true release-date
# discontinuity.
#
# This checks whether ESM-1b's IDR jump (+0.0545 full-history, from
# disorder_split_esm1b.py) survives on a narrower, post-convergence "before"
# window instead of the full 2010-2021 history.

IN_PATH = "data/clinvar_phase1_complete.csv"
CUTOFF = pd.Timestamp("2021-08-01")
WINDOW_START = pd.Timestamp("2019-01-01")


def main():
    print("Loading Phase 1 complete dataset...")
    df = pd.read_csv(IN_PATH, low_memory=False)
    df["LastEvaluated"] = pd.to_datetime(df["LastEvaluated"], errors="coerce")
    df["esm1b_score"] = -df["esm1b_llr"]
    df = df.dropna(subset=["esm1b_score", "label", "LastEvaluated", "disorder_class"])

    before = df[(df["LastEvaluated"] >= WINDOW_START) & (df["LastEvaluated"] < CUTOFF)]
    after = df[df["LastEvaluated"] >= CUTOFF]

    print(f"\n--- Narrow window check: {WINDOW_START.date()} to {CUTOFF.date()} vs {CUTOFF.date()}+ ---")
    print(f"Before ({WINDOW_START.date()} to 2021-08): {len(before)} variants")
    print(f"After (2021-08+):      {len(after)} variants")

    print("\n2x2 breakdown (ordered/IDR x narrow before/after):")
    results = {}
    for region_name in ["ordered", "disordered"]:
        for time_name, subset_all in [("Before(narrow)", before), ("After", after)]:
            subset = subset_all[subset_all["disorder_class"] == region_name]
            if len(subset) > 50 and subset["label"].nunique() == 2:
                auc = roc_auc_score(subset["label"], subset["esm1b_score"])
                results[(region_name, time_name)] = auc
                print(f"  {region_name} x {time_name}: n={len(subset)}, AUC={auc:.4f}")
            else:
                print(f"  {region_name} x {time_name}: n={len(subset)}, skipped (too small or single class)")

    if all(k in results for k in [("ordered", "Before(narrow)"), ("ordered", "After"),
                                    ("disordered", "Before(narrow)"), ("disordered", "After")]):
        ordered_jump = results[("ordered", "After")] - results[("ordered", "Before(narrow)")]
        idr_jump = results[("disordered", "After")] - results[("disordered", "Before(narrow)")]
        print(f"\nOrdered jump (narrow window): {ordered_jump:+.4f}")
        print(f"IDR jump (narrow window):     {idr_jump:+.4f}")
        print("\nCompare to full-history 2x2 from disorder_split_esm1b.py:")
        print("  Ordered jump (full history): +0.0300")
        print("  IDR jump (full history):     +0.0545")
        print("\nIf the narrow-window IDR jump is close to the full-history one, the signal")
        print("survives and isn't just pre-existing drift. If it shrinks toward the ordered")
        print("jump (or toward 0), the original 2x2 was mostly picking up the 2012-2021")
        print("convergence trend visible in the plot.")


if __name__ == "__main__":
    main()
