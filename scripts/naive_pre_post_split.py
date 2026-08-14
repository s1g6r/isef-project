import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt

# Session 5 (weekly_plan.md Week 5): naive pre/post split for AlphaMissense
# AND ESM-1b together -- first look at temporal signal across two independent
# predictors instead of one. "Naive" = split on ClinVar's own LastEvaluated
# date, same variable session 2 used for the AlphaMissense-only split. This
# is the baseline Week 6's RDD (local linear regression each side of the
# cutoff) will replace with a rigorous jump test.
#
# ESM-1b sign convention: esm1b_llr is scored so LOWER (more negative) = more
# pathogenic (Brandes et al. 2023 use a -7.5 cutoff) -- opposite of am_score
# and opposite of label (1 = pathogenic). Flip sign once, up front, so every
# AUC/comparison downstream treats "higher = more pathogenic" consistently
# across both predictors.

IN_PATH = "data/clinvar_phase1_complete.csv"
OUT_PLOT = "outputs/auc_over_time_naive_split.png"

PREDICTORS = [
    {"name": "AlphaMissense", "score_col": "am_score",
     "cutoff": pd.Timestamp("2023-09-19"), "release_label": "Sept 2023"},
    {"name": "ESM-1b", "score_col": "esm1b_score",
     "cutoff": pd.Timestamp("2021-08-01"), "release_label": "Aug 2021"},
]


def auc_or_nan(subset, score_col, min_n=50):
    if len(subset) > min_n and subset["label"].nunique() == 2:
        return roc_auc_score(subset["label"], subset[score_col])
    return np.nan


def main():
    print("Loading Phase 1 complete dataset...")
    df = pd.read_csv(IN_PATH, low_memory=False)
    df["LastEvaluated"] = pd.to_datetime(df["LastEvaluated"], errors="coerce")
    df["esm1b_score"] = -df["esm1b_llr"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, pred in zip(axes, PREDICTORS):
        name, score_col, cutoff = pred["name"], pred["score_col"], pred["cutoff"]
        print(f"\n{'=' * 60}\n{name} (release: {pred['release_label']})\n{'=' * 60}")

        valid = df.dropna(subset=[score_col, "label", "LastEvaluated"])
        print(f"Valid rows: {len(valid)} / {len(df)}")

        before = valid[valid["LastEvaluated"] < cutoff]
        after = valid[valid["LastEvaluated"] >= cutoff]

        auc_before = auc_or_nan(before, score_col)
        auc_after = auc_or_nan(after, score_col)
        jump = auc_after - auc_before

        print(f"Before {pred['release_label']}: n={len(before)}, AUC={auc_before:.4f}")
        print(f"After  {pred['release_label']}: n={len(after)}, AUC={auc_after:.4f}")
        print(f"Naive jump: {jump:+.4f}")

        # Full-history year-by-year trend, to check for a pre-existing
        # convergence artifact before trusting the jump above (session 2
        # lesson: AlphaMissense/pLDDT had a 2013-2021 convergence trend that
        # inflated part of the ordered-vs-IDR 2x2 jump).
        years = sorted(valid["LastEvaluated"].dt.year.dropna().unique())
        aucs = [auc_or_nan(valid[valid["LastEvaluated"].dt.year == yr], score_col, min_n=30)
                for yr in years]

        print("\nYear-by-year AUC (full history, for confound check):")
        for yr, a in zip(years, aucs):
            print(f"  {yr}: {'n/a' if np.isnan(a) else f'{a:.4f}'}")

        ax.plot(years, aucs, marker="o", color="steelblue")
        ax.axvline(x=cutoff.year + (cutoff.month - 1) / 12, color="black", linestyle="--",
                   linewidth=1.5, label=f"release ({pred['release_label']})")
        ax.set_xlabel("Year of ClinVar LastEvaluated")
        ax.set_ylabel("AUC")
        ax.set_title(f"{name} AUC vs ClinVar labels over time")
        ax.set_ylim(0.5, 1.0)
        ax.legend()

    plt.tight_layout()
    plt.savefig(OUT_PLOT, dpi=150)
    print(f"\nPlot saved to {OUT_PLOT}")

    print("\n" + "=" * 60)
    print("Reminder: this is the NAIVE pre/post split (LastEvaluated cutoff).")
    print("Before trusting either jump as evidence of circularity, check above")
    print("whether the year-by-year trend was already rising/falling smoothly")
    print("through the release year -- that would mean part of the naive jump")
    print("is pre-existing drift, not a release-date discontinuity. Week 6's")
    print("RDD (local linear regression each side of the cutoff) is the real")
    print("test for this.")


if __name__ == "__main__":
    main()
