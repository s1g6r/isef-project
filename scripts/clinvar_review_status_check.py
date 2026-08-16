import pandas as pd
from sklearn.metrics import roc_auc_score

# Session 11 wrap-up: one of the two open questions from today's ProteinGym
# writeup -- does restricting ClinVar to only its highest-confidence labels
# change the AUC gap to ProteinGym in an informative way?
#
# If restricting to "reviewed by expert panel" / "practice guideline" (the
# top two ClinVar review-confidence tiers) makes the ClinVar AUC noticeably
# HIGHER than the full-dataset AUC, that's consistent with the "ClinVar is
# an easier benchmark by composition" explanation from today -- expert-
# reviewed variants tend to be the clearest, most unambiguous cases, so an
# even bigger gap to ProteinGym's exhaustive, borderline-inclusive data
# would make sense without needing a leakage explanation at all.
# If the high-confidence AUC barely moves, that's more consistent with
# something else driving the gap -- though this check can't fully resolve
# leakage vs. composition either way on its own, it can only make one
# explanation more or less plausible.
#
# Matches review status by substring rather than an exact hardcoded list of
# categories, since the exact category strings weren't independently
# verified before writing this -- the printed value_counts() output lets
# that be checked directly against what actually matched.

CLINVAR_PATH = "data/clinvar_complete.csv"
HIGH_CONFIDENCE_PATTERNS = ["reviewed by expert panel", "practice guideline"]

# name, score column, orientation multiplier so higher = more pathogenic
# (am_score already is; esm1b_llr needs negating, same as everywhere else
# in this project)
PREDICTORS = [
    ("AlphaMissense", "am_score", 1, "0.7164"),
    ("ESM-1b", "esm1b_llr", -1, "0.6906"),
]


def main():
    df = pd.read_csv(CLINVAR_PATH, low_memory=False)

    print("Review status value counts (full dataset):")
    print(df["ReviewStatus"].value_counts())

    pattern = "|".join(HIGH_CONFIDENCE_PATTERNS)
    high_conf = df[df["ReviewStatus"].str.contains(pattern, case=False, na=False)]
    print(f"\n{len(high_conf):,} / {len(df):,} rows are high-confidence "
          f"({', '.join(HIGH_CONFIDENCE_PATTERNS)}).")

    print("\n" + "=" * 60)
    print("ClinVar AUC: full dataset vs. high-confidence-only")
    print("=" * 60)

    for name, score_col, orient, pg_auc in PREDICTORS:
        full_valid = df.dropna(subset=[score_col, "label"])
        full_auc = roc_auc_score(full_valid["label"], orient * full_valid[score_col])

        hc_valid = high_conf.dropna(subset=[score_col, "label"])
        if hc_valid["label"].nunique() < 2:
            print(f"\n{name}: not enough high-confidence data with both labels present (n={len(hc_valid)})")
            continue
        hc_auc = roc_auc_score(hc_valid["label"], orient * hc_valid[score_col])

        print(f"\n{name}:")
        print(f"  Full ClinVar AUC:          {full_auc:.4f} (n={len(full_valid):,})")
        print(f"  High-confidence-only AUC:  {hc_auc:.4f} (n={len(hc_valid):,})")
        print(f"  Difference:                {hc_auc - full_auc:+.4f}")
        print(f"  ProteinGym leak-free AUC:  {pg_auc}  (from session 11, for reference)")


if __name__ == "__main__":
    main()
