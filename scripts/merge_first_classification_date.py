import pandas as pd

# Builds the final joined dataset: one clean file with ClinVar labels +
# AlphaMissense scores + ESM-1b scores + pLDDT disorder annotations +
# first-classification dates, all in one place.
#
# Join key is chrom/pos/ref/alt, same key used for every join so far in this
# project. clinvar_am_joined_plddt_esm1b.csv's chrom/pos/ref/alt originate from
# ClinVar's PositionVCF/ReferenceAlleleVCF/AlternateAlleleVCF (via join_cv_am.py).
# first_classification_date.csv uses the same VCF-style columns for 2021+
# snapshots and the equivalent Start/ReferenceAllele/AlternateAllele columns for
# pre-2021 snapshots (see build_first_classification_date.py) -- equivalent for
# the single-nucleotide missense variants this project covers, so the keys line up.

IN_PATH = "data/clinvar_am_joined_plddt_esm1b.csv"
DATES_PATH = "data/first_classification_date.csv"
OUT_PATH = "data/clinvar_phase1_complete.csv"


def main():
    print("Loading main joined dataset...")
    df = pd.read_csv(IN_PATH, dtype={"chrom": str, "pos": str, "ref": str, "alt": str}, low_memory=False)
    print(f"Rows: {len(df)}")

    print("Loading first-classification-date table...")
    dates = pd.read_csv(DATES_PATH, dtype={"chrom": str, "pos": str, "ref": str, "alt": str})
    print(f"Rows: {len(dates)}")

    merged = df.merge(dates, on=["chrom", "pos", "ref", "alt"], how="left")
    n_matched = merged["first_seen_release"].notna().sum()
    print(f"\nVariants with a matched first-classification date: {n_matched} / {len(merged)} "
          f"({100 * n_matched / len(merged):.1f}%)")

    merged.to_csv(OUT_PATH, index=False)
    print(f"Saved to {OUT_PATH}")

    print("\nColumn completeness in the final dataset:")
    for col, label in [("label", "ClinVar labels"), ("am_score", "AlphaMissense scores"),
                        ("esm1b_llr", "ESM-1b scores"), ("disorder_class", "pLDDT disorder annotations"),
                        ("first_seen_release", "First-classification dates")]:
        pct = 100 * merged[col].notna().sum() / len(merged)
        print(f"  {label}: {pct:.1f}% populated")


if __name__ == "__main__":
    main()
