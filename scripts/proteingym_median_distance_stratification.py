import glob
import os

import pandas as pd
from sklearn.metrics import roc_auc_score

from proteingym_leak_free_analysis import (
    REFERENCE_PATH,
    MAPPING_PATH,
    SCORES_DIR,
    load_alphamissense_lookup,
)

# Session 12: session 11 found both predictors drop about 0.24 AUC on
# ProteinGym's leak-free data compared to ClinVar, but couldn't tell whether
# that's because of training-data leakage or because ClinVar's curated
# labels lean toward clear-cut cases while ProteinGym's exhaustive
# per-position mutagenesis includes a lot of borderline variants. The
# review-status check ruled out one version of the composition explanation
# (restricting ClinVar to its highest-confidence labels didn't help), but
# that only tested confidence, not how extreme the underlying biology is.
#
# This script tests the extremity angle directly: within each ProteinGym
# assay, how far is a variant's fitness score from that assay's own median?
# A variant right at the median is about as ambiguous as a fitness effect
# gets -- neither clearly damaging nor clearly tolerated. A variant way out
# in the tail is unambiguous. If predictor accuracy is much worse on the
# near-median variants than on the extreme ones, that's consistent with
# ClinVar's AUC being inflated mostly because it rarely contains anything
# that ambiguous in the first place -- a composition effect, not leakage.
# If accuracy stays roughly flat across the whole range, that's evidence
# against a simple composition explanation and leaves leakage (or
# something else) more in play.
#
# Percentile rank is used instead of raw distance from the median because
# assays differ wildly in what scale their fitness scores are measured on
# -- percentile makes "how extreme is this variant" comparable across all
# of them.

OUT_PATH = "data/proteingym_stratified.csv"
# Deciles instead of terciles: the tercile version showed AUC rising with
# distance from the median, but a rise-then-plateau pattern and a rise-all-
# the-way-to-the-edge pattern would have looked the same at only 3 points.
# Deciles show whether the AUC keeps climbing all the way to the extreme
# tail (composition explains most of the gap) or flattens out early
# (composition explains only part of it, before the tail hits some other
# ceiling -- e.g. leakage, or DMS fitness not mapping perfectly onto
# clinical pathogenicity).
N_STRATA = 10
STRATUM_LABELS = (["decile 1 (near median)"] + [f"decile {i}" for i in range(2, 10)]
                   + ["decile 10 (most extreme)"])


def load_proteingym_scores_with_dms_score(ref, entry_to_accession):
    human_ref = ref[ref["taxon"].str.lower() == "human"]
    uniprot_by_filename = {
        filename: entry_to_accession.get(entry_name)
        for filename, entry_name in zip(human_ref["DMS_filename"], human_ref["UniProt_ID"])
    }

    files = sorted(glob.glob(os.path.join(SCORES_DIR, "*.csv")))
    print(f"Found {len(files)} human-taxon score files.")

    rows = []
    for i, path in enumerate(files, 1):
        filename = os.path.basename(path)
        uniprot_id = uniprot_by_filename.get(filename)
        if uniprot_id is None:
            continue

        df = pd.read_csv(path, usecols=["mutant", "DMS_score", "DMS_score_bin", "ESM1b"],
                          low_memory=False)
        df = df[~df["mutant"].str.contains(":", na=False)].copy()
        if df.empty:
            continue

        df["uniprot_id"] = uniprot_id
        df["dms_filename"] = filename
        df["damaging"] = 1 - df["DMS_score_bin"]
        df["esm1b_oriented"] = -df["ESM1b"]
        # Percentile rank of this variant's raw fitness score within its
        # own assay -- 0.5 means right at the median, 0 or 1 means at
        # either extreme.
        df["dms_percentile"] = df.groupby("dms_filename")["DMS_score"].rank(pct=True)
        df["distance_from_median"] = (df["dms_percentile"] - 0.5).abs()

        rows.append(df[["uniprot_id", "dms_filename", "mutant", "damaging",
                         "esm1b_oriented", "distance_from_median"]])

        if i % 20 == 0:
            print(f"  ...{i}/{len(files)} files read")

    combined = pd.concat(rows, ignore_index=True)
    print(f"Total single-mutant rows across all human-taxon assays: {len(combined):,}")
    return combined


def auc_by_stratum(df, score_col, label, strata):
    valid = df.dropna(subset=[score_col, "damaging", "stratum"])
    print(f"\n{label}:")
    for stratum in strata:
        group = valid[valid["stratum"] == stratum]
        if group["damaging"].nunique() < 2:
            print(f"  {stratum:28s} n={len(group):,}  -- not enough of both classes to score")
            continue
        auc = roc_auc_score(group["damaging"], group[score_col])
        print(f"  {stratum:28s} n={len(group):,}  AUC={auc:.4f}")


def main():
    ref = pd.read_csv(REFERENCE_PATH)
    mapping = pd.read_csv(MAPPING_PATH)
    entry_to_accession = dict(zip(mapping["entry_name"], mapping["accession"]))

    am_lookup = load_alphamissense_lookup()
    combined = load_proteingym_scores_with_dms_score(ref, entry_to_accession)

    print("\nJoining AlphaMissense scores by (uniprot_id, protein_variant)...")
    combined = combined.merge(am_lookup, left_on=["uniprot_id", "mutant"],
                               right_on=["uniprot_id", "protein_variant"], how="left")
    n_am_matched = combined["am_pathogenicity"].notna().sum()
    print(f"  AlphaMissense matched: {n_am_matched:,} / {len(combined):,} "
          f"({100 * n_am_matched / len(combined):.1f}%)")

    # duplicates="drop" because pooling percentile-rank distances across 96
    # assays of very different sizes can put the same distance value right
    # on a bin boundary -- qcut would otherwise error out on that. Bin count
    # is checked afterward since dropping duplicate edges can leave fewer
    # than N_STRATA bins.
    raw_strata, bin_edges = pd.qcut(combined["distance_from_median"], N_STRATA,
                                     duplicates="drop", retbins=True)
    n_actual_strata = len(bin_edges) - 1
    if n_actual_strata == N_STRATA:
        combined["stratum"] = raw_strata.cat.rename_categories(STRATUM_LABELS)
        strata_order = STRATUM_LABELS
    else:
        print(f"Note: duplicate bin edges collapsed {N_STRATA} intended strata down to "
              f"{n_actual_strata} -- using the raw ranges as labels instead.")
        combined["stratum"] = raw_strata
        strata_order = list(raw_strata.cat.categories)

    combined.to_csv(OUT_PATH, index=False)
    print(f"Saved {OUT_PATH}")

    print("\n" + "=" * 60)
    print("AUC by distance from each assay's median fitness score")
    print("=" * 60)
    auc_by_stratum(combined, "esm1b_oriented", "ESM-1b", strata_order)
    auc_by_stratum(combined, "am_pathogenicity", "AlphaMissense", strata_order)


if __name__ == "__main__":
    main()
