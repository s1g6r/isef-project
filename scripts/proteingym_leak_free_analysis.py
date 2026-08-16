import glob
import os

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

# Session 11: H2 -- leak-free AUC for AlphaMissense and ESM-1b on
# ProteinGym's experimental DMS assays. This ground truth comes from lab
# measurements, not human clinical curation, so nothing here could have
# leaked into either predictor's training data the way ClinVar labels
# might have. If these leak-free AUCs are meaningfully lower than the
# same predictors' AUC on ClinVar, that's evidence the ClinVar numbers were
# inflated by something other than genuine predictive skill. If they're
# similar, that's evidence the predictors are genuinely good, not just
# circularly validated -- either answer is a real finding.
#
# Ground truth: DMS_score_bin (1 = high fitness/tolerated, 0 = low
# fitness/damaging -- each assay's own binarization, documented in the
# reference file). Flipped here to "damaging" = 1 - DMS_score_bin so it's
# oriented the same direction as "pathogenic" in the ClinVar analysis (1 =
# the predictor should flag this as bad).
#
# Score orientation: AlphaMissense's am_pathogenicity/am_score is already
# higher-is-more-damaging everywhere else in this project. ESM-1b's raw
# score is a log-likelihood ratio where MORE NEGATIVE means more damaging
# (same convention as esm1b_llr and the -7.5 threshold in rdd_analysis.py)
# -- negated here so higher-oriented-score = more damaging for both
# predictors, which is what roc_auc_score expects.
#
# Only single-substitution mutants are used (no ":" in the mutant string).
# AlphaMissense can only score single substitutions, so multi-mutant rows
# are dropped from both predictors' analysis to keep the comparison on
# exactly the same variant set.

REFERENCE_PATH = "data/proteingym/DMS_substitutions_reference.csv"
MAPPING_PATH = "data/proteingym/uniprot_entry_name_to_accession.csv"
SCORES_DIR = "data/proteingym/zero_shot_scores"
AM_PATH = "data/AlphaMissense_hg38.tsv.gz"
CLINVAR_PATH = "data/clinvar_complete.csv"
OUT_PATH = "data/proteingym_leak_free_joined.csv"


def load_alphamissense_lookup():
    print("Loading AlphaMissense genome-wide scores (large file, may take a minute)...")
    am = pd.read_csv(AM_PATH, sep="\t", comment="#",
                      names=["chrom", "pos", "ref", "alt", "genome", "uniprot_id",
                             "transcript_id", "protein_variant", "am_pathogenicity", "am_class"],
                      usecols=["uniprot_id", "protein_variant", "am_pathogenicity"],
                      dtype={"uniprot_id": str, "protein_variant": str, "am_pathogenicity": float})
    # Same multi-transcript duplication dbNSFP has -- a (uniprot, protein_variant)
    # pair can appear more than once across transcripts, take the max.
    am = am.groupby(["uniprot_id", "protein_variant"], as_index=False)["am_pathogenicity"].max()
    print(f"  {len(am):,} unique (uniprot, protein_variant) AlphaMissense entries loaded")
    return am


def load_proteingym_scores(ref, entry_to_accession):
    human_ref = ref[ref["taxon"].str.lower() == "human"]
    # ProteinGym's "UniProt_ID" reference column is actually the mnemonic
    # entry name (e.g. "PAI1_HUMAN"), not the accession AlphaMissense uses
    # as its join key -- translate through the mapping from
    # get_uniprot_mapping.py before building the filename lookup.
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

        df = pd.read_csv(path, usecols=["mutant", "DMS_score_bin", "ESM1b"], low_memory=False)
        df = df[~df["mutant"].str.contains(":", na=False)].copy()
        if df.empty:
            continue

        df["uniprot_id"] = uniprot_id
        df["dms_filename"] = filename
        df["damaging"] = 1 - df["DMS_score_bin"]
        df["esm1b_oriented"] = -df["ESM1b"]
        rows.append(df[["uniprot_id", "dms_filename", "mutant", "damaging", "esm1b_oriented"]])

        if i % 20 == 0:
            print(f"  ...{i}/{len(files)} files read")

    combined = pd.concat(rows, ignore_index=True)
    print(f"Total single-mutant rows across all human-taxon assays: {len(combined):,}")
    return combined


def clinvar_baseline_aucs():
    """Same two predictors, computed the same way (raw continuous score vs.
    binary label, no threshold), on the current ClinVar dataset -- so the
    leak-free comparison below is apples-to-apples instead of comparing
    against session 1's older, differently-computed AlphaMissense number."""
    print("\nLoading current ClinVar dataset for a same-method comparison baseline...")
    df = pd.read_csv(CLINVAR_PATH, low_memory=False)

    results = {}
    am_valid = df.dropna(subset=["am_score", "label"])
    results["AlphaMissense"] = (roc_auc_score(am_valid["label"], am_valid["am_score"]), len(am_valid))

    esm_valid = df.dropna(subset=["esm1b_llr", "label"])
    results["ESM-1b"] = (roc_auc_score(esm_valid["label"], -esm_valid["esm1b_llr"]), len(esm_valid))
    return results


def main():
    ref = pd.read_csv(REFERENCE_PATH)
    mapping = pd.read_csv(MAPPING_PATH)
    entry_to_accession = dict(zip(mapping["entry_name"], mapping["accession"]))

    am_lookup = load_alphamissense_lookup()
    combined = load_proteingym_scores(ref, entry_to_accession)

    print("\nJoining AlphaMissense scores by (uniprot_id, protein_variant)...")
    combined = combined.merge(am_lookup, left_on=["uniprot_id", "mutant"],
                               right_on=["uniprot_id", "protein_variant"], how="left")
    n_am_matched = combined["am_pathogenicity"].notna().sum()
    print(f"  AlphaMissense matched: {n_am_matched:,} / {len(combined):,} "
          f"({100 * n_am_matched / len(combined):.1f}%)")

    combined.to_csv(OUT_PATH, index=False)
    print(f"Saved {OUT_PATH}")

    clinvar_aucs = clinvar_baseline_aucs()

    print("\n" + "=" * 60)
    print("Leak-free vs. ClinVar AUC comparison (H2)")
    print("=" * 60)

    esm_valid = combined.dropna(subset=["esm1b_oriented", "damaging"])
    if esm_valid["damaging"].nunique() == 2:
        pg_auc = roc_auc_score(esm_valid["damaging"], esm_valid["esm1b_oriented"])
        cv_auc, cv_n = clinvar_aucs["ESM-1b"]
        print(f"ESM-1b:         ProteinGym (leak-free) AUC={pg_auc:.4f} (n={len(esm_valid):,})  "
              f"vs.  ClinVar AUC={cv_auc:.4f} (n={cv_n:,})  "
              f"diff={pg_auc - cv_auc:+.4f}")

    am_valid = combined.dropna(subset=["am_pathogenicity", "damaging"])
    if am_valid["damaging"].nunique() == 2:
        pg_auc = roc_auc_score(am_valid["damaging"], am_valid["am_pathogenicity"])
        cv_auc, cv_n = clinvar_aucs["AlphaMissense"]
        print(f"AlphaMissense:  ProteinGym (leak-free) AUC={pg_auc:.4f} (n={len(am_valid):,})  "
              f"vs.  ClinVar AUC={cv_auc:.4f} (n={cv_n:,})  "
              f"diff={pg_auc - cv_auc:+.4f}")


if __name__ == "__main__":
    main()
