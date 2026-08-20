import gzip
import re
import pandas as pd

# Session 17: the bonus test from the original plan, never implemented until
# now -- does restricting to ClinVar variants with documented non-computational
# evidence change anything about how these predictors are validated? This
# needed ClinVar's per-submission free-text data, which the project's main
# dataset never carried (variant_summary.txt.gz was filtered down in session 1
# and doesn't include VariationID's Description text -- that only lives in a
# separate per-submission file).
#
# Two ClinVar downloads make this tractable without touching the ~5GB nested
# XML files:
#   1. variant_summary.txt.gz (already on disk since session 1) has both
#      AlleleID and VariationID side by side, so it gives a free mapping
#      between the two without any new download.
#   2. submission_summary.txt.gz (387MB, much smaller than the full XML) has
#      one row per submission (SCV), keyed by VariationID, with a free-text
#      Description field. Checking real rows directly (not just the
#      documentation) showed submitters often write out literal ACMG evidence
#      codes in that field, e.g. "classified as Likely pathogenic based on
#      ACMG criteria: PVS1_vstrong, PM2_mod."
#
# Of the 28 standard ACMG/AMP criteria (PS1-4, PM1-6, PP1-5, BA1, BS1-4,
# BP1-7), exactly two -- PP3 and BP4 -- represent computational/in-silico
# prediction evidence. Every other code represents something else (functional
# studies, segregation, de novo occurrence, population frequency, etc.). A
# variant whose Description cites at least one non-PP3/BP4 code has
# documented non-computational evidence on record; this only works for
# submissions where the code was actually spelled out in the free text, which
# is a real subset, not all of them -- coverage is measured and reported
# below rather than assumed.

VARIANT_SUMMARY_PATH = "data/variant_summary.txt.gz"
SUBMISSION_SUMMARY_PATH = "data/submission_summary.txt.gz"
CLINVAR_PATH = "data/clinvar_complete.csv"
OUT_PATH = "data/clinvar_evidence_codes.csv"

# All 28 standard ACMG/AMP codes. PP3/BP4 are computational; everything else
# is not.
ACMG_CODES = [
    "PVS1", "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5",
    "BA1", "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
]
COMPUTATIONAL_CODES = {"PP3", "BP4"}
# \b treats underscore as a word character, so it silently misses ClinVar's
# common strength-modifier suffix style ("PM2_mod", "PVS1_vstrong") -- a
# submission writing exactly that would have produced zero matches. Using
# explicit lookaround on alphanumerics (not underscore) instead of \b so a
# trailing "_mod"/"_vstrong"/"_strong"/"_supporting" doesn't hide the code,
# while still not matching a code as a substring of a longer alnum run
# (e.g. "PP1" inside "PP12" or "XPP1").
CODE_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])(" + "|".join(ACMG_CODES) + r")(?![A-Za-z0-9])"
)


def build_allele_to_variation_map(needed_allele_ids):
    print("Building AlleleID -> VariationID map from variant_summary.txt.gz...")
    mapping = {}
    with gzip.open(VARIANT_SUMMARY_PATH, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        allele_idx = header.index("#AlleleID")
        variation_idx = header.index("VariationID")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            try:
                allele_id = int(fields[allele_idx])
            except ValueError:
                continue
            if allele_id in needed_allele_ids:
                mapping[allele_id] = fields[variation_idx]
    print(f"  Mapped {len(mapping):,} / {len(needed_allele_ids):,} needed AlleleIDs")
    return mapping


def extract_codes(description):
    if not description or description == "-":
        return set()
    return set(CODE_PATTERN.findall(description))


def scan_submission_summary(needed_variation_ids):
    print("\nStreaming submission_summary.txt.gz (filtering to needed VariationIDs)...")
    by_variation = {}
    n_rows = 0
    n_matched = 0
    with gzip.open(SUBMISSION_SUMMARY_PATH, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            n_rows += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            variation_id = fields[0]
            if variation_id not in needed_variation_ids:
                continue
            n_matched += 1
            description = fields[3]
            collection_method = fields[7]
            entry = by_variation.setdefault(variation_id, {
                "methods": set(), "codes": set(), "n_submissions": 0,
            })
            entry["methods"].add(collection_method)
            entry["codes"] |= extract_codes(description)
            entry["n_submissions"] += 1
            if n_rows % 2_000_000 == 0:
                print(f"  ...{n_rows:,} rows scanned, {n_matched:,} matched so far")
    print(f"Done: {n_rows:,} total submission rows scanned, {n_matched:,} matched a needed variant")
    return by_variation


def main():
    print("Loading main dataset to get needed AlleleIDs...")
    df = pd.read_csv(CLINVAR_PATH, usecols=["#AlleleID"], low_memory=False)
    needed_allele_ids = set(df["#AlleleID"].dropna().astype(int))
    print(f"Needed AlleleIDs: {len(needed_allele_ids):,}")

    allele_to_variation = build_allele_to_variation_map(needed_allele_ids)
    needed_variation_ids = set(allele_to_variation.values())

    by_variation = scan_submission_summary(needed_variation_ids)

    rows = []
    for allele_id, variation_id in allele_to_variation.items():
        entry = by_variation.get(variation_id)
        if entry is None:
            rows.append({"AlleleID": allele_id, "VariationID": variation_id,
                         "n_submissions": 0, "methods": "", "acmg_codes": "",
                         "has_noncomputational_evidence": None})
            continue
        codes = entry["codes"]
        non_comp_codes = codes - COMPUTATIONAL_CODES
        has_non_comp = len(non_comp_codes) > 0 if codes else None
        rows.append({
            "AlleleID": allele_id,
            "VariationID": variation_id,
            "n_submissions": entry["n_submissions"],
            "methods": "|".join(sorted(entry["methods"])),
            "acmg_codes": "|".join(sorted(codes)),
            "has_noncomputational_evidence": has_non_comp,
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT_PATH, index=False)
    print(f"\nSaved {OUT_PATH}")

    print("\nCoverage:")
    print(f"  Variants with >=1 submission found: {(out['n_submissions'] > 0).sum():,} / {len(out):,}")
    print(f"  Variants with >=1 ACMG code extracted from free text: {(out['acmg_codes'] != '').sum():,} / {len(out):,}")
    print(f"  Variants classified has_noncomputational_evidence=True:  {(out['has_noncomputational_evidence'] == True).sum():,}")
    print(f"  Variants classified has_noncomputational_evidence=False: {(out['has_noncomputational_evidence'] == False).sum():,}")
    print(f"  Variants with no code extracted (unknown):                {out['has_noncomputational_evidence'].isna().sum():,}")

    print("\nCollectionMethod value counts (may be multiple per variant, semicolon-joined methods set):")
    from collections import Counter
    method_counts = Counter()
    for methods in out["methods"]:
        if methods:
            for m in methods.split("|"):
                method_counts[m] += 1
    for method, count in method_counts.most_common():
        print(f"  {method:25s} {count:,}")


if __name__ == "__main__":
    main()
