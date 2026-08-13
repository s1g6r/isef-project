import os
import pandas as pd

# Week 3 goal: build a first-classification-date table.
#
# LastEvaluated (used everywhere so far) is the date of the MOST RECENT review,
# which can get bumped even for an old variant that was reclassified/re-reviewed
# for unrelated reasons. That's a confound risk: a "post-AlphaMissense" LastEvaluated
# date doesn't necessarily mean the variant was actually FIRST classified after
# AlphaMissense existed.
#
# Fix: walk every archival ClinVar snapshot in chronological order and record,
# per variant (chrom,pos,ref,alt), the EARLIEST snapshot it already appears in
# with a P/LP/B/LB call. That's an approximate first-classification date,
# bucketed to ~6-month resolution (the spacing of the snapshots we have).
#
# Output: data/first_classification_date.csv
#   columns: chrom, pos, ref, alt, first_seen_release
# This gets joined against clinvar_am_joined_plddt.csv (on chrom/pos/ref/alt,
# same key already used for the AlphaMissense join) in a later session, for
# H3 / confounder control -- not done in this script.

KEEP_SIG = {
    'Pathogenic', 'Likely pathogenic', 'Benign', 'Likely benign',
    'Pathogenic/Likely pathogenic', 'Benign/Likely benign',
}

# Chronological order. First element of each tuple is the label used in the
# output; second is the file path. The 2019-01 file lives outside data/archive/
# because it was downloaded in session 1 before this script's layout existed.
# The final entry is the current full release, standing in for "most recent".
SNAPSHOTS = [
    ("2018-01", "data/archive/variant_summary_2018-01.txt.gz"),
    ("2018-07", "data/archive/variant_summary_2018-07.txt.gz"),
    ("2019-01", "data/variant_summary_2019-01.txt.gz"),
    ("2019-07", "data/archive/variant_summary_2019-07.txt.gz"),
    ("2020-01", "data/archive/variant_summary_2020-01.txt.gz"),
    ("2020-07", "data/archive/variant_summary_2020-07.txt.gz"),
    ("2021-01", "data/archive/variant_summary_2021-01.txt.gz"),
    ("2021-07", "data/archive/variant_summary_2021-07.txt.gz"),
    ("2022-01", "data/archive/variant_summary_2022-01.txt.gz"),
    ("2022-07", "data/archive/variant_summary_2022-07.txt.gz"),
    ("2023-01", "data/archive/variant_summary_2023-01.txt.gz"),
    ("2023-07", "data/archive/variant_summary_2023-07.txt.gz"),
    ("2024-01", "data/archive/variant_summary_2024-01.txt.gz"),
    ("2024-07", "data/archive/variant_summary_2024-07.txt.gz"),
    ("2025-01", "data/archive/variant_summary_2025-01.txt.gz"),
    ("2025-07", "data/archive/variant_summary_2025-07.txt.gz"),
    ("2026-01", "data/archive/variant_summary_2026-01.txt.gz"),
    ("2026-07", "data/archive/variant_summary_2026-07.txt.gz"),
    ("current", "data/variant_summary.txt.gz"),
]

# NCBI added PositionVCF/ReferenceAlleleVCF/AlternateAlleleVCF starting with the
# 2021-01 release. Older releases (2018 through 2020-07) only have Start/
# ReferenceAllele/AlternateAllele -- not VCF-normalized, but equivalent to the
# VCF columns for single-nucleotide substitutions, which covers the missense
# variants this project cares about. Try the VCF columns first, fall back to
# the legacy ones.
VCF_COLS = ['Chromosome', 'PositionVCF', 'ReferenceAlleleVCF', 'AlternateAlleleVCF']
LEGACY_COLS = ['Chromosome', 'Start', 'ReferenceAllele', 'AlternateAllele']


def load_keys(path):
    """Return the set of (chrom,pos,ref,alt) keys present in this snapshot
    as GRCh38 P/LP/B/LB variants."""
    header = pd.read_csv(path, sep='\t', nrows=0).columns
    if all(c in header for c in VCF_COLS):
        pos_col, ref_col, alt_col = 'PositionVCF', 'ReferenceAlleleVCF', 'AlternateAlleleVCF'
    else:
        pos_col, ref_col, alt_col = 'Start', 'ReferenceAllele', 'AlternateAllele'

    usecols = ['Assembly', 'Chromosome', pos_col, ref_col, alt_col, 'ClinicalSignificance']
    df = pd.read_csv(path, sep='\t', usecols=usecols, dtype=str, low_memory=False)
    df = df[df['Assembly'] == 'GRCh38']
    df = df[df['ClinicalSignificance'].isin(KEEP_SIG)]
    keys = (
        df['Chromosome'].astype(str) + '|' +
        df[pos_col].astype(str) + '|' +
        df[ref_col].astype(str) + '|' +
        df[alt_col].astype(str)
    )
    return set(keys)


def main():
    seen = {}  # key -> first_seen_release label
    bucket_counts = {}

    for label, path in SNAPSHOTS:
        if not os.path.exists(path):
            print(f"  {label}: MISSING ({path}) -- skipping. Run get_archival_clinvar.py first.")
            continue
        print(f"  {label}: loading {path} ...")
        try:
            keys = load_keys(path)
        except Exception as e:
            print(f"    FAILED to parse: {e} -- skipping this snapshot")
            continue

        new_keys = 0
        for k in keys:
            if k not in seen:
                seen[k] = label
                new_keys += 1
        bucket_counts[label] = new_keys
        print(f"    {len(keys)} GRCh38 P/LP/B/LB variants in file, {new_keys} newly first-seen here")

    print(f"\nTotal unique variants with a first-seen release: {len(seen)}")
    print("\nFirst-seen counts by release (validation -- should roughly track ClinVar's growth):")
    for label, _ in SNAPSHOTS:
        if label in bucket_counts:
            print(f"  {label}: {bucket_counts[label]}")

    out = pd.DataFrame(
        [(k.split('|')[0], k.split('|')[1], k.split('|')[2], k.split('|')[3], v)
         for k, v in seen.items()],
        columns=['chrom', 'pos', 'ref', 'alt', 'first_seen_release']
    )
    out.to_csv('data/first_classification_date.csv', index=False)
    print("\nSaved to data/first_classification_date.csv")


if __name__ == "__main__":
    main()
