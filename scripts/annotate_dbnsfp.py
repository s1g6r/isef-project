import glob
import os
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

# Maps REVEL and Polyphen-2 scores (predictors #3 and #4) onto the main
# dataset from the dbNSFP v4.1a files downloaded by get_dbnsfp.py. Both come
# from the same source in one pass, since dbNSFP includes both columns.
#
# dbNSFP's primary chr/pos columns are hg38-based ("chr", "pos(1-based)"),
# matching the GRCh38 coordinates already used as the join key everywhere
# else in this project (ClinVar's PositionVCF/ReferenceAlleleVCF/
# AlternateAlleleVCF, via join_cv_am.py) -- no build conversion needed.
#
# dbNSFP quirk: a variant can map to multiple transcripts/isoforms, and some
# score columns pack all of them into one semicolon-separated string (e.g.
# "0.123;.;0.456"). Resolved by taking the max valid score per variant --
# the most-damaging prediction across transcripts, a common convention when
# collapsing transcript-level annotations to one variant-level call.
#
# Polyphen2 has two trained variants: HDIV (tuned for rare Mendelian disease
# variants) and HVAR (tuned to separate disease-causing from all remaining
# human variation). HVAR is the one Adzhubei et al. recommend for diagnostic-
# style classification, which is a closer match to what ClinVar P/LP/B/LB
# labels represent than HDIV -- used as the primary Polyphen2 column here,
# with HDIV also kept for comparison.

DBNSFP_DIR = "data/dbnsfp"
IN_PATH = "data/clinvar_phase1_complete.csv"
OUT_PATH = "data/clinvar_phase1_and_2_complete.csv"

USECOLS = ["#chr", "pos(1-based)", "ref", "alt",
           "REVEL_score", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score"]
SCORE_COLS = ["REVEL_score", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score"]
CHUNKSIZE = 1_000_000


def resolve_multi_one(value):
    """dbNSFP packs multi-transcript scores as 'v1;v2;.' -- take the max valid one.
    Only called on the (small) subset of values that actually contain ';' --
    see resolve_multi_column below."""
    parts = [p for p in value.split(";") if p not in (".", "")]
    return max(float(p) for p in parts) if parts else np.nan


def resolve_multi_column(s):
    """Vectorized fast path for the common case (single value, no ';'), with
    resolve_multi_one only applied to the rare multi-transcript rows. dbNSFP
    has ~80-90 million rows genome-wide, so a plain row-by-row .apply() over
    every row is slow -- pd.to_numeric handles the bulk of rows at C speed."""
    has_semi = s.str.contains(";", na=False)
    result = pd.to_numeric(s.where(~has_semi), errors="coerce")
    if has_semi.any():
        result.loc[has_semi] = s[has_semi].apply(resolve_multi_one)
    return result


def load_chrom_file(path, target_keys):
    """Streams the file in chunks instead of loading it whole. dbNSFP
    enumerates every possible amino acid substitution genome-wide (tens of
    millions of rows per chromosome file) -- we only need the ~213K variants
    already in our dataset, so each chunk is filtered down to matches before
    the next chunk is even read. Loading everything unfiltered previously
    blew up to 50+GB of memory and froze the machine."""
    matched_chunks = []
    n_chunks = 0
    for chunk in pd.read_csv(path, sep="\t", usecols=USECOLS, dtype=str,
                              chunksize=CHUNKSIZE, low_memory=False):
        n_chunks += 1
        chunk = chunk.rename(columns={"#chr": "chr", "pos(1-based)": "pos"})
        chunk["chr"] = chunk["chr"].str.replace("chr", "", regex=False)
        keys = chunk["chr"] + "|" + chunk["pos"] + "|" + chunk["ref"] + "|" + chunk["alt"]
        mask = keys.isin(target_keys)
        if mask.any():
            matched = chunk[mask].copy()
            for col in SCORE_COLS:
                matched[col] = resolve_multi_column(matched[col])
            matched_chunks.append(matched[["chr", "pos", "ref", "alt"] + SCORE_COLS])
        if n_chunks % 5 == 0:
            print(f"    ...{n_chunks * CHUNKSIZE:,} rows scanned so far")

    if not matched_chunks:
        return pd.DataFrame(columns=["chr", "pos", "ref", "alt"] + SCORE_COLS)
    result = pd.concat(matched_chunks, ignore_index=True)
    # A variant can appear on more than one row (not just multiple values
    # packed into one cell) if it maps to multiple gene/transcript entries --
    # collapse to one row per variant the same way, taking the max score.
    return result.groupby(["chr", "pos", "ref", "alt"], as_index=False)[SCORE_COLS].max()


def main():
    files = sorted(glob.glob(os.path.join(DBNSFP_DIR, "dbNSFP4.1a_variant.chr*.gz")))
    print(f"Found {len(files)} dbNSFP chromosome files.")
    if not files:
        print("No files found -- run get_dbnsfp.py first.")
        return

    print("Loading main dataset...")
    main_df = pd.read_csv(IN_PATH, dtype={"chrom": str, "pos": str, "ref": str, "alt": str}, low_memory=False)
    print(f"Rows: {len(main_df)}")

    target_keys = set(main_df["chrom"] + "|" + main_df["pos"] + "|" + main_df["ref"] + "|" + main_df["alt"])
    print(f"Unique variant keys to look up: {len(target_keys)}")

    chunks = []
    for i, path in enumerate(files, 1):
        print(f"  [{i}/{len(files)}] scanning {os.path.basename(path)}...")
        matched = load_chrom_file(path, target_keys)
        print(f"    matched {len(matched)} variants in this file")
        chunks.append(matched)
    dbnsfp = pd.concat(chunks, ignore_index=True)
    print(f"Total dbNSFP variants matched across all files: {len(dbnsfp)}")

    dbnsfp = dbnsfp.rename(columns={"chr": "chrom"})
    merged = main_df.merge(dbnsfp, on=["chrom", "pos", "ref", "alt"], how="left")

    for col in ["REVEL_score", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score"]:
        n_matched = merged[col].notna().sum()
        print(f"{col}: {n_matched} / {len(merged)} ({100 * n_matched / len(merged):.1f}%) populated")

    merged.to_csv(OUT_PATH, index=False)
    print(f"\nSaved to {OUT_PATH}")

    print("\nSanity-check AUCs vs published benchmarks:")
    for col, name, published in [
        ("REVEL_score", "REVEL", "0.90-0.96 (Ioannidis et al. 2016, varies by benchmark)"),
        ("Polyphen2_HVAR_score", "Polyphen2 (HVAR)", "~0.80-0.88 (Adzhubei et al. 2010, varies by benchmark)"),
        ("Polyphen2_HDIV_score", "Polyphen2 (HDIV)", "~0.85-0.92 (Adzhubei et al. 2010, varies by benchmark)"),
    ]:
        valid = merged.dropna(subset=[col, "label"])
        if len(valid) > 0 and valid["label"].nunique() == 2:
            auc = roc_auc_score(valid["label"], valid[col])
            print(f"  {name}: AUC={auc:.4f} (n={len(valid)}). Published: {published}")
        else:
            print(f"  {name}: not enough valid rows to compute AUC")


if __name__ == "__main__":
    main()
