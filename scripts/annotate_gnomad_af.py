import glob
import os
import numpy as np
import pandas as pd

# Session 10: adds gnomAD population allele frequency to the canonical
# dataset -- one of the three confounders the plan calls for (the other two,
# gene and review status, are already columns: GeneSymbol and ReviewStatus).
# No new download needed -- dbNSFP (pulled in session 7 for REVEL/Polyphen2)
# bundles gnomAD_exomes_AF and gnomAD_genomes_AF as extra columns in the
# exact same files already sitting in data/dbnsfp/, so this is just another
# pass over data already on disk.
#
# Allele frequency is a real potential confound here: rarer variants are
# generally harder for any predictor to get right (less population data
# backing them), and if ClinVar's variant mix shifted toward rarer or more
# common variants right around a predictor's release date for reasons that
# have nothing to do with the predictor, that alone could produce or mask an
# RDD jump.
#
# gnomAD_exomes_AF is the primary column used downstream -- exomes has a
# much bigger sample size for coding/missense variants specifically, which
# is all this dataset contains. gnomAD_genomes_AF is kept as a secondary
# column in case it's useful later.

DBNSFP_DIR = "data/dbnsfp"
IN_PATH = "data/clinvar_phase1_and_2_complete.csv"
OUT_PATH = "data/clinvar_complete.csv"

USECOLS = ["#chr", "pos(1-based)", "ref", "alt", "gnomAD_exomes_AF", "gnomAD_genomes_AF"]
AF_COLS = ["gnomAD_exomes_AF", "gnomAD_genomes_AF"]
CHUNKSIZE = 1_000_000


def resolve_multi_one(value):
    """Same multi-transcript packing dbNSFP uses for REVEL/Polyphen2 scores
    (session 7-8) -- take the max valid value out of a 'v1;v2;.' string."""
    parts = [p for p in value.split(";") if p not in (".", "")]
    return max(float(p) for p in parts) if parts else np.nan


def resolve_multi_column(s):
    has_semi = s.str.contains(";", na=False)
    result = pd.to_numeric(s.where(~has_semi), errors="coerce")
    if has_semi.any():
        result.loc[has_semi] = s[has_semi].apply(resolve_multi_one)
    return result


def load_chrom_file(path, target_keys):
    """Same memory-safe streaming pattern as annotate_dbnsfp.py (session 8)
    -- filter each chunk down to the target variants before accumulating
    anything, instead of loading a whole genome-wide file into memory."""
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
            for col in AF_COLS:
                matched[col] = resolve_multi_column(matched[col])
            matched_chunks.append(matched[["chr", "pos", "ref", "alt"] + AF_COLS])
        if n_chunks % 5 == 0:
            print(f"    ...{n_chunks * CHUNKSIZE:,} rows scanned so far")

    if not matched_chunks:
        return pd.DataFrame(columns=["chr", "pos", "ref", "alt"] + AF_COLS)
    result = pd.concat(matched_chunks, ignore_index=True)
    return result.groupby(["chr", "pos", "ref", "alt"], as_index=False)[AF_COLS].max()


def main():
    files = sorted(glob.glob(os.path.join(DBNSFP_DIR, "dbNSFP4.1a_variant.chr*.gz")))
    print(f"Found {len(files)} dbNSFP chromosome files.")
    if not files:
        print("No files found -- dbNSFP should already be downloaded from session 7 (data/dbnsfp/).")
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
    gnomad = pd.concat(chunks, ignore_index=True)
    print(f"Total gnomAD AF rows matched across all files: {len(gnomad)}")

    gnomad = gnomad.rename(columns={"chr": "chrom"})
    merged = main_df.merge(gnomad, on=["chrom", "pos", "ref", "alt"], how="left")

    for col in AF_COLS:
        n_matched = merged[col].notna().sum()
        print(f"{col}: {n_matched} / {len(merged)} ({100 * n_matched / len(merged):.1f}%) populated")

    merged.to_csv(OUT_PATH, index=False)
    print(f"\nSaved to {OUT_PATH}")
    print("This is now the canonical dataset going forward (supersedes clinvar_phase1_and_2_complete.csv).")


if __name__ == "__main__":
    main()
