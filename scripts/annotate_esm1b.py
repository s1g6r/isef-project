import io
import re
import zipfile
import pandas as pd
from sklearn.metrics import roc_auc_score

# Session 4: maps ESM-1b LLR scores onto clinvar_am_joined_plddt.csv, the same
# way annotate_plddt.py (session 2) mapped AlphaFold pLDDT onto the AlphaMissense
# join. Same pattern: group variants by protein, parse each protein's file once,
# not once per variant.
#
# Each protein's file inside ALL_hum_isoforms_ESM1b_LLR.zip is named
# {uniprot_accession}_LLR.csv: rows = the 20 amino acids (possible mutants),
# columns = every residue position labeled "{wildtype_AA} {position}" (with a
# space -- confirmed from the source of the ESM-1b web portal that serves this
# same file: it strips the space before using the label as a lookup key).
# variant strings in our data (aa_change, e.g. "N430S") decompose as
# wildtype+position ("N430") + mutant amino acid ("S") -- so the lookup is
# matrix.loc[mutant_aa, "N430"].

ZIP_PATH = "data/esm1b/ALL_hum_isoforms_ESM1b_LLR.zip"
IN_PATH = "data/clinvar_am_joined_plddt.csv"
OUT_PATH = "data/clinvar_am_joined_plddt_esm1b.csv"
MISSING_PATH = "data/esm1b_missing_proteins.txt"

AA_CHANGE_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def build_member_index(zf):
    """Map uniprot accession -> full path inside the zip, without scanning
    the whole namelist once per protein."""
    index = {}
    for name in zf.namelist():
        if name.endswith("_LLR.csv"):
            base = name.rsplit("/", 1)[-1]
            accession = base[:-len("_LLR.csv")]
            index[accession] = name
    return index


def load_matrix(zf, member_name):
    with zf.open(member_name) as f:
        df = pd.read_csv(io.BytesIO(f.read()), index_col=0)
    df.columns = [c.strip().replace(" ", "") for c in df.columns]
    return df


def main():
    print("Loading joined dataset...")
    variants = pd.read_csv(IN_PATH, dtype={"uniprot": str, "aa_change": str}, low_memory=False)
    print(f"Total variants: {len(variants)}")

    print(f"Opening {ZIP_PATH} (indexing members, not extracting)...")
    zf = zipfile.ZipFile(ZIP_PATH)
    member_index = build_member_index(zf)
    print(f"Proteins available in ESM-1b catalog: {len(member_index)}")

    unique_proteins = variants["uniprot"].dropna().unique()
    print(f"Unique proteins needed: {len(unique_proteins)}")

    scores = {}          # row index -> esm1b_llr
    missing_proteins = []
    parsed = 0

    for protein in unique_proteins:
        member_name = member_index.get(protein)
        if member_name is None:
            missing_proteins.append(protein)
            continue

        try:
            matrix = load_matrix(zf, member_name)
        except Exception as e:
            missing_proteins.append(f"{protein} (parse error: {e})")
            continue

        parsed += 1
        rows = variants.index[variants["uniprot"] == protein]
        for idx in rows:
            aa_change = variants.at[idx, "aa_change"]
            if not isinstance(aa_change, str):
                continue
            m = AA_CHANGE_RE.match(aa_change)
            if not m:
                continue
            wt, pos, mut = m.groups()
            col_key = f"{wt}{pos}"
            if col_key in matrix.columns and mut in matrix.index:
                scores[idx] = matrix.loc[mut, col_key]

        if parsed % 1000 == 0:
            print(f"  parsed {parsed}/{len(unique_proteins)} proteins, {len(scores)} scores mapped so far")

    zf.close()

    print(f"\nProteins parsed: {parsed}")
    print(f"Proteins missing from ESM-1b catalog: {len(missing_proteins)}")
    with open(MISSING_PATH, "w") as f:
        f.write("\n".join(str(p) for p in missing_proteins))
    print(f"Missing protein list saved to {MISSING_PATH}")

    variants["esm1b_llr"] = pd.Series(scores)
    n_scored = variants["esm1b_llr"].notna().sum()
    print(f"\nVariants with ESM-1b score: {n_scored} / {len(variants)} ({100 * n_scored / len(variants):.1f}%)")

    variants.to_csv(OUT_PATH, index=False)
    print(f"Saved to {OUT_PATH}")

    # Sanity check, same idea as session 1's AlphaMissense AUC reproduction:
    # ESM-1b's LLR is signed so that LOWER (more negative) = more pathogenic
    # (Brandes et al. 2023 use a -7.5 cutoff), which is the opposite convention
    # from am_score (higher = more pathogenic) and from ClinVar's label=1=pathogenic.
    # Flip the sign before computing AUC so higher = more pathogenic, matching
    # the rest of the pipeline.
    valid = variants.dropna(subset=["esm1b_llr", "label"])
    if len(valid) > 0 and valid["label"].nunique() == 2:
        auc = roc_auc_score(valid["label"], -valid["esm1b_llr"])
        print(f"\nOverall ESM-1b AUC vs ClinVar labels (n={len(valid)}): {auc:.4f}")
        print("Published benchmark (Brandes et al. 2023, ClinVar test set): 0.885")
    else:
        print("\nNot enough valid rows to compute a sanity-check AUC.")


if __name__ == "__main__":
    main()
