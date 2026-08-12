import os
import re
import pandas as pd

# Session 2: replace the invalid am_class=='ambiguous' IDR proxy (only caught 5.6%
# of variants, way below the expected ~30% disordered fraction) with real
# per-residue AlphaFold pLDDT scores.

CIF_DIR = "data/alphafold_cif"

print("Loading joined dataset...")
df = pd.read_csv("data/clinvar_am_joined.csv")
print(f"Total variants: {len(df)}")


def parse_pos(aa_change):
    # aa_change looks like "R330M" - pull out the residue number in the middle
    m = re.search(r"(\d+)", str(aa_change))
    return int(m.group(1)) if m else None


df["resnum"] = df["aa_change"].apply(parse_pos)


def parse_plddt(cif_path):
    # pLDDT is stored per-residue in the B_iso_or_equiv (B-factor) column of the
    # CA atom line. mmCIF _atom_site columns (fixed order for these AF files):
    # 0 group_PDB 1 id 2 type_symbol 3 label_atom_id ... 14 B_iso_or_equiv
    # ... 16 auth_seq_id (residue number matching UniProt numbering)
    lookup = {}
    with open(cif_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            fields = line.split()
            if len(fields) < 17 or fields[3] != "CA":
                continue
            try:
                resnum = int(fields[16])
                plddt = float(fields[14])
            except ValueError:
                continue
            lookup[resnum] = plddt
    return lookup


print("Annotating variants with pLDDT (one CIF parse per protein)...")
plddt_values = pd.Series(index=df.index, dtype="float64")
n_no_cif = 0
n_no_resnum = 0
n_no_match = 0

for uniprot, group in df.groupby("uniprot"):
    cif_path = f"{CIF_DIR}/{uniprot}.cif"
    if not os.path.exists(cif_path):
        n_no_cif += len(group)
        continue
    lookup = parse_plddt(cif_path)
    for idx, resnum in group["resnum"].items():
        if pd.isna(resnum):
            n_no_resnum += 1
            continue
        val = lookup.get(int(resnum))
        if val is None:
            n_no_match += 1
        else:
            plddt_values.at[idx] = val

df["plddt"] = plddt_values

print(f"\nVariants with no downloaded structure: {n_no_cif}")
print(f"Variants with unparseable aa_change: {n_no_resnum}")
print(f"Variants where residue number wasn't in the structure: {n_no_match}")
print(f"Variants with a real pLDDT value: {df['plddt'].notna().sum()} / {len(df)}")

# Standard pLDDT cutoffs used in Lin et al. 2025 and most IDR literature
def classify(p):
    if pd.isna(p):
        return None
    if p >= 70:
        return "ordered"
    if p >= 50:
        return "intermediate"
    return "disordered"

df["disorder_class"] = df["plddt"].apply(classify)

print("\nDisorder class counts:")
print(df["disorder_class"].value_counts())
print("\nDisorder class fractions (of annotated variants):")
print(df["disorder_class"].value_counts(normalize=True).round(4))
print("\nExpected proteome-wide disordered fraction is roughly ~30% - compare to the 'disordered' row above.")

df.to_csv("data/clinvar_am_joined_plddt.csv", index=False)
print("\nSaved to data/clinvar_am_joined_plddt.csv")
