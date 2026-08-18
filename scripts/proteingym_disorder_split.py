import glob
import os
import re

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

from proteingym_leak_free_analysis import (
    REFERENCE_PATH,
    MAPPING_PATH,
    SCORES_DIR,
    CLINVAR_PATH,
    load_alphamissense_lookup,
    load_proteingym_scores,
)

# Session 13: extend session 2's ordered-vs-IDR pLDDT split (originally built
# for ClinVar/AlphaMissense, then ClinVar/ESM-1b in session 5) onto the
# leak-free ProteinGym benchmark from sessions 11-12. Both predictors are
# known to score worse in disordered regions on ClinVar -- IDRs don't fold
# into a fixed structure, so structure-conditioned predictors have less to
# work with, and ClinVar itself under-curates IDR variants (few resolved
# structures, harder to classify pathogenicity), which is exactly the kind
# of composition effect sessions 11-12 were chasing on the DMS-score axis.
# This asks whether the ordered/IDR gap is a real property of the predictors
# (should reproduce on lab-measured fitness, not just clinical curation) or
# an artifact tied to how ClinVar itself is curated (should NOT reproduce on
# ProteinGym, since DMS assays exhaustively mutate every position in a
# protein regardless of how well-structured that region is).

CIF_DIR = "data/alphafold_cif"
OUT_PATH = "data/proteingym_plddt_joined.csv"


def parse_resnum(mutant):
    # ProteinGym's "mutant" strings are the same shape as ClinVar's aa_change
    # ("A673C" = residue 673) -- same regex as annotate_plddt.py's parse_pos.
    m = re.search(r"(\d+)", str(mutant))
    return int(m.group(1)) if m else None


def parse_plddt(cif_path):
    # Identical to annotate_plddt.py's parse_plddt: pLDDT lives in the
    # B_iso_or_equiv column of the CA atom line in AlphaFold DB mmCIF files.
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


def classify(p):
    if pd.isna(p):
        return None
    if p >= 70:
        return "ordered"
    if p >= 50:
        return "intermediate"
    return "disordered"


def clinvar_ordered_idr_baseline():
    """Same ordered/IDR AUC split, computed fresh on the current canonical
    ClinVar dataset, so the ProteinGym numbers below have an apples-to-apples
    comparison point instead of relying on older session 2/5 print output."""
    df = pd.read_csv(CLINVAR_PATH, low_memory=False)
    df = df.dropna(subset=["disorder_class"])
    results = {}
    for name, score_col, needs_flip in [("AlphaMissense", "am_score", False),
                                          ("ESM-1b", "esm1b_llr", True)]:
        sub = df.dropna(subset=[score_col, "label"])
        for region in ["ordered", "disordered"]:
            region_df = sub[sub["disorder_class"] == region]
            score = -region_df[score_col] if needs_flip else region_df[score_col]
            if len(region_df) > 50 and region_df["label"].nunique() == 2:
                results[(name, region)] = (roc_auc_score(region_df["label"], score), len(region_df))
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

    combined["resnum"] = combined["mutant"].apply(parse_resnum)

    print("\nAnnotating with pLDDT (one CIF parse per protein)...")
    plddt_values = pd.Series(index=combined.index, dtype="float64")
    n_no_cif, n_no_resnum, n_no_match = 0, 0, 0

    for uniprot_id, group in combined.groupby("uniprot_id"):
        cif_path = f"{CIF_DIR}/{uniprot_id}.cif"
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

    combined["plddt"] = plddt_values
    combined["disorder_class"] = combined["plddt"].apply(classify)

    print(f"\nVariants with no downloaded structure: {n_no_cif:,}")
    print(f"Variants with unparseable mutant string: {n_no_resnum:,}")
    print(f"Variants where residue number wasn't in the structure: {n_no_match:,}")
    print(f"Variants with a real pLDDT value: {combined['plddt'].notna().sum():,} / {len(combined):,}")
    print("\nDisorder class counts:")
    print(combined["disorder_class"].value_counts())
    print("\nDisorder class fractions (of annotated variants):")
    print(combined["disorder_class"].value_counts(normalize=True).round(4))

    combined.to_csv(OUT_PATH, index=False)
    print(f"\nSaved {OUT_PATH}")

    print("\n" + "=" * 70)
    print("Ordered vs. IDR AUC on ProteinGym (leak-free) vs. ClinVar baseline")
    print("=" * 70)

    clinvar_results = clinvar_ordered_idr_baseline()

    plot_data = {}
    for name, score_col, needs_flip in [("AlphaMissense", "am_pathogenicity", False),
                                          ("ESM-1b", "esm1b_oriented", False)]:
        sub = combined.dropna(subset=[score_col, "damaging", "disorder_class"])
        print(f"\n{name}:")
        pg_aucs = {}
        for region in ["ordered", "disordered"]:
            region_df = sub[sub["disorder_class"] == region]
            if len(region_df) > 50 and region_df["damaging"].nunique() == 2:
                auc = roc_auc_score(region_df["damaging"], region_df[score_col])
                pg_aucs[region] = auc
                cv_auc, cv_n = clinvar_results.get((name, region), (None, None))
                cv_str = f"AUC={cv_auc:.4f} (n={cv_n:,})" if cv_auc is not None else "n/a"
                print(f"  {region:12s} ProteinGym: n={len(region_df):,}  AUC={auc:.4f}   |   ClinVar: {cv_str}")
            else:
                print(f"  {region:12s} ProteinGym: n={len(region_df):,}  skipped (too small or single class)")
        if "ordered" in pg_aucs and "disordered" in pg_aucs:
            print(f"  ProteinGym ordered-IDR gap: {pg_aucs['ordered'] - pg_aucs['disordered']:+.4f}")
        if (name, "ordered") in clinvar_results and (name, "disordered") in clinvar_results:
            cv_gap = clinvar_results[(name, "ordered")][0] - clinvar_results[(name, "disordered")][0]
            print(f"  ClinVar ordered-IDR gap:    {cv_gap:+.4f}")
        plot_data[name] = {
            "ProteinGym ordered": pg_aucs.get("ordered"),
            "ProteinGym IDR": pg_aucs.get("disordered"),
            "ClinVar ordered": clinvar_results.get((name, "ordered"), (None,))[0],
            "ClinVar IDR": clinvar_results.get((name, "disordered"), (None,))[0],
        }

    print("\nBuilding plot...")
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    bar_labels = ["ClinVar\nordered", "ClinVar\nIDR", "ProteinGym\nordered", "ProteinGym\nIDR"]
    colors = ["steelblue", "coral", "steelblue", "coral"]
    for ax, name in zip(axes, ["AlphaMissense", "ESM-1b"]):
        vals = [plot_data[name]["ClinVar ordered"], plot_data[name]["ClinVar IDR"],
                plot_data[name]["ProteinGym ordered"], plot_data[name]["ProteinGym IDR"]]
        bars = ax.bar(bar_labels, vals, color=colors)
        for bar, v in zip(bars, vals):
            if v is not None:
                ax.text(bar.get_x() + bar.get_width() / 2, v + 0.01, f"{v:.3f}",
                         ha="center", va="bottom", fontsize=9)
        ax.axhline(y=0.5, color="gray", linestyle=":", linewidth=1)
        ax.set_title(name)
        ax.set_ylim(0.4, 1.0)
    axes[0].set_ylabel("AUC")
    fig.suptitle("Ordered vs. IDR AUC: ClinVar (curated) vs. ProteinGym (leak-free, exhaustive)")
    plt.tight_layout()
    plt.savefig("outputs/proteingym_disorder_split.png", dpi=150)
    print("Plot saved to outputs/proteingym_disorder_split.png")

    print("\nDone.")


if __name__ == "__main__":
    main()
