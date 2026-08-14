import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt

# Session 5 extension: ordered-vs-IDR split for ESM-1b, mirroring session 2's
# disorder_split_plddt.py analysis for AlphaMissense. disorder_class comes
# from AlphaFold pLDDT, a property of protein structure, not of either
# predictor -- so the same column applies unchanged here.
#
# ESM-1b sign convention: esm1b_llr is scored so LOWER (more negative) = more
# pathogenic -- opposite of am_score/label. Flip sign once, up front, so
# "higher = more pathogenic" holds everywhere below.

IN_PATH = "data/clinvar_phase1_complete.csv"
OUT_PLOT = "outputs/auc_over_time_esm1b_disorder.png"
CUTOFF = pd.Timestamp("2021-08-01")


def main():
    print("Loading Phase 1 complete dataset...")
    df = pd.read_csv(IN_PATH, low_memory=False)
    df["LastEvaluated"] = pd.to_datetime(df["LastEvaluated"], errors="coerce")
    df["esm1b_score"] = -df["esm1b_llr"]
    df = df.dropna(subset=["esm1b_score", "label", "LastEvaluated", "disorder_class"])
    print(f"Working rows: {len(df)}")

    ordered = df[df["disorder_class"] == "ordered"]
    idr = df[df["disorder_class"] == "disordered"]

    print(f"\nOrdered variants: {len(ordered)}")
    print(f"IDR variants:     {len(idr)}")

    auc_ordered = roc_auc_score(ordered["label"], ordered["esm1b_score"])
    auc_idr = roc_auc_score(idr["label"], idr["esm1b_score"])

    print(f"\nAUC (ordered): {auc_ordered:.4f}")
    print(f"AUC (IDR):     {auc_idr:.4f}")
    print(f"Gap:           {auc_ordered - auc_idr:.4f}")

    before = df[df["LastEvaluated"] < CUTOFF]
    after = df[df["LastEvaluated"] >= CUTOFF]

    print(f"\n--- Temporal split at ESM-1b release (Aug 2021) ---")
    print(f"Before: {len(before)} variants")
    print(f"After:  {len(after)} variants")
    for name, subset in [("Before", before), ("After", after)]:
        auc = roc_auc_score(subset["label"], subset["esm1b_score"])
        print(f"  AUC {name}: {auc:.4f}")

    # 2x2 breakdown: ordered/IDR x before/after -- this is the ESM-1b
    # equivalent of session 2's core H1 test for AlphaMissense.
    print(f"\n--- 2x2 breakdown (ordered/IDR x before/after) ---")
    results = {}
    for region_name, region_df in [("Ordered", ordered), ("IDR", idr)]:
        for time_name, is_after in [("Before", False), ("After", True)]:
            subset = region_df[region_df["LastEvaluated"] >= CUTOFF] if is_after \
                else region_df[region_df["LastEvaluated"] < CUTOFF]
            if len(subset) > 50 and subset["label"].nunique() == 2:
                auc = roc_auc_score(subset["label"], subset["esm1b_score"])
                results[(region_name, time_name)] = auc
                print(f"  {region_name} x {time_name}: n={len(subset)}, AUC={auc:.4f}")
            else:
                print(f"  {region_name} x {time_name}: n={len(subset)}, skipped (too small or single class)")

    if all(k in results for k in [("Ordered", "Before"), ("Ordered", "After"),
                                    ("IDR", "Before"), ("IDR", "After")]):
        ordered_jump = results[("Ordered", "After")] - results[("Ordered", "Before")]
        idr_jump = results[("IDR", "After")] - results[("IDR", "Before")]
        print(f"\nOrdered jump: {ordered_jump:+.4f}")
        print(f"IDR jump:     {idr_jump:+.4f}")

    # Plot AUC by year, ordered vs IDR -- same shape check session 2 used to
    # catch AlphaMissense/pLDDT's 2013-2021 convergence artifact.
    print("\nBuilding plot...")
    years = sorted(df["LastEvaluated"].dt.year.unique())
    aucs_ordered, aucs_idr = [], []
    for yr in years:
        o = ordered[ordered["LastEvaluated"].dt.year == yr]
        i = idr[idr["LastEvaluated"].dt.year == yr]
        aucs_ordered.append(roc_auc_score(o["label"], o["esm1b_score"])
                             if len(o) > 50 and o["label"].nunique() == 2 else np.nan)
        aucs_idr.append(roc_auc_score(i["label"], i["esm1b_score"])
                         if len(i) > 30 and i["label"].nunique() == 2 else np.nan)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(years, aucs_ordered, marker="o", label="Ordered (pLDDT >= 70)", color="steelblue")
    ax.plot(years, aucs_idr, marker="s", label="IDR (pLDDT < 50)", color="coral")
    ax.axvline(x=2021.58, color="black", linestyle="--", linewidth=1.5, label="ESM-1b release (Aug 2021)")
    ax.set_xlabel("Year of ClinVar LastEvaluated")
    ax.set_ylabel("AUC")
    ax.set_title("ESM-1b AUC vs ClinVar labels over time\nordered vs IDR (pLDDT)")
    ax.legend()
    ax.set_ylim(0.5, 1.0)
    plt.tight_layout()
    plt.savefig(OUT_PLOT, dpi=150)
    print(f"Plot saved to {OUT_PLOT}")
    print("\nDone.")


if __name__ == "__main__":
    main()
