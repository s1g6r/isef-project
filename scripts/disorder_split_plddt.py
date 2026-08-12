import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt

# Session 2: real ordered-vs-IDR split using AlphaFold pLDDT, replacing the
# invalid am_class=='ambiguous' proxy from session 1 (see disorder_split.py).

print("Loading pLDDT-annotated data...")
df = pd.read_csv("data/clinvar_am_joined_plddt.csv")
df["LastEvaluated"] = pd.to_datetime(df["LastEvaluated"], errors="coerce")
df = df.dropna(subset=["am_score", "label", "LastEvaluated", "disorder_class"])
print(f"Working rows: {len(df)}")

# Main comparison is ordered vs disordered (IDR). 'intermediate' (pLDDT 50-70)
# is excluded here since it's neither, same convention as the literature.
ordered = df[df["disorder_class"] == "ordered"]
idr = df[df["disorder_class"] == "disordered"]

print(f"\nOrdered variants: {len(ordered)}")
print(f"IDR variants:     {len(idr)}")

auc_ordered = roc_auc_score(ordered["label"], ordered["am_score"])
auc_idr = roc_auc_score(idr["label"], idr["am_score"])

print(f"\nAUC (ordered): {auc_ordered:.4f}")
print(f"AUC (IDR):     {auc_idr:.4f}")
print(f"Gap:           {auc_ordered - auc_idr:.4f}")
print("(Session 1 proxy gap was 0.3791, but the proxy only caught 5.6% of variants and was invalid)")

# Temporal split at AlphaMissense release date
cutoff = pd.Timestamp("2023-09-19")
before = df[df["LastEvaluated"] < cutoff]
after = df[df["LastEvaluated"] >= cutoff]

print(f"\n--- Temporal split at AlphaMissense release (Sept 19 2023) ---")
print(f"Before: {len(before)} variants")
print(f"After:  {len(after)} variants")
for name, subset in [("Before", before), ("After", after)]:
    auc = roc_auc_score(subset["label"], subset["am_score"])
    print(f"  AUC {name}: {auc:.4f}")

# 2x2 breakdown: ordered/IDR x before/after - this is the core H1 test
print(f"\n--- 2x2 breakdown (ordered/IDR x before/after) ---")
for region_name, region_df in [("Ordered", ordered), ("IDR", idr)]:
    for time_name, is_after in [("Before", False), ("After", True)]:
        subset = region_df[region_df["LastEvaluated"] >= cutoff] if is_after \
            else region_df[region_df["LastEvaluated"] < cutoff]
        if len(subset) > 50 and subset["label"].nunique() == 2:
            auc = roc_auc_score(subset["label"], subset["am_score"])
            print(f"  {region_name} x {time_name}: n={len(subset)}, AUC={auc:.4f}")
        else:
            print(f"  {region_name} x {time_name}: n={len(subset)}, skipped (too small or single class)")

# Plot AUC by year, ordered vs IDR
print("\nBuilding plot...")
years = sorted(df["LastEvaluated"].dt.year.unique())
aucs_ordered, aucs_idr = [], []

for yr in years:
    o = ordered[ordered["LastEvaluated"].dt.year == yr]
    i = idr[idr["LastEvaluated"].dt.year == yr]
    aucs_ordered.append(roc_auc_score(o["label"], o["am_score"]) if len(o) > 50 and o["label"].nunique() == 2 else np.nan)
    aucs_idr.append(roc_auc_score(i["label"], i["am_score"]) if len(i) > 30 and i["label"].nunique() == 2 else np.nan)

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(years, aucs_ordered, marker="o", label="Ordered (pLDDT >= 70)", color="steelblue")
ax.plot(years, aucs_idr, marker="s", label="IDR (pLDDT < 50)", color="coral")
ax.axvline(x=2023.72, color="black", linestyle="--", linewidth=1.5, label="AlphaMissense release (Sept 2023)")
ax.set_xlabel("Year of ClinVar LastEvaluated")
ax.set_ylabel("AUC")
ax.set_title("AlphaMissense AUC vs ClinVar labels over time\nordered vs IDR (real AlphaFold pLDDT)")
ax.legend()
ax.set_ylim(0.5, 1.0)
plt.tight_layout()
plt.savefig("outputs/auc_over_time_plddt.png", dpi=150)
print("Plot saved to outputs/auc_over_time_plddt.png")
print("\nDone.")
