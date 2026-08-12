import pandas as pd
from sklearn.metrics import roc_auc_score

# Sanity check on the session 2 finding: outputs/auc_over_time_plddt.png shows the
# ordered/IDR AUC lines converging gradually across 2013-2021, well before
# AlphaMissense's Sept 2023 release. That means the full-history "before" window
# used in disorder_split_plddt.py mixes in years where IDR performance was still
# catching up, which could fake a jump that isn't really about the release date.
#
# Fix for this check only: narrow "before" to 2020-01-01 through the release date
# (post-convergence, pre-release) and compare to "after" 2023-09-19+. If the IDR
# jump survives on this narrower, fairer window, that's a stronger signal. If it
# vanishes, the original 2x2 result was mostly picking up the pre-existing drift.

print("Loading pLDDT-annotated data...")
df = pd.read_csv("data/clinvar_am_joined_plddt.csv")
df["LastEvaluated"] = pd.to_datetime(df["LastEvaluated"], errors="coerce")
df = df.dropna(subset=["am_score", "label", "LastEvaluated", "disorder_class"])

cutoff = pd.Timestamp("2023-09-19")
window_start = pd.Timestamp("2020-01-01")

before = df[(df["LastEvaluated"] >= window_start) & (df["LastEvaluated"] < cutoff)]
after = df[df["LastEvaluated"] >= cutoff]

print(f"\n--- Narrow window check: 2020-01-01 to {cutoff.date()} vs {cutoff.date()}+ ---")
print(f"Before (2020-2023): {len(before)} variants")
print(f"After (2023+):      {len(after)} variants")

print("\n2x2 breakdown (ordered/IDR x narrow before/after):")
results = {}
for region_name in ["ordered", "disordered"]:
    for time_name, subset_all in [("Before(2020-2023)", before), ("After(2023+)", after)]:
        subset = subset_all[subset_all["disorder_class"] == region_name]
        if len(subset) > 50 and subset["label"].nunique() == 2:
            auc = roc_auc_score(subset["label"], subset["am_score"])
            results[(region_name, time_name)] = auc
            print(f"  {region_name} x {time_name}: n={len(subset)}, AUC={auc:.4f}")
        else:
            print(f"  {region_name} x {time_name}: n={len(subset)}, skipped (too small or single class)")

if all(k in results for k in [("ordered", "Before(2020-2023)"), ("ordered", "After(2023+)"),
                                ("disordered", "Before(2020-2023)"), ("disordered", "After(2023+)")]):
    ordered_jump = results[("ordered", "After(2023+)")] - results[("ordered", "Before(2020-2023)")]
    idr_jump = results[("disordered", "After(2023+)")] - results[("disordered", "Before(2020-2023)")]
    print(f"\nOrdered jump (narrow window): {ordered_jump:+.4f}")
    print(f"IDR jump (narrow window):     {idr_jump:+.4f}")
    print("\nCompare to full-history 2x2 from disorder_split_plddt.py:")
    print("  Ordered jump (full history): -0.0007")
    print("  IDR jump (full history):     +0.0212")
    print("\nIf the narrow-window IDR jump is close to the full-history one, the signal")
    print("survives and isn't just pre-existing drift. If it shrinks toward ~0, the")
    print("original 2x2 was mostly picking up the 2013-2021 convergence trend.")
