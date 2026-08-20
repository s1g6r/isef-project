import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

# Session 17: the actual H3 test, now that evidence-code data exists
# (get_clinvar_evidence_codes.py). Sessions 6-10 already found no sudden
# accuracy jump at either predictor's release date, so this isn't about
# release dates -- it's a different, more direct cut at the same underlying
# circularity question: does AlphaMissense/ESM-1b agree unusually well with
# ClinVar specifically on variants whose classification was based only on
# PP3/BP4 (the two ACMG criteria that represent computational/in-silico
# prediction), compared to variants with documented non-computational
# evidence (functional studies, segregation, population data, etc.) on
# record?
#
# Important interpretive limit, worth stating up front rather than
# discovering after computing the numbers: ClinVar's free text usually
# doesn't say *which* computational tool a submitter's PP3/BP4 call was
# based on. A variant flagged BP4-only here could have been classified using
# SIFT or an older tool years before AlphaMissense existed. So this tests
# something related but not identical to "AlphaMissense specifically caused
# this" -- it tests whether being classified via computational evidence in
# general lines up with today's AlphaMissense/ESM-1b agreement.

CLINVAR_PATH = "data/clinvar_complete.csv"
EVIDENCE_PATH = "data/clinvar_evidence_codes.csv"
OUT_PLOT = "outputs/clinvar_evidence_auc_comparison.png"

cv = pd.read_csv(CLINVAR_PATH, low_memory=False)
ev = pd.read_csv(EVIDENCE_PATH)

df = cv.merge(ev, left_on="#AlleleID", right_on="AlleleID", how="inner")
print(f"Merged rows: {len(df):,} (from {len(cv):,} ClinVar rows, {len(ev):,} evidence rows)")


def bucket_label(v):
    if v == True:
        return "non-computational evidence documented"
    if v == False:
        return "computational (PP3/BP4) only"
    return "no ACMG code extracted (unknown)"


df["evidence_bucket"] = df["has_noncomputational_evidence"].apply(bucket_label)

print("\nLabel mix within each evidence bucket (pathogenic fraction):")
for bucket in ["computational (PP3/BP4) only", "non-computational evidence documented",
               "no ACMG code extracted (unknown)"]:
    sub = df[df["evidence_bucket"] == bucket]
    if len(sub) > 0:
        print(f"  {bucket:45s} n={len(sub):,}  pathogenic={sub['label'].mean():.3f}")

print("\n" + "=" * 78)
print("AlphaMissense / ESM-1b AUC by evidence bucket")
print("=" * 78)

results = {}
for name, score_col, needs_flip in [("AlphaMissense", "am_score", False),
                                      ("ESM-1b", "esm1b_llr", True)]:
    print(f"\n{name}:")
    results[name] = {}
    for bucket in ["computational (PP3/BP4) only", "non-computational evidence documented",
                   "no ACMG code extracted (unknown)"]:
        sub = df[df["evidence_bucket"] == bucket].dropna(subset=[score_col, "label"])
        if len(sub) > 50 and sub["label"].nunique() == 2:
            score = -sub[score_col] if needs_flip else sub[score_col]
            auc = roc_auc_score(sub["label"], score)
            results[name][bucket] = (auc, len(sub))
            print(f"  {bucket:45s} n={len(sub):,}  AUC={auc:.4f}")
        else:
            print(f"  {bucket:45s} n={len(sub):,}  skipped (too small or single class)")

print("\nDifference (computational-only minus non-computational-evidence):")
for name in ["AlphaMissense", "ESM-1b"]:
    r = results[name]
    if "computational (PP3/BP4) only" in r and "non-computational evidence documented" in r:
        diff = r["computational (PP3/BP4) only"][0] - r["non-computational evidence documented"][0]
        print(f"  {name}: {diff:+.4f}")

print("\nBuilding plot...")
fig, ax = plt.subplots(figsize=(8, 5.5))
buckets = ["non-computational evidence documented", "computational (PP3/BP4) only"]
bucket_short = ["Non-computational\nevidence documented", "Computational\n(PP3/BP4) only"]
x = range(len(buckets))
width = 0.35
am_vals = [results["AlphaMissense"].get(b, (None, None))[0] for b in buckets]
esm_vals = [results["ESM-1b"].get(b, (None, None))[0] for b in buckets]

bars1 = ax.bar([i - width / 2 for i in x], am_vals, width, label="AlphaMissense", color="steelblue")
bars2 = ax.bar([i + width / 2 for i in x], esm_vals, width, label="ESM-1b", color="coral")
for bars in (bars1, bars2):
    for bar, v in zip(bars, [am_vals[0], am_vals[1]] if bars is bars1 else [esm_vals[0], esm_vals[1]]):
        if v is not None:
            ax.text(bar.get_x() + bar.get_width() / 2, v + 0.01, f"{v:.3f}",
                     ha="center", va="bottom", fontsize=9)

ax.axhline(y=0.5, color="gray", linestyle=":", linewidth=1)
ax.set_xticks(list(x))
ax.set_xticklabels(bucket_short)
ax.set_ylabel("AUC")
ax.set_ylim(0.4, 1.05)
ax.set_title("ClinVar AUC by type of documented evidence (n=209,931 variants)")
ax.legend(loc="lower center")
plt.tight_layout()
plt.savefig(OUT_PLOT, dpi=300)
print(f"Saved {OUT_PLOT}")
print("\nDone.")
