import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

# Session 16: session 12's decile-stratification result (does the ClinVar-vs-
# ProteinGym AUC gap trace to composition-difficulty) was only ever printed
# to console, no figure. Building one now for the poster, reusing session
# 12's already-saved stratified data instead of rerunning the full pipeline.

IN_PATH = "data/proteingym_stratified.csv"
CLINVAR_PATH = "data/clinvar_complete.csv"
OUT_PATH = "outputs/proteingym_median_distance_deciles.png"

df = pd.read_csv(IN_PATH, low_memory=False)
cv = pd.read_csv(CLINVAR_PATH, low_memory=False)

esm_cv_auc = roc_auc_score(cv.dropna(subset=["esm1b_llr", "label"])["label"],
                            -cv.dropna(subset=["esm1b_llr", "label"])["esm1b_llr"])
am_cv_auc = roc_auc_score(cv.dropna(subset=["am_score", "label"])["label"],
                           cv.dropna(subset=["am_score", "label"])["am_score"])

strata_order = sorted(df["stratum"].dropna().unique(), key=lambda s: str(s))
# stratum labels from session 12 were "decile 1 (near median)" ... "decile 10
# (most extreme)" - sort numerically by the leading digit, not alphabetically
strata_order = sorted(df["stratum"].dropna().unique(),
                       key=lambda s: int(str(s).split()[1]))

esm_aucs, am_aucs, deciles = [], [], []
for i, stratum in enumerate(strata_order, 1):
    sub = df[df["stratum"] == stratum]
    esm_sub = sub.dropna(subset=["esm1b_oriented", "damaging"])
    am_sub = sub.dropna(subset=["am_pathogenicity", "damaging"])
    esm_aucs.append(roc_auc_score(esm_sub["damaging"], esm_sub["esm1b_oriented"]))
    am_aucs.append(roc_auc_score(am_sub["damaging"], am_sub["am_pathogenicity"]))
    deciles.append(i)

print("ESM-1b by decile:", [round(a, 4) for a in esm_aucs])
print("AlphaMissense by decile:", [round(a, 4) for a in am_aucs])

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.plot(deciles, am_aucs, marker="o", label="AlphaMissense (ProteinGym)", color="steelblue")
ax.plot(deciles, esm_aucs, marker="s", label="ESM-1b (ProteinGym)", color="coral")
ax.axhline(y=am_cv_auc, color="steelblue", linestyle="--", linewidth=1.2,
           label=f"AlphaMissense ClinVar AUC ({am_cv_auc:.3f})")
ax.axhline(y=esm_cv_auc, color="coral", linestyle="--", linewidth=1.2,
           label=f"ESM-1b ClinVar AUC ({esm_cv_auc:.3f})")
ax.axhline(y=0.5, color="gray", linestyle=":", linewidth=1)
ax.set_xlabel("Decile: distance from assay median fitness score\n(1 = most ambiguous, 10 = most extreme/unambiguous)")
ax.set_ylabel("AUC")
ax.set_xticks(deciles)
ax.set_ylim(0.45, 1.0)
ax.set_title("Composition-difficulty explains about half the ClinVar-ProteinGym gap")
ax.legend(fontsize=9, loc="lower right")
plt.tight_layout()
plt.savefig(OUT_PATH, dpi=300)
print(f"Saved {OUT_PATH}")
