import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

# Session 16: session 11's leak-free ClinVar-vs-ProteinGym comparison was only
# ever printed to console, no figure. Building one now for the poster,
# reusing the already-saved join from session 11 instead of re-loading the
# 64-million-row AlphaMissense genome-wide file.

LEAK_FREE_PATH = "data/proteingym_leak_free_joined.csv"
CLINVAR_PATH = "data/clinvar_complete.csv"
OUT_PATH = "outputs/proteingym_leak_free_comparison.png"

pg = pd.read_csv(LEAK_FREE_PATH, low_memory=False)
cv = pd.read_csv(CLINVAR_PATH, low_memory=False)

esm_pg = pg.dropna(subset=["esm1b_oriented", "damaging"])
esm_pg_auc = roc_auc_score(esm_pg["damaging"], esm_pg["esm1b_oriented"])

am_pg = pg.dropna(subset=["am_pathogenicity", "damaging"])
am_pg_auc = roc_auc_score(am_pg["damaging"], am_pg["am_pathogenicity"])

esm_cv = cv.dropna(subset=["esm1b_llr", "label"])
esm_cv_auc = roc_auc_score(esm_cv["label"], -esm_cv["esm1b_llr"])

am_cv = cv.dropna(subset=["am_score", "label"])
am_cv_auc = roc_auc_score(am_cv["label"], am_cv["am_score"])

print(f"AlphaMissense: ClinVar={am_cv_auc:.4f}  ProteinGym={am_pg_auc:.4f}")
print(f"ESM-1b:        ClinVar={esm_cv_auc:.4f}  ProteinGym={esm_pg_auc:.4f}")

fig, ax = plt.subplots(figsize=(7, 5.5))
labels = ["AlphaMissense", "ESM-1b"]
clinvar_vals = [am_cv_auc, esm_cv_auc]
pg_vals = [am_pg_auc, esm_pg_auc]

x = range(len(labels))
width = 0.35
bars1 = ax.bar([i - width / 2 for i in x], clinvar_vals, width, label="ClinVar (curated)", color="steelblue")
bars2 = ax.bar([i + width / 2 for i in x], pg_vals, width, label="ProteinGym (leak-free)", color="coral")

for bars in (bars1, bars2):
    for bar in bars:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.01, f"{h:.3f}", ha="center", va="bottom", fontsize=10)

ax.axhline(y=0.5, color="gray", linestyle=":", linewidth=1)
ax.set_xticks(list(x))
ax.set_xticklabels(labels)
ax.set_ylabel("AUC")
ax.set_ylim(0.4, 1.05)
ax.set_title("ClinVar accuracy does not transfer to leak-free lab data")
ax.legend(loc="lower center")
plt.tight_layout()
plt.savefig(OUT_PATH, dpi=300)
print(f"Saved {OUT_PATH}")
