# ISEF Lab Notebook
**Project:** Quantifying Circular Evidence in ClinVar Variant Classification  
**Student:** Sagar Raut  
**Started:** August 11, 2026  

> **ISEF rule:** Every entry must be dated and signed. Keep a physical copy of this notebook as well - print and sign each session. Digital notebook is a backup, not a replacement.

---

═══════════════════════════════════════════  
**DATE:** August 11, 2026  
**SESSION:** #1  
**TIME:** 9:47 AM - 11:53 AM  
═══════════════════════════════════════════  

**GOAL FOR TODAY:**  
Set up the data pipeline and confirm the ClinVar-AlphaMissense join works correctly before beginning any novel analysis.

**BACKGROUND / REASONING:**  
The project tests whether variant effect predictors (VEPs) like AlphaMissense are circularly evaluated - labs use them to classify variants, then researchers benchmark them against those same classifications. Before testing this hypothesis, I needed to confirm I could correctly join the two core datasets on genomic coordinates.

**WHAT I DID:**  
1. Installed Python environment on MacBook, confirmed Python 3.9.6  
2. Downloaded `ClinVar variant_summary.txt.gz` from NCBI FTP → 9,036,351 total rows  
3. Downloaded `AlphaMissense_hg38.tsv.gz` from Zenodo/DeepMind → 71,697,556 rows  
4. Wrote `scripts/peek.py` → confirmed ClinVar column names; identified two issues: (a) both GRCh37 and GRCh38 in same file, (b) ClinicalSignificance contains compound labels like "Pathogenic/Likely pathogenic"  
5. Wrote `scripts/filter_grch38.py` → filtered to GRCh38 only, kept P/LP/B/LB only, parsed LastEvaluated as date, mapped to binary label (1=pathogenic, 0=benign) → 1,728,357 variants; 1,713,413 with valid date  
6. Wrote `scripts/join_cv_am.py` → identified chromosome format mismatch (ClinVar: "1", AlphaMissense: "chr1"), fixed by stripping "chr" prefix, joined on chrom + pos + ref + alt → 212,903 variants joined  

**DATA / RESULTS:**  
| Metric | Value |
|---|---|
| Total ClinVar rows | 9,036,351 |
| GRCh38 only | 4,484,398 |
| After P/LP/B/LB filter | 1,728,357 |
| With valid date | 1,713,413 |
| After join with AlphaMissense | 212,903 |
| Pathogenic (label=1) | 69,288 |
| Benign (label=0) | 143,615 |
| AlphaMissense AUC vs ClinVar | **0.9592** |
| Published AUC (literature) | ~0.90-0.93 |

**OBSERVATIONS:**  
AUC of 0.9592 is ~0.03 above the published benchmark of ~0.90-0.93. Published numbers were measured on held-out sets designed to reduce label leakage. My measurement is on all of ClinVar with no such controls. This excess is a preliminary signal consistent with circular evidence inflating apparent performance. **Recorded before any novel analysis was run.**

**PROBLEMS ENCOUNTERED:**  
1. Tried to run Python code directly in zsh terminal → `zsh: parse error near ')'`. Fix: save all code as .py files, run with `python3 filename.py`  
2. pandas not installed → ran `pip3 install pandas`  
3. AlphaMissense file had no header row - first data row read as column names → fixed with `header=None` and manually assigned column names  
4. Chromosome format mismatch → fixed with `.str.replace('chr', '')`  
5. Data files accidentally committed to git before .gitignore was set up → removed with `git rm -r --cached data/`  

**ADDITIONAL WORK COMPLETED THIS SESSION:**  
After the initial pipeline, ran first disorder split and temporal analysis.

7. Wrote `scripts/disorder_split.py` → split variants by disorder proxy (am_class == ambiguous), ran temporal split at AlphaMissense release date (Sept 19, 2023), produced first figure `outputs/auc_over_time.png`
8. Attempted MobiDB disorder annotation download via API → returned 0 rows, API endpoint appears broken. Will switch to AlphaFold pLDDT approach next session.

**FULL DATA / RESULTS:**  
| Metric | Value |
|---|---|
| Total ClinVar rows | 9,036,351 |
| GRCh38 only | 4,484,398 |
| After P/LP/B/LB filter | 1,728,357 |
| With valid date | 1,713,413 |
| After join with AlphaMissense | 212,903 |
| Pathogenic (label=1) | 69,288 |
| Benign (label=0) | 143,615 |
| AlphaMissense overall AUC | **0.9592** |
| Published AUC (literature) | ~0.90-0.93 |
| Working rows after date filter | 208,824 |
| Ordered variants (proxy) | 197,056 |
| IDR proxy variants (ambiguous class) | 11,768 |
| AUC ordered | 0.9668 |
| AUC IDR proxy | 0.5878 |
| Gap ordered vs IDR | 0.3791 |
| Variants before Sept 19 2023 | 77,593 |
| Variants after Sept 19 2023 | 131,231 |
| AUC before AlphaMissense release | 0.9586 |
| AUC after AlphaMissense release | 0.9608 |
| Temporal gap (overall) | 0.0022 |
| Ordered x Before AUC | 0.9646 |
| Ordered x After AUC | 0.9678 |
| IDR proxy x Before AUC | 0.5796 |
| IDR proxy x After AUC | 0.5927 |

**OBSERVATIONS:**  
Overall AUC of 0.9592 exceeds published benchmark (~0.90-0.93), consistent with circular evidence inflating performance. Recorded before novel analysis.

The IDR proxy (am_class == ambiguous) captures only 11,768 of 208,824 variants (5.6%), far below the expected ~30% of the proteome that is intrinsically disordered. This indicates the proxy is invalid - ambiguous AlphaMissense class does not correspond to structural disorder. The ordered/IDR AUC gap of 0.38 is directionally consistent with published literature but cannot be interpreted until real disorder annotations are applied.

The temporal gap (0.0022) is small overall. However, the IDR proxy shows a proportionally larger relative increase (0.5796 to 0.5927) compared to ordered regions (0.9646 to 0.9678). This is directionally consistent with the hypothesis but inconclusive until real disorder annotations replace the proxy.

First figure produced: `outputs/auc_over_time.png` - AUC plotted by year with AlphaMissense release marked.

**PROBLEMS ENCOUNTERED:**  
1. Tried to run Python code directly in zsh terminal → `zsh: parse error near ')'`. Fix: save all code as .py files, run with `python3 filename.py`  
2. pandas not installed → ran `pip3 install pandas`  
3. AlphaMissense file had no header row - first data row read as column names → fixed with `header=None` and manually assigned column names  
4. Chromosome format mismatch (ClinVar "1" vs AlphaMissense "chr1") → fixed with `.str.replace('chr', '')`  
5. Data files accidentally committed to git before .gitignore was set up → removed with `git filter-branch`, then force pushed  
6. MobiDB API returned empty dataframe despite correct column headers → endpoint broken, switching to AlphaFold pLDDT next session  

**NEXT SESSION GOAL:**  
Download AlphaFold pLDDT scores for human proteome, annotate each variant as ordered (pLDDT >= 70) or disordered (pLDDT < 50), rerun the temporal split with real disorder calls, and check whether the temporal signal becomes visible.

**QUESTIONS TO RESEARCH:**  
- What exactly is pLDDT and why is <50 the standard IDR cutoff?  
- What is a regression discontinuity design mathematically?  
- What is PP3/BP4 in the ACMG variant classification guidelines?  
- Why does the ambiguous class in AlphaMissense not correspond to structural disorder?  

**Signed:** _________________ &emsp; **Date:** August 11, 2026

---

