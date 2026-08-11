# ISEF Lab Notebook
**Project:** Quantifying Circular Evidence in ClinVar Variant Classification  
**Student:** Sagar Raut  
**Started:** August 11, 2026  

> **ISEF rule:** Every entry must be dated and signed. Keep a physical copy of this notebook as well — print and sign each session. Digital notebook is a backup, not a replacement.

---

## How to use this file
- Add a new session block at the bottom every time you work on the project
- Never edit past entries — if you made an error, note it in the next session
- Commit this file to GitHub at the end of every session (`git add notebook.md && git commit -m "notebook: session #X"`)
- Print and sign physically at least once a week

---

═══════════════════════════════════════════  
**DATE:** August 11, 2026  
**SESSION:** #1  
**TIME:** [fill in start] – [fill in end]  
═══════════════════════════════════════════  

**GOAL FOR TODAY:**  
Set up the data pipeline and confirm the ClinVar–AlphaMissense join works correctly before beginning any novel analysis.

**BACKGROUND / REASONING:**  
The project tests whether variant effect predictors (VEPs) like AlphaMissense are circularly evaluated — labs use them to classify variants, then researchers benchmark them against those same classifications. Before testing this hypothesis, I needed to confirm I could correctly join the two core datasets on genomic coordinates.

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
| Published AUC (literature) | ~0.90–0.93 |

**OBSERVATIONS:**  
AUC of 0.9592 is ~0.03 above the published benchmark of ~0.90–0.93. Published numbers were measured on held-out sets designed to reduce label leakage. My measurement is on all of ClinVar with no such controls. This excess is a preliminary signal consistent with circular evidence inflating apparent performance. **Recorded before any novel analysis was run.**

**PROBLEMS ENCOUNTERED:**  
1. Tried to run Python code directly in zsh terminal → `zsh: parse error near ')'`. Fix: save all code as .py files, run with `python3 filename.py`  
2. pandas not installed → ran `pip3 install pandas`  
3. AlphaMissense file had no header row — first data row read as column names → fixed with `header=None` and manually assigned column names  
4. Chromosome format mismatch → fixed with `.str.replace('chr', '')`  
5. Data files accidentally committed to git before .gitignore was set up → removed with `git rm -r --cached data/`  

**NEXT SESSION GOAL:**  
Add real disorder annotations using AlphaFold pLDDT scores, then run temporal split (pre/post Sept 19, 2023) to get first preview of the RDD signal.

**QUESTIONS TO RESEARCH:**  
- What exactly is pLDDT and why is <50 the standard IDR cutoff?  
- What is a regression discontinuity design mathematically?  
- What is PP3/BP4 in the ACMG variant classification guidelines?  

**Signed:** _________________ &emsp; **Date:** August 11, 2026

---

═══════════════════════════════════════════  
**DATE:** [fill in]  
**SESSION:** #2  
**TIME:** [fill in start] – [fill in end]  
═══════════════════════════════════════════  

**GOAL FOR TODAY:**  
[fill in]

**BACKGROUND / REASONING:**  
[fill in]

**WHAT I DID:**  
1.   
2.   
3.   

**DATA / RESULTS:**  
| Metric | Value |
|---|---|
| | |

**OBSERVATIONS:**  
[fill in]

**PROBLEMS ENCOUNTERED:**  
[fill in]

**NEXT SESSION GOAL:**  
[fill in]

**QUESTIONS TO RESEARCH:**  
- 

**Signed:** _________________ &emsp; **Date:** [fill in]

---

<!-- COPY THE SESSION BLOCK ABOVE FOR EACH NEW SESSION -->
