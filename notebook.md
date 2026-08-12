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

**Signed:** Sagar Raut &emsp; **Date:** August 11, 2026

---

═══════════════════════════════════════════  
**DATE:** August 12, 2026  
**SESSION:** #2  
**TIME:** 9:05 AM - 10:51 AM  
═══════════════════════════════════════════  

**GOAL FOR TODAY:**  
Replace the invalid IDR proxy from session 1 with real disorder annotations (AlphaFold pLDDT), rerun the ordered-vs-IDR and temporal analysis, and sanity-check the result before trusting it.

**BACKGROUND / REASONING:**  
Session 1 ended with the am_class=='ambiguous' IDR proxy only capturing 5.6% of variants (expected ~30%), meaning it was not actually measuring structural disorder. The plan called for real AlphaFold pLDDT scores this session (Week 2 of weekly_plan.md).

**WHAT I DID:**  
1. Re-ran `scripts/get_disorder.py` (MobiDB download) to double check — confirmed the API is still broken: it returns HTTP 200 with correct column headers but zero data rows. Abandoned MobiDB permanently.
2. Modified `scripts/join_cv_am.py` to also keep the `uniprot` and `aa_change` columns from AlphaMissense (needed to map a variant to a residue position on the protein) → reran, reproduced the same 212,903 joined variants and 0.9592 AUC as session 1, confirming the join logic itself was unaffected.
3. Found unique UniProt accessions needed: 15,144.
4. Wrote `scripts/get_alphafold_plddt.py` → downloaded AlphaFold CIF structure files (one per protein) from the EBI AlphaFold DB (`AF-{accession}-F1-model_v6.cif`), concurrent + resumable → 14,930 downloaded (98.6%), 214 with no AlphaFold model (isoforms/obsolete accessions, expected), 0 network errors.
5. Wrote `scripts/annotate_plddt.py` → parsed per-residue pLDDT from the B-factor column of each CIF file (one parse per protein, not per variant), mapped to each variant via (uniprot, residue number parsed from aa_change), classified ordered (pLDDT >= 70) / intermediate (50-70) / disordered (<50) → produced `data/clinvar_am_joined_plddt.csv`.
6. Wrote `scripts/disorder_split_plddt.py` → real ordered-vs-IDR AUC split, temporal split, and 2x2 breakdown using the pLDDT annotations → produced `outputs/auc_over_time_plddt.png`.
7. Inspected the figure directly: the ordered/IDR AUC lines converge gradually across 2013-2021, well before the Sept 2023 AlphaMissense release — meaning a full-history pre/post comparison could be misleadingly picking up that pre-existing drift instead of a real discontinuity.
8. Wrote `scripts/temporal_window_check.py` to test this directly: narrowed the "before" window to 2020-2023 (post-convergence, pre-release) instead of full history, recomputed the 2x2.

**DATA / RESULTS:**  
| Metric | Value |
|---|---|
| Unique UniProt accessions needed | 15,144 |
| AlphaFold structures downloaded | 14,930 (98.6%) |
| No AlphaFold model available | 214 |
| Variants with real pLDDT annotation | 183,869 / 212,903 (86.4%) |
| Disorder class: ordered | 60.6% |
| Disorder class: disordered (IDR) | 31.1% (expected ~30% - validates the annotation) |
| Disorder class: intermediate | 8.3% |
| AUC ordered (real pLDDT) | 0.9541 |
| AUC IDR (real pLDDT) | 0.9444 |
| Gap ordered vs IDR | 0.0096 |
| AUC before AlphaMissense release (full history) | 0.9608 |
| AUC after AlphaMissense release (full history) | 0.9630 |
| Ordered x Before (full history) | 0.9538 |
| Ordered x After (full history) | 0.9531 |
| IDR x Before (full history) | 0.9333 |
| IDR x After (full history) | 0.9545 |
| IDR jump, full history | +0.0212 |
| Ordered x Before (narrow window 2020-2023) | 0.9606 |
| Ordered x After (narrow window) | 0.9531 |
| IDR x Before (narrow window 2020-2023) | 0.9504 |
| IDR x After (narrow window) | 0.9545 |
| IDR jump, narrow window | +0.0042 |

**OBSERVATIONS:**  
The pLDDT annotation pipeline is validated: 31.1% of annotated variants fall in disordered regions, closely matching the expected ~30% proteome-wide disordered fraction (the old proxy only hit 5.6% and was invalid). This replaces the session 1 proxy result entirely.

The full-history 2x2 breakdown initially looked like evidence for H1 (IDR AUC jumped +0.0212 after the AlphaMissense release while ordered stayed flat), but the figure shows the ordered/IDR AUC lines were already converging gradually between 2013 and 2021 - years before the release date. Restricting the "before" window to 2020-2023 (post-convergence, pre-release) to remove that drift shrank the IDR jump to +0.0042, essentially noise. This means the naive pre/post comparison does NOT currently show real evidence for H1 - the earlier-looking jump was mostly an artifact of pre-existing drift in ClinVar review quality over time, not something caused by AlphaMissense specifically.

This is not evidence against H1 - it means a naive difference-in-means test is the wrong tool, which is exactly why the project design already called for a formal regression discontinuity (local linear regression fit on each side of the cutoff, bandwidth-limited near the cutoff, testing for a jump in the fitted intercept) rather than a simple pre/post split. That is Phase 2 of weekly_plan.md (weeks 5-9), not yet built.

**PROBLEMS ENCOUNTERED:**  
1. MobiDB API confirmed permanently broken (returns correct headers, zero rows) - do not retry in future sessions.
2. `clinvar_am_joined.csv` from session 1 did not retain the `uniprot`/`aa_change` columns needed for residue-level mapping - had to rerun the join step with those columns added.
3. AlphaFold DB serves current structures at API path `AF-{accession}-F1-model_v6.cif`, not `v4` as originally assumed (v4 URLs return 404).
4. Naive full-history pre/post comparison was misleading due to a pre-existing convergence trend in the data (2013-2021) unrelated to the AlphaMissense release - caught this by visually inspecting the figure before trusting the number, then confirmed with the narrow-window check.

**NEXT SESSION GOAL:**  
Per weekly_plan.md Week 3 (Aug 25-31): download archival ClinVar releases at 6-month intervals (2018-2026) and build a first-classification-date table, needed for H3 and for confounder control.

**QUESTIONS TO RESEARCH:**  
- What is local linear regression in the context of regression discontinuity design, and how is the bandwidth around the cutoff chosen?  
- How do I formally test whether a discontinuity is statistically significant (not just visually apparent)?  
- Why would ClinVar review quality for IDR variants have been improving specifically between 2013 and 2021 - is there a known cause (e.g. ACMG guideline changes in 2015)?  

**Signed:** Sagar Raut &emsp; **Date:** August 12, 2026

---

