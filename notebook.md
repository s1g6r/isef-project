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
7. Wrote `scripts/disorder_split.py` → split variants by disorder proxy (am_class == ambiguous), ran temporal split at AlphaMissense release date (Sept 19, 2023), produced first figure `outputs/auc_over_time.png`  
8. Attempted MobiDB disorder annotation download via API → returned 0 rows, API endpoint appears broken. Will switch to AlphaFold pLDDT approach next session.  

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
| AlphaMissense overall AUC | **0.9592** |
| Published AUC (literature) | ~0.90-0.93 |
| Joined variants with valid LastEvaluated date | 208,824 |
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
AUC of 0.9592 is ~0.03 above the published benchmark of ~0.90-0.93. Published numbers were measured on held-out sets designed to reduce label leakage. My measurement is on all of ClinVar with no such controls. This excess is a preliminary signal consistent with circular evidence inflating apparent performance. **Recorded before any novel analysis was run.**

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

═══════════════════════════════════════════  
**DATE:** August 13, 2026  
**SESSION:** #3  
**TIME:** 9:18 AM - 10:03 AM  
═══════════════════════════════════════════  

**GOAL FOR TODAY:**  
Download archival ClinVar releases at 6-month intervals and build a first-classification-date table, per Week 3 of weekly_plan.md.

**BACKGROUND / REASONING:**  
Every analysis so far has used `LastEvaluated`, the date of a variant's most recent review. That's a confound risk: a variant can get its `LastEvaluated` bumped by an unrelated re-review long after it was actually first classified, which would make it look "post-AlphaMissense" in a naive temporal split even though its original classification predates AlphaMissense entirely. A first-classification-date table (earliest date a variant appears with a P/LP/B/LB call, not its most recent review) is needed to control for this later, and is also required for H3.

**WHAT I DID:**  
1. Checked NCBI's ClinVar FTP archive structure - found 2018-2024 releases live under `archive/{year}/`, 2025+ live flat under `archive/`, with monthly resolution available going back to 2018.
2. Wrote `scripts/get_archival_clinvar.py` - downloads 18 archival snapshots at 6-month intervals (2018-01 through 2026-07), resumable, concurrent, stdlib-only. Ran it: 18/18 succeeded, ~2.5GB total.
3. Wrote `scripts/build_first_classification_date.py` - walks all snapshots chronologically (plus the existing `variant_summary_2019-01.txt.gz` and the current full release), records the earliest snapshot each (chrom,pos,ref,alt) variant appears in as GRCh38 P/LP/B/LB.
4. First run failed: `Usecols do not match columns` for the 2018-01 through 2020-07 files. Diagnosed by comparing column headers across files directly - found NCBI only added `PositionVCF`/`ReferenceAlleleVCF`/`AlternateAlleleVCF` starting with the 2021-01 release; older files only have `Start`/`ReferenceAllele`/`AlternateAllele`.
5. Fixed `build_first_classification_date.py` to auto-detect which column set a given file has (VCF-style vs legacy) and use the matching one - the legacy columns aren't VCF-normalized but are equivalent for single-nucleotide substitutions, which is all AlphaMissense/missense variants are anyway.
6. Reran successfully - produced `data/first_classification_date.csv`.

**DATA / RESULTS:**  
| Metric | Value |
|---|---|
| Archival snapshots downloaded | 18/18 |
| Total size on disk | ~2.5 GB |
| Total unique variants with a first-seen release | 1,886,425 |
| First-seen, 2018-01 | 178,089 |
| First-seen, 2020-01 | 170,556 |
| First-seen, 2022-07 | 189,043 |
| First-seen, 2023-07 | 203,708 |
| First-seen, 2024-07 | 311,575 |
| First-seen, current (most recent) | 10,623 |

**OBSERVATIONS:**  
The total (1,886,425) is larger than the current file's filtered count (1,714,672) - expected, since this is a union across 8 years and includes variants that were P/LP/B/LB at some point but have since been reclassified to VUS or removed. Every July snapshot shows substantially more newly-first-seen variants than the following January (e.g. 2024-07: 311,575 vs 2024-01: 47,028), consistently across both the legacy-column era and the VCF-column era - this looks like a real submission-cadence pattern (large mid-year batches from clinical labs) rather than a parsing artifact, since it isn't tied to the column-format switch.

This table is not yet joined into the main analysis dataset - building it was the whole Week 3 goal. The join (on chrom/pos/ref/alt) happens when confounder control is actually needed, in Phase 2.

**PROBLEMS ENCOUNTERED:**  
1. `Usecols do not match columns` on 2018-2020 archival files - NCBI changed the ClinVar schema over time, adding VCF-normalized position/allele columns starting 2021-01. Fixed by auto-detecting column set per file rather than assuming one fixed schema.

**NEXT SESSION GOAL:**  
Per weekly_plan.md Week 4 (Sep 1-7): add ESM-1b scores (Meta, Aug 2021) to the joined dataset - predictor #2 of the 4-predictor RDD design. Need to find a precomputed genome-wide ESM-1b variant-effect score source (similar to how AlphaMissense_hg38.tsv.gz was used) rather than running the model locally.

**QUESTIONS TO RESEARCH:**  
- Is there a precomputed, downloadable genome-wide ESM-1b missense score set, or does it need to be run from the model directly?  
- Does ESM-1b use UniProt-based coordinates like AlphaMissense, or a different scheme that needs its own join logic?  

**Signed:** Sagar Raut &emsp; **Date:** August 13, 2026

---

═══════════════════════════════════════════  
**DATE:** August 13, 2026  
**SESSION:** #4  
**TIME:** 10:12 AM - 11:23 AM  
═══════════════════════════════════════════  

**GOAL FOR TODAY:**  
Add ESM-1b scores (predictor #2 of the 4-predictor design) to the joined dataset, then close out Phase 1 by merging in the first-classification-date table built in session 3.

**BACKGROUND / REASONING:**  
weekly_plan.md Week 4 calls for ESM-1b scores, confirming the join works the same way AlphaMissense's did. Separately, Phase 1's exit check (weekly_plan.md line 19) requires one clean dataset with all five pieces - ClinVar labels, AlphaMissense scores, ESM-1b scores, pLDDT disorder annotations, and first-classification dates - before moving to Phase 2. Session 3 built the first-classification-date table but hadn't merged it in yet.

**WHAT I DID:**  
1. Researched ESM-1b (Brandes et al. 2023, Nature Genetics) distribution options - found there is no single bulk genome-wide flat file like AlphaMissense's; the official distribution is an interactive web portal (Hugging Face Space).
2. Found that the portal's own Hugging Face Space repo hosts the complete precomputed catalog as a public, no-registration file: `ALL_hum_isoforms_ESM1b_LLR.zip` (1.34 GB). Confirmed this by reading the portal's own `app.py` source to understand the exact format: one CSV per UniProt accession, structured as a position x amino-acid LLR score matrix.
3. Also looked at dbNSFP (which bundles ESM-1b, REVEL, and PolyPhen-2 together in genomic coordinates) as an alternative - ruled it out for today because it requires institutional-email registration and is a ~50-80GB download. Noting it as a likely future source for REVEL/PolyPhen-2.
4. Wrote `scripts/get_esm1b_scores.py` (downloads the zip) and `scripts/annotate_esm1b.py` (maps scores onto the existing `uniprot`/`aa_change` columns - same per-protein-parsed-once pattern as session 2's `annotate_plddt.py`). The existing `aa_change` format from the AlphaMissense join happened to match ESM-1b's variant-key format exactly, so no coordinate translation was needed.
5. Ran the download and annotation: 212,676/212,903 variants scored (99.9%), only 21/15,144 proteins missing from the catalog.
6. Computed a sanity-check AUC - had to flip the sign first, since ESM-1b's LLR score uses lower (more negative) = more pathogenic, the opposite convention from am_score.
7. Wrote `scripts/merge_first_classification_date.py` to join session 3's first-classification-date table into the main dataset on chrom/pos/ref/alt, and ran it.

**DATA / RESULTS:**  
| Metric | Value |
|---|---|
| ESM-1b source | Hugging Face Space public file (ALL_hum_isoforms_ESM1b_LLR.zip, 1.34 GB) |
| Proteins in ESM-1b catalog | 42,286 |
| Proteins needed | 15,144 |
| Proteins matched | 15,123 (99.86%) |
| Variants with ESM-1b score | 212,676 / 212,903 (99.9%) |
| ESM-1b AUC vs ClinVar (sign-flipped) | 0.9299 |
| Published ESM-1b AUC (Brandes et al. 2023) | 0.885 |
| Variants matched to first-classification date | 212,903 / 212,903 (100%) |
| Phase 1 exit check - ClinVar labels | 100% |
| Phase 1 exit check - AlphaMissense scores | 100% |
| Phase 1 exit check - ESM-1b scores | 99.9% |
| Phase 1 exit check - pLDDT disorder | 86.4% |
| Phase 1 exit check - first-classification dates | 100% |

**OBSERVATIONS:**  
ESM-1b's AUC exceeds its published benchmark by +0.045 (0.9299 vs 0.885), closely paralleling AlphaMissense's own excess (0.9592 vs ~0.90-0.93 published). These are two independent predictors, trained completely differently (AlphaMissense uses structure+evolution-based deep learning, ESM-1b is a protein language model with no explicit structural input), both showing elevated performance specifically against ClinVar. That's a meaningfully stronger circularity signal than either predictor alone would be - harder to explain away as one model's individual quirk. Phase 1 is now formally complete per weekly_plan.md's exit criteria (all five pieces present in one dataset, `data/clinvar_phase1_complete.csv`).

**PROBLEMS ENCOUNTERED:**  
None major this session - the coordinate-format alignment between AlphaMissense's `aa_change` and ESM-1b's variant-key convention worked without needing a translation step, though this was confirmed by checking rather than assumed, and won't necessarily hold for REVEL/PolyPhen-2 later.

**NEXT SESSION GOAL:**  
Phase 1 is complete (all 4 weeks done). Per weekly_plan.md, Week 5 (Sep 8-14) begins Phase 2: run a naive pre/post split for AlphaMissense (Sept 2023) and ESM-1b (Aug 2021) together - the first look at temporal signal across two predictors instead of one.

**QUESTIONS TO RESEARCH:**  
- Does ESM-1b's own temporal split (before/after Aug 2021) show the same pre-existing convergence-trend confound found for AlphaMissense in session 2, or something different?  
- Where can REVEL and PolyPhen-2 genome-wide precomputed scores be obtained without dbNSFP's registration friction?  

**Signed:** Sagar Raut &emsp; **Date:** August 13, 2026

---

═══════════════════════════════════════════  
**DATE:** August 14, 2026  
**SESSION:** #5  
**TIME:** 10:18 AM - 11:10 AM  
═══════════════════════════════════════════  

**GOAL FOR TODAY:**  
Run a naive pre/post split for AlphaMissense and ESM-1b together (Week 5 of weekly_plan.md), then extend the same ordered-vs-IDR rigor session 2 applied to AlphaMissense to ESM-1b as well, checking for the same pre-existing convergence-trend confound before trusting any jump.

**BACKGROUND / REASONING:**  
Phase 1 finished ahead of schedule (sessions 3-4, same day), so Phase 2 could start immediately. Week 5's task is a first, naive look at whether either predictor's accuracy jumps at its own release date. Session 2 already showed that a naive full-history pre/post comparison can be fooled by a slow pre-existing convergence trend unrelated to the release date - that check needed to be run for ESM-1b too, not just assumed to be clean.

**WHAT I DID:**  
1. Wrote `scripts/naive_pre_post_split.py` - naive pre/post AUC split for AlphaMissense (cutoff Sept 19, 2023) and ESM-1b (cutoff Aug 1, 2021) together on `data/clinvar_phase1_complete.csv`, with ESM-1b's LLR sign flipped up front so higher = more pathogenic for both predictors. Also plots full-history year-by-year AUC for each. Ran it.
2. Wrote `scripts/disorder_split_esm1b.py` - mirrors session 2's `disorder_split_plddt.py` exactly, but for ESM-1b: ordered-vs-IDR AUC, 2x2 breakdown at the Aug 2021 cutoff, year-by-year ordered/IDR plot. Ran it.
3. Inspected `outputs/auc_over_time_esm1b_disorder.png` directly before trusting the 2x2 numbers - found the same warning sign as session 2: ordered and IDR AUC lines converge gradually from about 2012 through 2021, well before the Aug 2021 release.
4. Wrote `scripts/temporal_window_check_esm1b.py` - mirrors session 2's `temporal_window_check.py`: narrowed the "before" window to 2019-01-01 through the release date (post-most-of-the-convergence, pre-release) and recomputed the 2x2. Ran it.

**DATA / RESULTS:**  
| Metric | Value |
|---|---|
| AlphaMissense: valid rows | 208,824 / 212,903 |
| AlphaMissense: AUC before Sept 2023 | 0.9586 |
| AlphaMissense: AUC after Sept 2023 | 0.9608 |
| AlphaMissense: naive overall jump | +0.0022 |
| ESM-1b: valid rows | 208,597 / 212,903 |
| ESM-1b: AUC before Aug 2021 | 0.9110 |
| ESM-1b: AUC after Aug 2021 | 0.9351 |
| ESM-1b: naive overall jump | +0.0241 |
| ESM-1b: ordered variants / IDR variants (disorder split) | 108,437 / 56,527 |
| ESM-1b: overall AUC ordered / IDR | 0.9154 / 0.9106 |
| ESM-1b: 2x2 full history - Ordered x Before / After | 0.8896 / 0.9196 |
| ESM-1b: 2x2 full history - IDR x Before / After | 0.8681 / 0.9226 |
| ESM-1b: full-history jump - Ordered / IDR | +0.0300 / +0.0545 |
| ESM-1b: 2x2 narrow window (2019-01 to 2021-08) - Ordered x Before / After | 0.8910 / 0.9196 |
| ESM-1b: 2x2 narrow window - IDR x Before / After | 0.8938 / 0.9226 |
| ESM-1b: narrow-window jump - Ordered / IDR | +0.0286 / +0.0288 |
| For comparison, AlphaMissense (session 2) - full-history IDR jump | +0.0212 |
| For comparison, AlphaMissense (session 2) - narrow-window IDR jump | +0.0042 |

**OBSERVATIONS:**  
AlphaMissense's naive overall pre/post jump is essentially zero (+0.0022) - consistent with session 2's finding that its full-history IDR-specific jump was a convergence-trend artifact, not a real discontinuity. This session's overall split confirms the same picture from a different angle: nothing moves much around the AlphaMissense release date once you look at the whole population.

ESM-1b is different. Its naive overall jump (+0.0241) is real - and unlike AlphaMissense, it survives the narrow-window check: the disorder-stratified IDR jump shrinks from +0.0545 (full history) to +0.0288 (narrow window), but does not disappear. The key finding is that the narrow-window IDR jump (+0.0288) lands almost exactly on top of the narrow-window ordered jump (+0.0286) - meaning ordered and disordered regions are moving together by about the same amount. The full-history 2x2's apparent "IDR jumps more than ordered" pattern was mostly the pre-existing 2012-2021 convergence trend, same as AlphaMissense's case - but underneath that artifact, ESM-1b still has a genuine, uniform ~2.9% AUC increase around its own release date that AlphaMissense does not have.

This means the disorder-specific sub-hypothesis (H1: circularity worse in IDRs) is not currently supported by either predictor's naive 2x2 once the convergence confound is removed. But ESM-1b's overall (non-disorder-specific) jump is a real, still-unexplained signal worth carrying into the formal RDD in Week 6 - it just isn't the IDR-concentrated effect the raw numbers first suggested.

**PROBLEMS ENCOUNTERED:**  
1. Nearly reported ESM-1b's full-history IDR jump (+0.0545) as evidence for H1 without first checking for the same convergence-trend confound session 2 found for AlphaMissense - caught this by inspecting `outputs/auc_over_time_esm1b_disorder.png` before writing up the 2x2 numbers, same discipline as session 2.
2. Unlike AlphaMissense's case, the confound check did not fully explain away ESM-1b's jump - it only removed the part that looked disorder-specific. Needed to slow down and compare narrow-window ordered vs. IDR jumps directly (0.0286 vs 0.0288) rather than assuming "narrowed the window" automatically meant "signal debunked."

**NEXT SESSION GOAL:**  
Per weekly_plan.md Week 6 (Sep 15-21): implement the actual regression discontinuity design (local linear regression fit on each side of the cutoff, bandwidth-limited near the cutoff, testing for a jump in the fitted intercept) for both AlphaMissense and ESM-1b - the real hypothesis test the naive splits in sessions 2 and 5 were building toward.

**QUESTIONS TO RESEARCH:**  
- Why would ESM-1b show a real overall jump at its release while AlphaMissense shows none at all - is that a property of how each model was trained/released, or something about ClinVar's review practices specifically around 2021?  
- For the Week 6 RDD, what bandwidth should be used around each cutoff, and does using `first_seen_release` (6-month resolution) instead of `LastEvaluated` (continuous) as the running variable change the answer?

**Signed:** Sagar Raut &emsp; **Date:** August 14, 2026

---

