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
Session 1 ended with the am_class=='ambiguous' IDR proxy only capturing 5.6% of variants (expected ~30%), meaning it was not actually measuring structural disorder. I needed a real source of structural disorder annotations before the ordered-vs-IDR comparison could mean anything - AlphaFold's per-residue pLDDT confidence score was the obvious candidate.

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

The full-history 2x2 breakdown initially looked like evidence for what I'm calling H1 - my core hypothesis, that a predictor's ClinVar agreement should jump right at its own release date, and more so in disordered protein regions where reviewers have less other evidence to lean on. (IDR AUC jumped +0.0212 after the AlphaMissense release while ordered stayed flat.) But the figure shows the ordered/IDR AUC lines were already converging gradually between 2013 and 2021 - years before the release date. Restricting the "before" window to 2020-2023 (post-convergence, pre-release) to remove that drift shrank the IDR jump to +0.0042, essentially noise. So the naive pre/post comparison does NOT currently show real evidence for H1 - the earlier-looking jump was mostly an artifact of pre-existing drift in ClinVar review quality over time, not something caused by AlphaMissense specifically.

This isn't evidence against H1 either - it just means a naive difference-in-means test is the wrong tool for this question. What I actually need is a formal regression discontinuity design: local linear regression fit separately on each side of the cutoff, bandwidth-limited near the cutoff, testing for a jump in the fitted intercept, instead of a simple pre/post split. I haven't built that yet - that's the real next phase of this project.

**PROBLEMS ENCOUNTERED:**  
1. MobiDB API confirmed permanently broken (returns correct headers, zero rows) - do not retry in future sessions.
2. `clinvar_am_joined.csv` from session 1 did not retain the `uniprot`/`aa_change` columns needed for residue-level mapping - had to rerun the join step with those columns added.
3. AlphaFold DB serves current structures at API path `AF-{accession}-F1-model_v6.cif`, not `v4` as originally assumed (v4 URLs return 404).
4. Naive full-history pre/post comparison was misleading due to a pre-existing convergence trend in the data (2013-2021) unrelated to the AlphaMissense release - caught this by visually inspecting the figure before trusting the number, then confirmed with the narrow-window check.

**NEXT SESSION GOAL:**  
Download archival ClinVar releases at 6-month intervals (2018-2026) and build a first-classification-date table. I need this for confounder control, and also for a bonus test later I'm calling H3: whether the release-date jump disappears once I restrict to variants that have documented non-computational evidence in ClinVar, instead of just AI-predictor scores.

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
Download archival ClinVar releases at 6-month intervals and build a first-classification-date table.

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

This table isn't joined into the main analysis dataset yet - just building it was the goal for today. The join (on chrom/pos/ref/alt) can happen later, whenever confounder control is actually needed.

**PROBLEMS ENCOUNTERED:**  
1. `Usecols do not match columns` on 2018-2020 archival files - NCBI changed the ClinVar schema over time, adding VCF-normalized position/allele columns starting 2021-01. Fixed by auto-detecting column set per file rather than assuming one fixed schema.

**NEXT SESSION GOAL:**  
Add ESM-1b scores (Meta, Aug 2021) to the joined dataset - the second of four predictors I'm testing. Need to find a precomputed genome-wide ESM-1b variant-effect score source (similar to how AlphaMissense_hg38.tsv.gz was used) rather than running the model locally.

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
Add ESM-1b scores (the second of four predictors) to the joined dataset, then finish merging in the first-classification-date table built in session 3, so the core dataset finally has everything it needs.

**BACKGROUND / REASONING:**  
The next predictor to add is ESM-1b, and I want to confirm the join works the same way AlphaMissense's did. I also want one clean dataset with all five pieces I need before I can start the real analysis - ClinVar labels, AlphaMissense scores, ESM-1b scores, pLDDT disorder annotations, and first-classification dates. Session 3 built the first-classification-date table but I hadn't merged it in yet.

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
| ClinVar labels populated | 100% |
| AlphaMissense scores populated | 100% |
| ESM-1b scores populated | 99.9% |
| pLDDT disorder annotations populated | 86.4% |
| First-classification dates populated | 100% |

**OBSERVATIONS:**  
ESM-1b's AUC exceeds its published benchmark by +0.045 (0.9299 vs 0.885), closely paralleling AlphaMissense's own excess (0.9592 vs ~0.90-0.93 published). These are two independent predictors, trained completely differently (AlphaMissense uses structure+evolution-based deep learning, ESM-1b is a protein language model with no explicit structural input), both showing elevated performance specifically against ClinVar. That's a meaningfully stronger circularity signal than either predictor alone would be - harder to explain away as one model's individual quirk. The core dataset now has everything I need in one place (`data/clinvar_phase1_complete.csv`) - all five pieces present, so I can move on to the actual temporal analysis.

**PROBLEMS ENCOUNTERED:**  
None major this session - the coordinate-format alignment between AlphaMissense's `aa_change` and ESM-1b's variant-key convention worked without needing a translation step, though this was confirmed by checking rather than assumed, and won't necessarily hold for REVEL/PolyPhen-2 later.

**NEXT SESSION GOAL:**  
The core dataset is finally complete, so it's time to start the actual temporal analysis. First step: run a naive pre/post split for AlphaMissense (Sept 2023) and ESM-1b (Aug 2021) together - a first look at the temporal signal across two predictors instead of one.

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
Run a naive pre/post split for AlphaMissense and ESM-1b together, then extend the same ordered-vs-IDR rigor session 2 applied to AlphaMissense to ESM-1b as well, checking for the same pre-existing convergence-trend confound before trusting any jump.

**BACKGROUND / REASONING:**  
The core dataset finished ahead of schedule (sessions 3-4, same day), so I could start the actual temporal analysis right away. This is a first, naive look at whether either predictor's accuracy jumps at its own release date. Session 2 already showed that a naive full-history pre/post comparison can be fooled by a slow pre-existing convergence trend unrelated to the release date - I needed to check for that with ESM-1b too, not just assume it was clean because AlphaMissense wasn't.

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

This means H1 - circularity being worse specifically in IDRs - isn't currently supported by either predictor's naive 2x2, once the convergence confound is removed. But ESM-1b's overall jump, which isn't disorder-specific, is real and still unexplained. That's worth carrying into the actual RDD next, even though it isn't the IDR-concentrated effect the raw numbers first suggested.

**PROBLEMS ENCOUNTERED:**  
1. Nearly reported ESM-1b's full-history IDR jump (+0.0545) as evidence for H1 without first checking for the same convergence-trend confound session 2 found for AlphaMissense - caught this by inspecting `outputs/auc_over_time_esm1b_disorder.png` before writing up the 2x2 numbers, same discipline as session 2.
2. Unlike AlphaMissense's case, the confound check did not fully explain away ESM-1b's jump - it only removed the part that looked disorder-specific. Needed to slow down and compare narrow-window ordered vs. IDR jumps directly (0.0286 vs 0.0288) rather than assuming "narrowed the window" automatically meant "signal debunked."

**NEXT SESSION GOAL:**  
Implement the actual regression discontinuity design (local linear regression fit on each side of the cutoff, bandwidth-limited near the cutoff, testing for a jump in the fitted intercept) for both AlphaMissense and ESM-1b. This is the real test the naive splits in sessions 2 and 5 were only ever a warm-up for.

**QUESTIONS TO RESEARCH:**  
- Why would ESM-1b show a real overall jump at its release while AlphaMissense shows none at all - is that a property of how each model was trained/released, or something about ClinVar's review practices specifically around 2021?  
- For the RDD, what bandwidth should be used around each cutoff, and does using `first_seen_release` (6-month resolution) instead of `LastEvaluated` (continuous) as the running variable change the answer?

**Signed:** Sagar Raut &emsp; **Date:** August 14, 2026

---

═══════════════════════════════════════════  
**DATE:** August 15, 2026  
**SESSION:** #6  
**TIME:** 9:12 AM - 10:03 AM  
═══════════════════════════════════════════  

**GOAL FOR TODAY:**  
Implement the actual regression discontinuity design (RDD) for AlphaMissense and ESM-1b - the real hypothesis test the naive splits in sessions 2 and 5 were only ever leading up to.

**BACKGROUND / REASONING:**  
Sessions 2 and 5 both showed that comparing flat before/after averages can be misleading - a slow pre-existing trend can fake a jump, or hide a real one, depending on where you draw the "before" window. What's actually needed is a method that only looks at data close to the release date and fits a separate trend on each side, so a real discontinuity can be told apart from an ongoing drift. That's what RDD is for.

**WHAT I DID:**  
1. Wrote `scripts/rdd_analysis.py`. Since AUC is an aggregate across many variants and isn't something you can regress a single row on, I used a per-variant outcome instead: whether the predictor's own published call matched the ClinVar label (1/0) - AlphaMissense's own `am_class` categories (excluding "ambiguous," same as excluding "intermediate" pLDDT elsewhere), and ESM-1b's own published -7.5 LLR threshold. Fit one OLS regression with a before/after dummy and its interaction with time, on data windowed around each release date - the coefficient on the dummy is the jump estimate, with heteroskedasticity-robust standard errors computed by hand (no new packages needed).
2. Ran it at three different bandwidths (1, 2, and 3 years) for each predictor instead of picking one window and trusting it, since RDD estimates are known to be sensitive to bandwidth choice.
3. Hit a wall of RuntimeWarnings on every run ("divide by zero"/"overflow encountered in matmul"). Instead of just suppressing them, checked whether the actual numbers were wrong first - wrote the same regression by hand with no BLAS calls at all and compared it against the warned version row by row. They matched exactly. Turned out to be a known cosmetic bug in how numpy's Apple Accelerate backend handles certain matrix-vector shapes on M-series Macs, not a real numerical failure. Suppressed the specific warning once I'd confirmed it was safe to.
4. Looked at `outputs/rdd_plots.png` (monthly-binned scatter with the two fitted local-linear lines) before trusting any of the jump numbers.

**DATA / RESULTS:**  
| Metric | Value |
|---|---|
| AlphaMissense RDD, bandwidth 1yr: n, jump, p-value | n=50,928, jump=+0.0046, p=0.353 |
| AlphaMissense RDD, bandwidth 2yr: n, jump, p-value | n=112,992, jump=+0.0017, p=0.619 |
| AlphaMissense RDD, bandwidth 3yr: n, jump, p-value | n=168,250, jump=+0.0037, p=0.208 |
| ESM-1b RDD, bandwidth 1yr: n, jump, p-value | n=21,658, jump=+0.0237, p=0.0024 |
| ESM-1b RDD, bandwidth 2yr: n, jump, p-value | n=51,540, jump=-0.0026, p=0.677 |
| ESM-1b RDD, bandwidth 3yr: n, jump, p-value | n=87,154, jump=+0.0019, p=0.724 |

**OBSERVATIONS:**  
AlphaMissense shows no real discontinuity at any bandwidth - all three estimates are small and nowhere near significant (p between 0.21 and 0.62), and in the plot the fitted lines from each side meet almost seamlessly right at the release date. This is a clean, robust null result, consistent with everything found for this predictor so far.

ESM-1b turned out more complicated than the naive check in session 5 suggested. Only the narrowest window (1 year) shows a statistically significant jump (+0.0237, p=0.0024) - widen it to 2 or 3 years and the effect shrinks to nearly zero and even flips sign. The plot explains why: P(correct) climbs steadily for about two years before the release, peaks right around the release date, then declines afterward. That's a rise-then-fall shape, not a clean sustained jump, and it means this result isn't robust to bandwidth choice - one of the most common reasons an RDD finding turns out weaker than it first looks. I'm not treating the +0.0237 number as a confirmed finding. The honest read is that something was happening around ESM-1b's release, but it doesn't hold up as a stable discontinuity the way a real one should.

**PROBLEMS ENCOUNTERED:**  
1. The script threw a wall of RuntimeWarnings on every regression call. Rather than assume they were harmless or panic that the numbers were wrong, verified directly: wrote the same regression by hand with no BLAS calls and compared results row by row - they matched exactly. Confirmed this as a known benign numpy/Apple Accelerate quirk on M-series Macs before suppressing it.
2. Almost reported ESM-1b's 1-year-bandwidth jump as the finding without checking the other two bandwidths first - it only looked clean because that happened to be the window I plotted. Checking 2yr and 3yr caught this before it went in as a real result.

**NEXT SESSION GOAL:**  
Source REVEL and PolyPhen-2, the last two of the four predictors, so all four can eventually go through the same RDD. Also worth revisiting the ESM-1b rise-then-fall pattern more carefully instead of dropping it - maybe by looking into what was actually happening in ClinVar review practices around 2020-2022.

**QUESTIONS TO RESEARCH:**  
- Is there a plausible non-circularity explanation for ESM-1b's rise-then-fall pattern (e.g. a temporary spike in reliance on it right when it was novel, followed by a correction), or is this more likely just noise from a smaller "before" sample?  
- Is there a more principled, data-driven way to choose an RDD bandwidth instead of round numbers like 1/2/3 years, so the choice isn't arbitrary?

**Signed:** Sagar Raut &emsp; **Date:** August 15, 2026

---

═══════════════════════════════════════════  
**DATE:** August 15, 2026  
**SESSION:** #7  
**TIME:** 10:11 AM - 11:26 AM  
═══════════════════════════════════════════  

**GOAL FOR TODAY:**  
Source the last two predictors, REVEL and PolyPhen-2, and go back to the open question session 6 left about ESM-1b's release-date pattern instead of just leaving it flagged.

**BACKGROUND / REASONING:**  
Session 6's RDD left two things unresolved: REVEL and PolyPhen-2 still hadn't been sourced, and ESM-1b's release-date jump wasn't robust across bandwidths, with a rise-then-fall shape that deserved a real explanation rather than just getting noted as "not robust" and set aside.

**WHAT I DID:**  
1. Looked into where REVEL and PolyPhen-2 actually come from before assuming dbNSFP's registration-gated site was the only option. Found REVEL has its own clean direct download (526MB, no registration). Found PolyPhen-2's own bulk file (WHRESS) too, but it uses RefSeq protein-relative coordinates instead of genomic coordinates, which would need its own translation step - the same kind of friction I'd already avoided once before with ESM-1b.
2. Found a better option: dbNSFP has a full academic version (v4.1a) mirrored on Zenodo that doesn't require registration. There's also a "commercial-safe" version that quietly strips out Polyphen2, REVEL, and CADD for licensing reasons - checked the record description carefully to make sure I had the right one. Its coordinates are already chr/pos/ref/alt in hg38, the same join key used for every other predictor in this project, so one download solves both remaining predictors at once.
3. Wrote `scripts/get_dbnsfp.py` and started the download (~25GB across 24 per-chromosome files). Still running as of this entry.
4. Wrote `scripts/annotate_dbnsfp.py` to join REVEL and Polyphen2 scores in once the download finishes - not run yet.
5. While that downloaded, went back to session 6's open question. Checked whether AlphaMissense's retrospective accuracy shows the same rise-then-fall shape over the same calendar months as ESM-1b - if it does, that would prove the pattern isn't specific to ESM-1b at all, since AlphaMissense didn't exist until 2023.
6. Found that both predictors peak in the exact same month: August 2021, with a 0.60 correlation across 5 years of overlapping months. AlphaMissense's actual release date (Sept 2023) has nothing to do with August 2021, which rules out "ESM-1b's release caused this" as the explanation.
7. Dug into what actually changed in ClinVar that month: a 4.6x volume surge, and the label mix flipping from mostly-pathogenic to mostly-benign, spread across many genes rather than one big single-submission dump.

**DATA / RESULTS:**  
| Metric | Value |
|---|---|
| REVEL source found | own direct download, 526MB, no registration (not used - see below) |
| PolyPhen-2 source found (WHRESS) | 8.2GB, RefSeq protein coordinates, not used |
| dbNSFP version/source chosen | v4.1a via Zenodo, no registration, ~25.2GB, covers REVEL + Polyphen2 both |
| dbNSFP download status | in progress as of this entry |
| AlphaMissense peak month (calendar 2019-2023) | August 2021, P(correct)=0.9682 |
| ESM-1b peak month (calendar 2019-2023) | August 2021, P(correct)=0.9357 |
| Correlation between AM and ESM-1b monthly correctness | 0.601 (60 overlapping months) |
| Row count, Jul-Aug 2021 vs Jul-Aug 2020 | 2,643 vs 573 (4.6x) |
| Label balance, Jul-Aug 2021 (pathogenic / benign) | 26.5% / 73.5% |
| Label balance, Jul-Aug 2020 (pathogenic / benign) | 70.0% / 30.0% |
| Top single gene's share of the Jul-Aug 2021 batch | GATA2, 0.8% of the batch |

**OBSERVATIONS:**  
The shared August 2021 peak across two unrelated predictors - one of which didn't exist yet - rules out circularity around ESM-1b's own release as the explanation for session 6's jump. What's actually happening looks like a data-composition confound: ClinVar received an unusually large, mostly-benign batch of new entries in July-August 2021, and a benign-heavy batch is easier for almost any reasonable predictor to score correctly, since rejecting a clear benign variant is generally an easier call than confirming a pathogenic one. That would lift both predictors' apparent accuracy that month for reasons that have nothing to do with either tool specifically. This also connects back to session 3's finding that ClinVar gets large, predictable batch submissions every July - this looks like an unusually large instance of that same known pattern, not something ESM-1b caused. The gene spread (no single gene above 1% of the batch) argues against this being one lab's single gene-panel dump.

This changes how I should read session 6's ESM-1b result: the +0.0237 jump at the 1-year bandwidth isn't just "not robust to bandwidth choice" in some generic statistical sense anymore - there's now a specific, plausible, and testable non-circularity explanation for it. That's a stronger conclusion than flagging instability and moving on.

**PROBLEMS ENCOUNTERED:**  
1. Would have defaulted to dbNSFP's official registration-gated site for REVEL/PolyPhen-2, per what was noted as the plan back in session 4. Checked for alternatives first instead - found REVEL's own registration-free download, and more importantly a full-academic Zenodo mirror of dbNSFP itself, avoiding the registration wait entirely.
2. Didn't accept session 6's "not robust to bandwidth choice" as the final word on ESM-1b - the rise-then-fall shape was specific enough to be testable, and using AlphaMissense as a natural control (since it didn't exist in 2021) turned a vague caveat into a concrete explanation.

**NEXT SESSION GOAL:**  
Once the dbNSFP download finishes, run `scripts/annotate_dbnsfp.py` to join REVEL and Polyphen2 scores onto the main dataset, closing out predictors #3 and #4. Then all four predictors are ready for the same RDD treatment session 6 built for AlphaMissense and ESM-1b.

**QUESTIONS TO RESEARCH:**  
- Is there a way to identify which submitters were behind the Jul-Aug 2021 batch, to confirm the benign-variant-batch explanation more directly instead of inferring it from the label mix alone?  
- Does the same kind of submission-batch confound show up around REVEL's (2016) and PolyPhen-2's (2010) release dates too - should every predictor's RDD get this same "test with an unrelated predictor as a control" check by default going forward?

**Signed:** Sagar Raut &emsp; **Date:** August 15, 2026

---

═══════════════════════════════════════════  
**DATE:** August 15, 2026  
**SESSION:** #8  
**TIME:** 11:31 AM - 1:24 PM  
═══════════════════════════════════════════  

**GOAL FOR TODAY:**  
Finish joining REVEL and Polyphen2 scores once the dbNSFP download completed, then extend the RDD to all four predictors so today closes with the same real discontinuity test run on everything instead of just the first two.

**BACKGROUND / REASONING:**  
Session 7 left the dbNSFP download running and the join script untested. Once the download finished, running the join was the last step before all four predictors would be in one place and ready for the RDD built in session 6.

**WHAT I DID:**  
1. Ran `scripts/annotate_dbnsfp.py` on the completed download - it crashed immediately with a memory error. The dbNSFP files enumerate every possible amino acid substitution genome-wide (tens of millions of rows per chromosome), and I'd read each one fully into memory as string data before filtering - that blew up to 50+GB of RAM and froze the machine. Had to force-kill the process.
2. Rewrote the join to stream each file in chunks and filter down to only the ~213K variants already in the dataset as it goes, instead of loading everything first. Also fixed a column-name bug caught on the first run (`chr` vs the actual `#chr` header) - checked the real file header directly instead of guessing again.
3. Reran the fixed version successfully - matched 209,796 of ~209,931 unique variant keys, producing `data/clinvar_phase1_and_2_complete.csv`.
4. Computed sanity-check AUCs for REVEL and both Polyphen2 variants against their published benchmarks, same pattern as every previous predictor.
5. Extended `scripts/rdd_analysis.py` to all four predictors (REVEL cutoff Oct 2016, Polyphen-2 cutoff Feb 2010), using each score's own 0.5 threshold as the predicted-pathogenic call, since both REVEL and Polyphen2 are explicitly designed as pathogenicity probabilities rather than having a published three-way class like AlphaMissense.
6. Ran it and looked at the resulting 2x2 plot before trusting any of the numbers - the REVEL and Polyphen-2 panels looked visibly different from AlphaMissense and ESM-1b's.

**DATA / RESULTS:**  
| Metric | Value |
|---|---|
| Unique variant keys matched in dbNSFP | 209,796 / 209,931 |
| REVEL_score populated | 211,325 / 212,903 (99.3%) |
| Polyphen2_HVAR_score populated | 197,133 / 212,903 (92.6%) |
| Polyphen2_HDIV_score populated | 197,133 / 212,903 (92.6%) |
| REVEL AUC vs ClinVar | 0.9670 (published range: 0.90-0.96 - above it) |
| Polyphen2 (HVAR) AUC vs ClinVar | 0.9169 (published range: ~0.80-0.88 - well above it) |
| Polyphen2 (HDIV) AUC vs ClinVar | 0.8973 (published range: ~0.85-0.92 - inside it) |
| REVEL RDD jump, bandwidth 1yr / 2yr / 3yr | -0.0137 (p=0.58) / -0.0271 (p=0.084) / +0.0517 (p=0.00005) |
| PolyPhen-2 RDD jump, bandwidth 1yr / 2yr / 3yr | +0.0114 (p=0.89) / +0.0082 (p=0.88) / +0.0419 (p=0.35) |
| PolyPhen-2 RDD sample size, bandwidth 1yr / 2yr / 3yr | n=359 / n=782 / n=1,385 |

**OBSERVATIONS:**  
REVEL and Polyphen2 (HVAR) both join AlphaMissense and ESM-1b in showing AUC excess over their own published benchmarks - REVEL at 0.9670 against a published 0.90-0.96, and Polyphen2-HVAR at 0.9169 against a published ~0.80-0.88, a substantial gap. That's four architecturally unrelated predictors (structure+evolution deep learning, a protein language model, a random-forest ensemble, and a much older single classifier) all showing the same pattern - stronger evidence than any one of them alone. Polyphen2-HDIV is the one exception, falling inside its published range, which is a real and specific nuance worth keeping rather than glossing over (HDIV is tuned for a different classification task than HVAR).

The RDD results for REVEL and Polyphen-2 tell a very different, and just as important, story. REVEL's jump estimate flips sign depending on bandwidth (negative at 1yr and 2yr, only turning positive and "significant" at 3yr) - a jump that isn't stable in sign as the window changes is not something to trust, the same lesson session 6 already taught with ESM-1b's bandwidth sensitivity, just more extreme here. Polyphen-2 is worse: none of the three bandwidths come close to significance, and the sample sizes are tiny (as low as 359 total variants). The reason is structural, not a mistake: Polyphen-2 was published in 2010, and ClinVar itself didn't launch until 2013 - there's barely any real "before" data to compare against, since most of the comparison window predates ClinVar's own existence. That's not something a bigger bandwidth or a different threshold would fix.

**PROBLEMS ENCOUNTERED:**  
1. First run of `annotate_dbnsfp.py` used 50+GB of RAM and froze the machine, requiring a force-kill (`kill -9`). Root cause: loaded every row of every dbNSFP chromosome file into memory before filtering, when only ~213K of the tens of millions of rows per file were actually needed. Fixed by streaming each file in chunks and filtering as it reads, which also made the whole run much faster.
2. Guessed a column name (`chr`) from documentation instead of checking the actual file - the real header was `#chr`. Caught it from the error message and verified the real header directly in the downloaded file before fixing it, rather than guessing again.
3. Several dbNSFP download interruptions this session from switching networks mid-download (home fiber going out and back, hotspot switches for testing) - each one silently stalled the in-flight files rather than erroring out cleanly, since `urllib.request.urlretrieve` doesn't detect a dead connection on its own. Had to check file sizes over time to confirm a stall rather than trust that the process was still working.
4. REVEL's RDD estimate is not stable across bandwidths, and Polyphen-2's RDD can't really be run at all given how little ClinVar data predates its 2010 release. Both are being reported as inconclusive/not usable rather than stretched into a finding either way.

**NEXT SESSION GOAL:**  
Run the permutation test (randomizing the cutoff date and re-running the RDD many times) for AlphaMissense and ESM-1b, the two predictors where the RDD is actually usable - this is the rigorous significance check the OLS p-values have been standing in for since session 6. Also worth deciding what, if anything, can be done for REVEL (maybe a wider bandwidth or a different running variable) versus just documenting Polyphen-2 as untestable with this design and moving on.

**QUESTIONS TO RESEARCH:**  
- Is there a standard way researchers handle an RDD when the "before" period barely has any data, or is "this predictor can't be tested this way" simply the correct conclusion?  
- For REVEL specifically, would using `first_seen_release` instead of `LastEvaluated` as the running variable change the sign-flipping behavior, since the small pre-2016 sample might be disproportionately affected by the LastEvaluated-bump confound session 3 was built to control for?

**Signed:** Sagar Raut &emsp; **Date:** August 15, 2026

---

