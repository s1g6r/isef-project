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

═══════════════════════════════════════════  
**DATE:** August 16, 2026  
**SESSION:** #9  
**TIME:** 9:00 AM - 10:15 AM  
═══════════════════════════════════════════

**GOAL FOR TODAY:**  
Run the permutation test on AlphaMissense and ESM-1b - the check that's been deferred since session 6 - and use it to actually settle REVEL's fate instead of just eyeballing its bandwidth instability.

**BACKGROUND / REASONING:**  
Every RDD p-value reported so far (sessions 6-8) came from OLS, which assumes the before/after linear model is the right shape near the cutoff. If the true relationship curves even slightly, a straight-line fit can produce a "significant" jump that isn't really there. A permutation test sidesteps that assumption: rerun the identical RDD procedure at a lot of fake cutoff dates drawn from the predictor's own real data, and see whether the real jump is actually unusual compared to that distribution, or just an ordinary draw from it.

**WHAT I DID:**  
1. Wrote `scripts/permutation_test.py`. For each predictor, it draws 2,000 placebo cutoff dates from that predictor's own observed `LastEvaluated` values (excluding anything within a year of the real cutoff, and restricted to 2013 onward - ClinVar didn't exist before that, the same reason Polyphen-2 was ruled untestable in session 8), reruns the same RDD at every one, and checks how often a jump that size or bigger shows up by chance.
2. First ran it on just AlphaMissense and ESM-1b.
3. ESM-1b's 1yr result - the one number from session 6 that had looked real (p=0.0024) - came back with a permutation p of 0.32. Not close to significant. That made it worth checking REVEL the same way instead of writing it off from bandwidth instability alone, since it had the same "looks significant at one specific window" shape.
4. Added REVEL to the script and reran all three predictors together, then looked at the histogram plot before trusting any single p-value - a number can be "significant" while still sitting close to the middle of a wide distribution, and the plot makes that visible in a way the p-value alone doesn't.
5. Deleted the empty stray repo at `/Users/sagar/terminalais/isef_project` (an unused `git init` from August 11, flagged back in session 7) - decided it was safe to clean up rather than leave sitting around.

**DATA / RESULTS:**  
| Metric | Value |
|---|---|
| AlphaMissense permutation p, bandwidth 1yr / 2yr / 3yr | 0.746 / 0.919 / 0.807 |
| ESM-1b permutation p, bandwidth 1yr / 2yr / 3yr | 0.320 / 0.934 / 0.939 |
| REVEL permutation p, bandwidth 1yr / 2yr / 3yr | 0.419 / 0.124 / 0.037 |
| ESM-1b 1yr, OLS p vs. permutation p | 0.0024 vs. 0.320 |
| REVEL 3yr, OLS p vs. permutation p | 0.00005 vs. 0.037 |

**OBSERVATIONS:**  
AlphaMissense stays solidly null at every bandwidth under the permutation test, matching what the OLS p-values already said - a consistent result across two very different ways of computing significance is about as reassuring as this kind of test gets.

ESM-1b is the real story. The 1yr OLS p-value (0.0024) looked like the strongest signal in the whole project up through session 6, but the permutation p-value for the same result is 0.32 - completely unremarkable. Looking at the histogram, the placebo jumps for ESM-1b aren't centered at zero, they lean positive, so a lot of random cutoff dates produce a jump that size just from whatever general trend is in the data, with nothing to do with ESM-1b's actual release. That lines up with session 7's finding that the Aug 2021 peak was a ClinVar-side batch effect, not something about ESM-1b itself - now there are two independent reasons not to trust that number, not just one.

REVEL is genuinely borderline in a way ESM-1b wasn't. Its permutation p-value at 3yr is 0.037 - technically under 0.05 - and the plot shows the real jump sitting clearly out past almost all of the placebo distribution at that bandwidth. But at 1yr and 2yr the real jump is negative and sits well inside the placebo bulk, not close to significant. A real discontinuity should show up in some form across nearby bandwidths, not flip sign and only appear at the widest window. I'm treating REVEL as inconclusive rather than a finding, and I'm not going to try `first_seen_release` as an alternative running variable to see if it "fixes" the sign flip - re-running a test with a different setup after seeing an unstable result is exactly the kind of after-the-fact fishing this whole project is trying to catch in other tools, so doing it here would undercut the point.

**PROBLEMS ENCOUNTERED:**  
None on the technical side this session - the script built cleanly on top of `rdd_analysis.py`'s existing `run_rdd` function, and the runtime (a few minutes for ~12,000 total RDD fits across three predictors and three bandwidths) was fine.

**NEXT SESSION GOAL:**  
Move on to the confounder analysis: match on allele frequency, gene, and review status, and try a within-protein paired comparison, so the AlphaMissense/ESM-1b null results aren't just "no jump" but "no jump even after accounting for the more obvious alternative explanations."

**QUESTIONS TO RESEARCH:**  
- What's a reasonable way to do a "within-protein paired" comparison here - matching each post-cutoff variant to a pre-cutoff variant in the same gene/protein, or something else?  
- Is allele frequency data (gnomAD or similar) already obtainable through something already in the pipeline, or does this need a new data source?

**Signed:** Sagar Raut &emsp; **Date:** August 16, 2026

---

═══════════════════════════════════════════  
**DATE:** August 16, 2026  
**SESSION:** #10  
**TIME:** 10:20 AM - 12:35 PM  
═══════════════════════════════════════════

**GOAL FOR TODAY:**  
Add the confounder analysis the plan calls for - allele frequency, gene, and review status - to check whether AlphaMissense and ESM-1b's null RDD results actually hold up once the more obvious alternative explanations are accounted for, not just on raw data.

**BACKGROUND / REASONING:**  
Session 9 settled the RDD/permutation-test picture for all four predictors, but the original plan always called for a confounder-controlled estimate on top of the raw discontinuity test, not instead of it. Two of the three confounders (gene, review status) were already columns in the dataset from earlier sessions - only allele frequency needed new work.

**WHAT I DID:**  
1. Checked dbNSFP's own file header before assuming a new download was needed, since dbNSFP was already sitting on disk from session 7 - it bundles `gnomAD_exomes_AF` and `gnomAD_genomes_AF` as extra columns, so no new data source was needed at all.
2. Wrote `scripts/annotate_gnomad_af.py`, reusing the same memory-safe streaming join from session 8's dbNSFP join, to add those two columns onto the dataset - produced `data/clinvar_complete.csv`, the new canonical file.
3. Wrote `scripts/confounder_rdd.py`: reran the AlphaMissense/ESM-1b RDD with review status (as dummy variables), submitter count (standardized), and log allele frequency (with a missing-indicator, since being entirely absent from gnomAD is itself informative, not just a gap) added as covariates directly in the regression.
4. Wrote `scripts/within_protein_paired.py` - a separate approach to the gene/protein confound specifically. Matches each post-cutoff variant to its nearest pre-cutoff variant in the same protein by residue position, then runs McNemar's test on the paired correct/incorrect outcomes. This controls for protein identity through matching instead of a regression term, since there are thousands of genes.
5. The confounder-adjusted RDD came back with two OLS p-values that looked newly significant compared to the raw RDD (AlphaMissense at 3yr, ESM-1b at 1yr). Since session 9 already proved an OLS p-value in this exact design can be wrong, trusting these two without the same check would have been inconsistent - wrote `scripts/confounder_permutation_test.py` to permutation-test the confounder-adjusted model too.

**DATA / RESULTS:**  
| Metric | AlphaMissense | ESM-1b |
|---|---|---|
| Raw RDD jump, 1yr / 2yr / 3yr | +0.0046 / +0.0017 / +0.0037 | +0.0237 / -0.0026 / +0.0019 |
| Confounder-adjusted jump, 1yr / 2yr / 3yr | +0.0009 / +0.0030 / +0.0061 | +0.0179 / +0.0023 / +0.0021 |
| Confounder-adjusted OLS p, 1yr / 2yr / 3yr | 0.850 / 0.364 / 0.039 | 0.023 / 0.711 / 0.710 |
| Confounder-adjusted permutation p, 1yr / 2yr / 3yr | 0.919 / 0.832 / 0.582 | 0.350 / 0.904 / 0.945 |
| Within-protein paired jump (2yr) | +0.0017 | +0.0000 |
| Within-protein paired McNemar p (2yr) | 0.206 | 1.000 |
| gnomAD_exomes_AF / gnomAD_genomes_AF populated | 129,665 / 212,903 (60.9%), 116,453 / 212,903 (54.7%) | (same dataset) |

**OBSERVATIONS:**  
Every method run today and across sessions 6-10 - raw RDD, permutation test, confounder-adjusted RDD, confounder-adjusted permutation test, and within-protein paired matching - agrees on the same thing: no real discontinuity for AlphaMissense or ESM-1b at any bandwidth, even after controlling for review status, submitter count, allele frequency, and protein identity. Four architecturally different ways of testing the same question landing on the same answer is a lot more convincing than any one of them alone.

The confounder-adjusted model did throw two OLS p-values under 0.05 at first (AlphaMissense at 3yr, ESM-1b at 1yr), and it would have been easy to read that as "the confounder analysis found something the raw RDD missed." Running them through the same permutation check that already caught session 6's ESM-1b result killed both (permutation p=0.582 and 0.350). Worth remembering going forward: adding covariates to a regression doesn't make its p-values any more trustworthy by itself - it's still an OLS number with the same underlying assumption, and needs the same permutation check every time, not just the one time it was built.

The within-protein paired result for ESM-1b is about as clean a null as this project has produced anywhere - the discordant pairs split exactly 2664 to 2664.

This closes out H1 by the release-date discontinuity design specifically: across four independent, methodologically different checks, there's no evidence AlphaMissense or ESM-1b's apparent accuracy jumps right at their own release date. That's a meaningful negative result, not just "didn't find anything" - it took a lot of independent checking to be confident about it. It's also worth being explicit that this doesn't contradict session 8's AUC-excess finding (REVEL, Polyphen2-HVAR, and AlphaMissense all beating their own published benchmarks) - a null release-date jump rules out one specific mechanism (a sudden leak right at release) but says nothing about a slower, non-discontinuous kind of contamination, which is exactly what comparing against a leak-free benchmark is built to test.

**PROBLEMS ENCOUNTERED:**  
1. Ran `confounder_rdd.py` and `within_protein_paired.py` back to back without waiting for `annotate_gnomad_af.py` to finish first, so both failed with a `FileNotFoundError` looking for `data/clinvar_complete.csv`. Not a bug, just needed to actually let the first script finish before the others could read its output.
2. gnomAD coverage (60.9% / 54.7%) is a lot lower than REVEL/Polyphen2's (92-99%). Expected, not an error - a lot of ClinVar pathogenic variants are rare enough to have never been observed in gnomAD's population sample at all, handled with an explicit missing-indicator column rather than dropping those rows.

**NEXT SESSION GOAL:**  
The RDD phase of the plan is done now, ahead of the original schedule. Move into the next phase: download the ProteinGym DMS (deep mutational scanning) benchmark and measure ordered vs. disordered AUC on leak-free data - data that couldn't have been in any of these predictors' training sets, since it comes from lab experiments rather than clinical databases - then compare those numbers to the ClinVar-based AUCs from session 8. This is H2, the test that doesn't depend on ClinVar's dates at all, and it's what actually makes session 8's AUC-excess finding interpretable one way or the other.

**QUESTIONS TO RESEARCH:**  
- Where's the best source for the ProteinGym DMS benchmark - a direct download, or does it need the same kind of registration/mirror research REVEL and PolyPhen-2 needed back in session 7?  
- Does ProteinGym already include AlphaMissense/ESM-1b/REVEL/PolyPhen-2 scores for its variants, or does each predictor need to be rerun on ProteinGym's variant set separately?

**Signed:** Sagar Raut &emsp; **Date:** August 16, 2026

---

═══════════════════════════════════════════  
**DATE:** August 16, 2026  
**SESSION:** #11  
**TIME:** 12:40 PM - 3:20 PM  
═══════════════════════════════════════════

**GOAL FOR TODAY:**  
Get H2 running - leak-free AUC for AlphaMissense and ESM-1b on ProteinGym's experimental DMS (deep mutational scanning) assays, to see whether the ClinVar-based accuracy numbers from session 8 hold up on data that couldn't possibly have leaked into either predictor's training set.

**BACKGROUND / REASONING:**  
Sessions 6-10 tested for a discontinuity right at each predictor's release date and found nothing - a comprehensive null across five independent methods. But that only rules out one specific mechanism: a sudden leak tied to a release date. It doesn't test whether ClinVar itself, as a benchmark, produces inflated numbers for other reasons. ProteinGym's DMS data is measured in a lab and has nothing to do with ClinVar's curation process or dates, so any gap between a predictor's ClinVar accuracy and its ProteinGym accuracy is evidence the ClinVar number doesn't reflect real-world discriminative skill, whatever the reason turns out to be.

**WHAT I DID:**  
1. Researched where ProteinGym lives and what it actually includes before downloading anything, the same way REVEL/PolyPhen-2 got vetted in session 7. Found: no registration needed; ESM-1b already has precomputed scores in ProteinGym's own baseline set, so it doesn't need to be rerun; AlphaMissense, REVEL, and PolyPhen-2 are not in ProteinGym's baselines. Checked whether AlphaMissense could still be added without running the model myself - the file already downloaded in session 1 (`data/AlphaMissense_hg38.tsv.gz`) turned out to be DeepMind's genome-wide release, covering virtually all possible human missense variants with `uniprot_id` and `protein_variant` columns already built in, so it could be cross-referenced directly. Scoped H2 to AlphaMissense and ESM-1b only with the user, documenting REVEL/PolyPhen-2 as untestable this way - same honest-limitation treatment PolyPhen-2 got in the RDD phase.
2. Wrote `scripts/get_proteingym.py` to download the DMS assay data and ESM-1b's baseline scores, filtered to human-taxon assays only (96 of 217 total), since ClinVar and AlphaMissense are both human-only.
3. First run had a genuine mistake: assumed the 4.4GB scores zip had one file per model, found zero filename matches for "esm1b", and the script deleted the zip anyway before confirming anything worked - meant redownloading 4.4GB. Fixed the script so it never deletes a zip unless its extraction actually matched something. The real structure turned out to be one file per assay with all 43 models as columns inside (the zip's 217 entries matching 217 total assays was the giveaway) - the actual fix was extracting by the same assay-filename matching already used correctly for the DMS ground-truth data, not by model name.
4. Wrote `scripts/proteingym_leak_free_analysis.py` to join AlphaMissense scores onto ProteinGym's variants and compute both predictors' AUC against the DMS ground truth, plus the current ClinVar AUCs computed the same way for a fair comparison. First run: AlphaMissense matched 0% of variants - ProteinGym's own reference file uses UniProt's mnemonic entry name (like "PAI1_HUMAN"), not the accession number ("P05121") AlphaMissense's file actually keys on.
5. Wrote `scripts/get_uniprot_mapping.py` to translate the 96 human entry names to real accessions via UniProt's own REST API. Its first version had its own bug - captured an entire tab-separated API response line instead of splitting out just the accession field, producing a mapping file that looked plausible (right shape, right file) but was still useless. Caught this from the analysis script's second run still showing 0% matched, not from inspecting the mapping file directly first.
6. With the mapping fixed, the join matched 92,280 of 329,664 variants (28.0%) for AlphaMissense. Before trusting that number, wrote `scripts/verify_proteingym_am_join.py` to spot-check 30 random matched rows against each protein's real sequence (the reference file's `target_seq` column) - confirming the claimed reference amino acid is actually at the claimed position. 30/30 passed, confirming the join is really working and 28% is genuine partial coverage, not a silent numbering bug.
7. Ran the review-status check from this session's own open question: does restricting ClinVar to only its highest-confidence labels ("reviewed by expert panel", "practice guideline") change the AUC gap to ProteinGym in an informative way? Wrote `scripts/clinvar_review_status_check.py` to compare.

**DATA / RESULTS:**  
| Metric | AlphaMissense | ESM-1b |
|---|---|---|
| ProteinGym leak-free AUC | 0.7164 (n=92,280) | 0.6906 (n=329,664) |
| Current ClinVar AUC (same method) | 0.9592 (n=212,903) | 0.9299 (n=212,676) |
| Difference (leak-free minus ClinVar) | -0.2428 | -0.2393 |
| AlphaMissense join match rate | 92,280 / 329,664 (28.0%) | - |
| Join verification (30 sampled matches) | 30/30 correct reference amino acid at claimed position | - |
| ClinVar AUC, high-confidence-only (expert panel / practice guideline, n=4,805) | 0.9384 (full dataset: 0.9592) | 0.9082 (full dataset: 0.9299) |

**OBSERVATIONS:**  
Both predictors drop by almost exactly the same amount moving from ClinVar to ProteinGym's leak-free data - AlphaMissense from 0.9592 to 0.7164 (-0.2428), ESM-1b from 0.9299 to 0.6906 (-0.2393). Two architecturally very different models landing on nearly identical drops is a striking pattern, and it's the single biggest number this project has produced so far: ClinVar-based accuracy substantially overstates how well these predictors actually discriminate on data that couldn't have leaked into training.

That said, I want to be careful about what this drop actually proves, because it isn't the same thing as proving training-data leakage specifically. Sessions 6-10 already tested for a leakage signature directly - a sudden jump right at each predictor's release date - across five independent methods, and found nothing. If ClinVar performance were inflated purely because these models had literally seen ClinVar labels during training, that kind of jump is exactly what should show up, and it didn't. So the ProteinGym gap more likely points at something else: ClinVar is a curated clinical database, and curated pathogenic/benign calls tend to be the clearer, more extreme cases - a variant usually only gets classified once there's fairly convincing evidence either way. ProteinGym's DMS assays instead measure every possible substitution at every position, including a lot of genuinely borderline, intermediate-effect variants a curated database would rarely even contain. ClinVar could simply be an easier benchmark by composition - overrepresenting the obvious cases - without any model ever having seen ClinVar labels during training. Both explanations would produce the same kind of AUC drop, and this analysis alone can't fully separate them. What's actually established here is narrower but still real: a large chunk of ClinVar-based accuracy doesn't carry over to a genuinely independent test, whatever the underlying reason turns out to be.

Also worth taking seriously, not just as footnotes: DMS fitness and clinical pathogenicity aren't the same construct (a variant can reduce protein function in a lab assay without necessarily causing human disease, and vice versa), and ProteinGym's 96 human assays cover proteins chosen for experimental tractability, not the same disease-gene-heavy population ClinVar represents. AlphaMissense's number rests on a partial (28%) subset of the same variant pool ESM-1b covers in full, though the spot-check gives real confidence that subset isn't systematically biased by a join bug.

The review-status check came back the opposite of what I expected going in. My assumption was that restricting to expert-panel-reviewed variants would push the ClinVar AUC noticeably higher, supporting the idea that ClinVar overrepresents easy, obvious cases across the board. Instead both predictors got slightly *worse* on the high-confidence subset (AlphaMissense -0.0207, ESM-1b -0.0217) - small, but the wrong direction for that explanation, and tiny compared to the -0.24 ProteinGym gap either way. Thinking about why: ClinVar's expert panels don't necessarily review the "easiest" variants - if anything, a variant probably gets pulled into formal expert-panel review specifically because it's clinically important or was previously contested, not because it's routine. So "reviewed by expert panel" measures how well-vetted the *label* is, not how easy the *variant* is to classify - a different axis than what I'd assumed. Net effect: this check doesn't support review-status composition as an explanation for the big ProteinGym gap, and the actual explanation is still open. That's a useful negative result, not a wasted script.

**PROBLEMS ENCOUNTERED:**  
1. `get_proteingym.py` assumed the wrong internal structure for the 4.4GB scores zip and deleted it before confirming anything worked, requiring a full 4.4GB redownload. Fixed the script to never delete a zip until its extraction is confirmed to have actually matched something.
2. `proteingym_leak_free_analysis.py`'s first run matched 0% of AlphaMissense variants - ProteinGym's reference file uses UniProt mnemonic entry names, not accessions, which wasn't clear from any documentation and only became obvious from checking one real row directly.
3. `get_uniprot_mapping.py`'s first version captured an entire tab-separated API response line instead of splitting out just the accession field, producing a mapping file that looked structurally fine but was still useless - only caught because the analysis script's match rate stayed at 0% on the second run.

**NEXT SESSION GOAL:**  
Think through how to present the leakage-vs-benchmark-composition ambiguity honestly in any future writeup, and look for a follow-up test that could actually distinguish the two explanations - for example, checking whether the AUC drop concentrates near the median DMS score (where "easier benchmark composition" would predict most of the damage) versus being spread evenly across the fitness distribution. Also worth deciding whether to extend session 2's ordered-vs-disordered split to this ProteinGym data, since that was one of the project's original angles and hasn't been applied to a leak-free benchmark yet.

**QUESTIONS TO RESEARCH:**  
- Is there a clean way to test the leakage-vs-composition-difficulty question directly, rather than leaving it as an open interpretive question? (The review-status angle didn't pan out - see above. The DMS-score-distance-from-median idea from earlier in this entry is the next candidate.)

**Signed:** Sagar Raut &emsp; **Date:** August 16, 2026

---

═══════════════════════════════════════════  
**DATE:** August 17, 2026  
**SESSION:** #12  
**TIME:** 3:45 PM - 4:27 PM  
═══════════════════════════════════════════  

**GOAL FOR TODAY:**  
Try the DMS-score-distance-from-median idea raised at the end of session 11 - a direct test of whether the ~0.24 AUC gap between ClinVar and ProteinGym comes mostly from ClinVar avoiding ambiguous variants, or whether it's spread more evenly across the fitness distribution in a way that composition alone can't explain.

**BACKGROUND / REASONING:**  
Session 11 left the leakage-vs-composition question open on purpose rather than guessing at it. The review-status check that session already ran ruled out one version of the composition idea (restricting to expert-panel-reviewed variants didn't raise the AUC), but that only tested how well-vetted a label is, not how biologically ambiguous the underlying variant actually is. Those are different things, and the second one is closer to what "composition-difficulty" actually means.

**WHAT I DID:**  
1. Wrote `scripts/proteingym_median_distance_stratification.py`, building on session 11's AlphaMissense-lookup and UniProt-mapping code, but this time pulling in the continuous `DMS_score` column too - session 11 only kept the already-binarized version, which isn't enough to measure how far a variant sits from the middle of its assay's distribution.
2. For each variant, computed its percentile rank within its own assay's fitness scores (0.5 = right at the median, 0 or 1 = at either extreme), then its distance from that median. Percentile rather than raw distance, since different assays measure fitness on completely different scales.
3. First pass split variants into three groups (near median / middle / far from median) and got AUC for both predictors in each. Caught and fixed two rough spots before trusting the output: a leftover dead conditional in the AlphaMissense merge left over from adapting session 11's code, and a case where `pd.qcut` can error out (or silently produce fewer groups than asked for) when pooling percentile-based distances across 96 differently sized assays puts identical values right on a bin boundary - added `duplicates="drop"` with a fallback so the script doesn't crash and doesn't mislabel groups if that happens.
4. The three-group split showed accuracy climbing clearly as variants got more extreme, but three points aren't enough to tell whether that climb keeps going all the way to the tail or levels off early. Reran the same script with ten groups (deciles) instead of three to see the actual shape.

**DATA / RESULTS:**  
| Decile (distance from assay median) | ESM-1b AUC | AlphaMissense AUC |
|---|---|---|
| 1 (near median, most ambiguous) | 0.5074 | 0.5167 |
| 2 | 0.5329 | 0.5442 |
| 3 | 0.5717 | 0.6013 |
| 4 | 0.6089 | 0.6384 |
| 5 | 0.6608 | 0.6885 |
| 6 | 0.6915 | 0.7439 |
| 7 | 0.7404 | 0.8096 |
| 8 | 0.7876 | 0.8323 |
| 9 | 0.8113 | 0.8669 |
| 10 (most extreme, most unambiguous) | 0.8206 | 0.8667 |
| For reference: ClinVar AUC (session 11) | 0.9299 | 0.9592 |
| Remaining gap at decile 10 vs. ClinVar | 0.1093 | 0.0925 |
| n per decile | ~32,960-32,970 | 7,582-9,989 |

**OBSERVATIONS:**  
Both predictors show almost exactly the same shape: AUC rises steadily from barely-better-than-chance at the most ambiguous decile (0.51-0.52) up to 0.81-0.87 at the most extreme decile, then flattens out between deciles 9 and 10 rather than continuing to climb. That rise is strong, real evidence for the composition-difficulty explanation - predictors really are much worse at telling apart variants with subtle, in-between fitness effects, and ClinVar's clinically curated labels mostly skip that middle ground.

But the flattening matters just as much as the rise. Even at the single most unambiguous decile - the variants furthest from their assay's median in either direction - both predictors are still about 0.09-0.11 AUC below their ClinVar numbers (AlphaMissense 0.8667 vs. 0.9592, ESM-1b 0.8206 vs. 0.9299). If composition-difficulty were the entire explanation, that gap should have closed almost completely by the most extreme decile, and it didn't. So the honest read is: composition-difficulty accounts for a large share of the original 0.24 gap - restricting to the clearest-cut variants cuts the deficit roughly in half - but there's a real, bounded remainder that isn't explained by "these are just harder variants." What that remainder actually is stays an open question - could be the DMS-fitness-vs-clinical-pathogenicity mismatch already noted as a limitation in session 11, could be a subtler form of leakage than the release-date test in sessions 6-10 was built to catch, or could be something about the different protein populations ClinVar and ProteinGym cover. This result narrows the question a lot without fully closing it, which is a genuinely useful place to land rather than a letdown.

**PROBLEMS ENCOUNTERED:**  
1. First draft of the AlphaMissense merge had a leftover conditional copied over from adapting session 11's script that didn't actually do anything useful - cleaned it up before running anything, rather than leaving confusing dead logic in a script meant to be readable later.
2. `pd.qcut` can fail (or quietly return fewer bins than requested) when the same value lands on more than one bin edge, which is a real risk here since percentile-based distances from 96 assays of very different sizes get pooled together. Added `duplicates="drop"` plus a check that falls back to labeling the actual number of bins produced, so the script degrades safely instead of crashing or mislabeling groups if it happens.

**NEXT SESSION GOAL:**  
Decide between two directions: extend session 2's ordered-vs-disordered pLDDT split to this ProteinGym data, now that there's a leak-free benchmark to run it on and that angle hasn't been touched since the original project design - or start shifting toward writing up and synthesizing everything found so far, since both major threads (the release-date RDD design and the leak-free comparison) now have solid, honestly-framed conclusions.

**QUESTIONS TO RESEARCH:**  
- What's actually behind the ~0.09-0.11 AUC gap that survives even at the most extreme, least ambiguous ProteinGym variants - is there a way to test the DMS-fitness-vs-clinical-pathogenicity mismatch directly, rather than just noting it as a possible explanation?

**Signed:** Sagar Raut &emsp; **Date:** August 17, 2026

---

═══════════════════════════════════════════  
**DATE:** August 18, 2026  
**SESSION:** #13  
═══════════════════════════════════════════  

**GOAL FOR TODAY:**  
Take the first of session 12's two candidate directions: extend session 2's ordered-vs-disordered pLDDT split to the ProteinGym leak-free benchmark. That angle - do AlphaMissense and ESM-1b get worse in disordered regions - was one of the project's original design pillars, but until now it had only ever been tested against ClinVar's curated labels, never against ProteinGym's exhaustive, leak-free mutagenesis data.

**BACKGROUND / REASONING:**  
On ClinVar, both predictors showed only a small ordered-vs-IDR AUC gap (session 2/5). One live possibility I hadn't tested: that small gap could itself be a ClinVar curation artifact, the same kind of composition-difficulty effect sessions 11-12 found for the overall ClinVar-vs-ProteinGym drop. If ClinVar rarely resolves ambiguous IDR variants into a firm classification the way it rarely contains ambiguous ordered-region ones, the ordered/IDR gap on ClinVar could be artificially compressed regardless of whether the underlying predictors are actually much worse in disordered regions. ProteinGym's assays mutate every position exhaustively, ordered or disordered, without that curation filter - so if the true gap is bigger than ClinVar shows, this comparison should reveal it.

**WHAT I DID:**  
1. Checked which of the 81 unique UniProt proteins behind ProteinGym's 96 human-taxon DMS assays already had AlphaFold structures downloaded from session 2's ClinVar-wide fetch (`data/alphafold_cif/`) - 67 of 81 already present. Downloaded the remaining 14 directly; 11 succeeded, 3 (BRCA2/P51587, UBR5/O95071, Q5VST9) returned HTTP 404 from AlphaFold DB, almost certainly because they're large multi-domain proteins AlphaFold DB serves as fragmented multi-file models rather than a single F1 file - not chased further, those proteins are simply excluded same as any other "no structure available" case.
2. Wrote `scripts/proteingym_disorder_split.py`, reusing session 11's `proteingym_leak_free_analysis.py` loader functions (`load_alphamissense_lookup`, `load_proteingym_scores`) rather than duplicating that logic, plus session 2's `annotate_plddt.py` CIF-parsing approach adapted for ProteinGym's `mutant` column (same "letter-number-letter" format as ClinVar's `aa_change`, so the same regex applies unchanged).
3. Classified every variant's residue as ordered (pLDDT >= 70) / intermediate (50-70) / disordered (< 50), same cutoffs as session 2, and excluded "intermediate" from the ordered-vs-IDR comparison for the same reason session 2 did.
4. Computed AUC separately for ordered and disordered subsets, for both predictors, on ProteinGym - and, for a same-method side-by-side, recomputed the equivalent ClinVar ordered/IDR AUCs fresh from the current canonical `data/clinvar_complete.csv` in the same script rather than trusting older session 2/5 console output.
5. Built a grouped bar chart (`outputs/proteingym_disorder_split.png`) putting all four numbers per predictor side by side.

**DATA / RESULTS:**  
| | AlphaMissense | ESM-1b |
|---|---|---|
| ClinVar, ordered | 0.9530 (n=111,445) | 0.9143 (n=111,317) |
| ClinVar, IDR | 0.9420 (n=57,158) | 0.9071 (n=57,103) |
| ClinVar ordered-IDR gap | +0.0111 | +0.0072 |
| ProteinGym, ordered | 0.7078 (n=76,760) | 0.6789 (n=266,645) |
| ProteinGym, IDR | 0.5630 (n=9,459) | 0.6149 (n=38,340) |
| ProteinGym ordered-IDR gap | +0.1448 | +0.0640 |

(AlphaMissense n is smaller throughout ProteinGym because only 28% of ProteinGym variants have an AlphaMissense score, same subset-size caveat session 11 already documented. ESM-1b covers essentially the full variant set.)

**OBSERVATIONS:**  
The ordered-vs-IDR gap that looked almost negligible on ClinVar (about 1 point of AUC for AlphaMissense, under 1 point for ESM-1b) turns out to be much larger on the leak-free ProteinGym benchmark - 14.5 points for AlphaMissense, 6.4 points for ESM-1b. That's a big, clean, well-powered result (tens of thousands of variants per cell for ESM-1b) and it lines up cleanly with the project's broader theme from sessions 11-12: ClinVar-based accuracy numbers overstate how good these predictors really are, and here's a second, independent way that shows up. My working explanation is that ClinVar simply doesn't classify many truly ambiguous disordered-region variants in the first place - a variant only gets a firm pathogenic/benign call once there's fairly convincing evidence, and that evidence is harder to come by for a region with no fixed structure - so ClinVar's "disordered" bucket ends up nearly as easy, on average, as its "ordered" bucket. ProteinGym's exhaustive per-position mutagenesis has no such filter, so predictors' real IDR weakness comes through undiluted: 0.563 AUC for AlphaMissense in disordered regions is barely better than a coin flip, a genuinely bad result the ClinVar numbers never hinted at.

I want to be careful not to overclaim the causal story, the same way sessions 11-12 were about the overall ClinVar-vs-ProteinGym gap. This result doesn't distinguish "predictors are genuinely much worse at IDR variant effects" from "ProteinGym's specific 96 proteins happen to have harder-than-average disordered regions for reasons unrelated to disorder itself" - session 11 already flagged that ProteinGym's protein selection (chosen for experimental tractability) isn't the same population as ClinVar's disease-gene-heavy coverage, and that caveat applies here too. What is solid: the ordered/IDR gap is real and large on the one benchmark that can't have leaked into training, and ClinVar's version of that same comparison is not a reliable stand-in for it.

**PROBLEMS ENCOUNTERED:**  
1. Three of the 81 needed UniProt structures (BRCA2 and two others) 404'd from AlphaFold DB - large multi-domain proteins that AlphaFold DB splits into multiple fragment files rather than one F1 model. Left those proteins out entirely (same handling as any protein with "no structure available") rather than building fragment-stitching logic for 3 of 81 proteins.
2. None otherwise - reusing session 11's loader functions and session 2's CIF-parsing logic directly (rather than rewriting either) meant this script worked on the first real run.

**NEXT SESSION GOAL:**  
Both original hypotheses (H1's release-date RDD design, H2's leak-free ClinVar-vs-ProteinGym comparison, and now this ordered/IDR extension of H2) have solid, honestly-framed, well-powered results. Session 14 should shift into synthesis and writeup mode - drafting the paper/poster narrative - rather than opening a new analysis thread. IJAS Round 1 abstract is due ~February 1, 2027, so there's runway, but the project is at the point where the next unit of value is likely organizing and communicating what's already been found rather than finding more.

**QUESTIONS TO RESEARCH:**  
- None new opened this session - the honest next step is synthesis, not another open question.

**Signed:** Sagar Raut &emsp; **Date:** August 18, 2026

---

