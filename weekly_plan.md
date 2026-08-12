# ISEF Weekly Plan
**Project:** Circular Evidence in ClinVar Variant Classification  
**Start:** August 11, 2026  
**IJAS Round 1 Abstract Due:** ~February 1, 2027  
**Regional Fair:** ~Mid-February 2027  

---

## Phase 1 — Pipeline and Baseline (Weeks 1-4)
*Goal: Get the core data working and reproduce published numbers. No novel claims yet.*

| Week | Dates | Task | Done? |
|---|---|---|---|
| 1 | Aug 11-17 | Download ClinVar + AlphaMissense, build join pipeline, reproduce AUC gap | YES - AUC=0.9592 |
| 2 | Aug 18-24 | Download AlphaFold pLDDT, annotate variants as ordered/disordered, rerun disorder split with real annotations | |
| 3 | Aug 25-31 | Download archival ClinVar releases (2018, 2019, 2020, 2021, 2022, 2023, 2024 at 6-month intervals), build first-classification-date table | |
| 4 | Sep 1-7 | Add ESM-1b scores to joined dataset, confirm ESM-1b join works the same way as AlphaMissense | |

**Phase 1 exit check:** You should have one clean dataset with ClinVar labels, AlphaMissense scores, ESM-1b scores, pLDDT disorder annotations, and first-classification dates. If you don't have all five, do not move to Phase 2.

---

## Phase 2 — Core RDD Analysis (Weeks 5-9)
*Goal: Run the actual hypothesis test. This is the heart of the project.*

| Week | Dates | Task | Done? |
|---|---|---|---|
| 5 | Sep 8-14 | Run naive pre/post split for AlphaMissense (Sept 2023) and ESM-1b (Aug 2021) — get first look at temporal signal | |
| 6 | Sep 15-21 | Implement regression discontinuity design (RDD) — local linear regression on either side of cutoff, test for jump | |
| 7 | Sep 22-28 | Run RDD for all 4 predictors: AlphaMissense (2023), ESM-1b (2021), REVEL (2016), PolyPhen-2 (2010) | |
| 8 | Sep 29 - Oct 5 | Confounder analysis — match on allele frequency, gene, review status. Within-protein paired analysis | |
| 9 | Oct 6-12 | Permutation test on the discontinuity — confirm the jumps are not due to chance | |

**Phase 2 exit check:** You have an RDD result for 4 predictors, a confounder-controlled estimate, and a p-value. You know whether H1 is supported.

---

## Phase 3 — ProteinGym Comparison and Tool Build (Weeks 10-14)
*Goal: H2 analysis on leak-free data, plus the reusable artifact that separates this from a reanalysis.*

| Week | Dates | Task | Done? |
|---|---|---|---|
| 10 | Oct 13-19 | Download ProteinGym DMS benchmark, measure ordered vs IDR AUC on leak-free data, compare to ClinVar numbers | |
| 11 | Oct 20-26 | Build circularity risk score: per-variant score combining classification date, review status, disorder status | |
| 12 | Oct 27 - Nov 2 | Run circularity risk score over all of ClinVar, produce decontaminated benchmark subset | |
| 13 | Nov 3-9 | Upload decontaminated benchmark to Zenodo, get DOI — this is your reusable artifact | |
| 14 | Nov 10-16 | H3 if time: parse ClinVar XML for evidence codes, test whether temporal jump disappears in non-computational evidence subset | |

**Phase 3 exit check:** You have a Zenodo dataset with a DOI. The project is now a tool, not just a finding.

---

## Phase 4 — Writing and Preparation (Weeks 15-24)
*Goal: Convert results into a submittable, defensible ISEF entry.*

| Week | Dates | Task | Done? |
|---|---|---|---|
| 15 | Nov 17-23 | Write abstract (250 words max for IJAS) — get feedback from teacher sponsor | |
| 16 | Nov 24-30 | Write full IJAS research paper — introduction and methods sections | |
| 17 | Dec 1-7 | Write results and discussion sections | |
| 18 | Dec 8-14 | Make all final figures (publication quality — not matplotlib defaults) | |
| 19 | Dec 15-21 | Full paper draft done, send to teacher sponsor and professor (if you have one) for critique | |
| 20 | Dec 22-28 | Holiday buffer — do not plan real work here, just revisions if feedback comes in | |
| 21 | Dec 29 - Jan 4 | Revise paper based on feedback | |
| 22 | Jan 5-11 | Build poster and quad chart | |
| 23 | Jan 12-18 | Record 2-minute project video | |
| 24 | Jan 19-25 | Submit IJAS Round 1 package — paper, poster, quad chart, video | |
| Buffer | Jan 26 - Feb 1 | Final fixes, confirm submission received, prepare for regional poster session | |

---

## Cold Emails — When to Send

**Start: End of Week 2 (around August 24)**

Do not email before you have a result to attach. The result you have now (AUC=0.9592 vs published 0.90-0.93) is not quite enough on its own. Wait until you have the real pLDDT disorder split done — that gives you two numbers to put in the email body and a figure to attach, which is what gets replies.

**What triggers sending:** when you finish Week 2 and have a real ordered vs disordered AUC split using pLDDT (not the ambiguous proxy), you have enough to write a compelling cold email. That is approximately August 24.

**Who to email:** faculty at Northwestern (Feinberg genetics, comp bio), UChicago (Human Genetics dept), UIC bioinformatics, Argonne National Lab, and the corresponding authors of the papers in your lit review — especially the BMC Genomics 2025 paper (Lin et al.) since your project directly extends their finding.

**How many:** send 25-30 in one batch. Do not trickle them out one at a time.

**What to attach:** one figure (the disorder split plot), one page max describing the project and the specific question you want their input on.

**Timeline after sending:** expect replies within 2 weeks. If nobody replies by September 7, send a second batch to different targets.

**Realistic outcome:** 3-5 replies, 1-2 real conversations, ideally one person willing to meet monthly and critique your analysis. Even just getting one email critique of your RDD design from a real researcher is worth putting in your paper's acknowledgments and mentioning in the judge interview.

---

## Key Dates Summary

| Date | Milestone |
|---|---|
| Aug 24 | Cold emails go out (after real pLDDT results) |
| Sep 7 | Phase 1 complete — full dataset ready |
| Oct 12 | Phase 2 complete — RDD results in hand |
| Nov 16 | Phase 3 complete — Zenodo dataset live |
| Dec 19 | Full paper draft done |
| Jan 25 | IJAS submission |
| ~Feb 1 | IJAS Round 1 abstract deadline |
| ~Feb 15 | Regional fair poster session |

---

## Warning Signs You Are Behind

- If you reach October without a working RDD result, cut H3 entirely and focus on H1+H2
- If you reach November without a professor contact, that's fine — the project does not require one
- If you reach December without a Zenodo dataset, that is a problem — prioritize it over H3
- If you reach January with no paper draft, skip H2 and write up H1 only — a clean single finding is better than a rushed multi-finding paper
