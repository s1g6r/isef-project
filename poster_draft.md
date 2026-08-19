# Poster Content Draft — Session 16

*This is a content/layout plan for the physical poster, not the poster itself — meant to be transferred into PowerPoint/Canva/Illustrator. Text is written short and punchy for poster reading distance (a judge reads this standing several feet away, not sitting with a paper). Full explanations and hedging live in `paper_draft.md`; this pulls the condensed, poster-appropriate version of the same verified numbers.*

**Suggested layout:** standard ISEF/IJAS tri-fold or single-board layout, left-to-right, top-to-bottom reading order. Six sections below map to six poster panels.

---

## Panel 1 — Title / Hook

**Title (large, top of poster):**
Does ClinVar Circularly Validate the AI Tools Used to Build It?

**Subtitle (smaller, under title):**
Testing four missense variant effect predictors for training-data leakage and benchmark inflation

**Author line:** Sagar Raut · Computational Biology and Bioinformatics

**Hook box (optional, large pull-quote style, top corner):**
> "AI predicts which genetic mutations are dangerous. Doctors use that prediction as evidence when writing the official answer key. Then the AI gets graded against that same answer key. I tested whether that's a problem — and found something more specific than a yes-or-no."

---

## Panel 2 — The Question / Background

**Headline:** AI helped write the test it's being graded on

**Body (keep to ~80-100 words):**
Doctors classify genetic mutations as dangerous ("pathogenic") or harmless ("benign") and record that call in ClinVar, a shared public database. AI tools like AlphaMissense now predict pathogenicity instantly, and clinicians increasingly cite those predictions as evidence when writing a ClinVar classification. The field then measures how good these AI tools are by checking them against ClinVar. If the tool's own prediction shaped the label it's later graded on, that's circular — and the reported accuracy might not reflect real-world performance.

**Callout stat box:**
**~30%** of the human proteome has no fixed 3D shape (intrinsically disordered regions) — exactly where doctors have the least other evidence and may lean hardest on AI predictions.

**Small diagram suggestion:** simple circular arrow diagram — "AI predicts → clinician cites prediction as evidence → ClinVar label written → AI graded against that label" — looping back to the start.

---

## Panel 3 — Approach / Methods

**Headline:** Two independent tests, not one

**Test 1 box — "Does accuracy jump at release?"**
Each of 4 AI tools (AlphaMissense, ESM-1b, REVEL, PolyPhen-2) has a known release date. A tool can't have influenced a ClinVar label written before it existed. If circularity is real, accuracy should jump right at that date — tested with a regression discontinuity design (zoom in tight around the release date, check for a sudden break, not just a slow drift).

**Test 2 box — "Does accuracy hold up on leak-free data?"**
Re-score the same tools against ProteinGym — lab experiments (deep mutational scanning) that measure a mutation's real effect directly, with no clinician and no ClinVar involved at any point. This data could not have leaked into training or influenced any label.

**Pipeline strip (small, bottom of panel):** ClinVar (212,903 variants) + AlphaMissense + ESM-1b + REVEL + PolyPhen-2 + AlphaFold structure data (14,930 proteins) → one joined dataset → both tests run on top of it.

---

## Panel 4 — Result 1: No Sudden Jump at Release

**Headline:** Five independent checks, same answer: no leak

**Figure:** `outputs/rdd_plots.png` (all 4 predictors, 2-year window, dashed line = release date)

**Caption / key numbers (short):**
- AlphaMissense: no jump at any window tested (largest estimate +0.005, not significant)
- ESM-1b: an apparent jump (p=0.002) disappeared under a permutation test (p=0.32) — traced to an unrelated ClinVar submission surge the same month, not the tool's release
- REVEL & PolyPhen-2: not stable enough / not enough pre-release ClinVar data to test conclusively
- Checked 5 ways: raw regression, permutation test, confounder-adjusted regression, confounder-adjusted permutation test, within-protein paired matching — **all agree: no discontinuity**

**Callout stat box:**
**0 of 4** predictors show a real accuracy jump at their own release date.

---

## Panel 5 — Result 2: ClinVar Accuracy Doesn't Transfer

**Headline:** But ClinVar accuracy overstates real performance by a lot

**Figure:** `outputs/proteingym_leak_free_comparison.png`

**Caption:**
Both AlphaMissense and ESM-1b drop by almost exactly the same amount — about 0.24 AUC — moving from ClinVar to lab-measured, leak-free data. Two very differently-built AI models landing on nearly the same drop is a strong, consistent signal.

**Callout stat box (large, high-impact):**
**−0.24 AUC** — the accuracy drop for both predictors on data that could not have leaked into training

**Follow-up figure:** `outputs/proteingym_median_distance_deciles.png`

**Follow-up caption:**
Restricting to only the most clear-cut lab-measured variants closes about half the gap — ClinVar really does contain easier cases. But even the clearest-cut cases still trail ClinVar's number by 0.09–0.11 AUC. **Composition explains about half. Something else explains the rest.**

---

## Panel 6 — Result 3, Discussion, and What's Next

**Headline:** ClinVar hides a real weakness in disordered regions

**Figure:** `outputs/proteingym_disorder_split.png`

**Caption:**
The gap between how well these tools do in structured vs. disordered protein regions looks tiny on ClinVar (about 1 AUC point) — but 6 to 14 times bigger on leak-free data. AlphaMissense's accuracy on disordered-region variants (0.56 AUC) is barely better than a coin flip, something ClinVar's numbers never showed.

**Bottom-line box (largest text on this panel):**
No evidence of a sudden training-data leak at release. But strong, repeated, independent evidence that ClinVar-based accuracy does not transfer to real, lab-measured ground truth — especially in disordered regions, exactly where clinicians have the least other evidence to fall back on.

**Limitations (small print, bottom corner):**
REVEL/PolyPhen-2 not conclusively tested · lab fitness ≠ clinical pathogenicity · ProteinGym's proteins aren't the same population as ClinVar's · ~0.09–0.11 AUC gap still unexplained

**Future work (small print):**
Parsing ClinVar's full evidence-code data to test whether the gap disappears for variants with documented non-computational evidence · extending the leak-free comparison to REVEL and PolyPhen-2

**QR code placeholder:** link to `github.com/s1g6r/isef-project`

---

## Notes for building the physical poster

- All four figures above are saved in `outputs/` at 300 DPI, sized for a printed poster board rather than a screen.
- Word counts above are deliberately tight for poster reading distance — resist the urge to paste in paragraphs from `paper_draft.md` directly, they're too dense for a poster.
- Color scheme used throughout the existing figures is steelblue / coral — worth carrying that same two-color scheme into any new poster graphics (the intro diagram, panel backgrounds) for visual consistency.
