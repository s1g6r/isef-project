import numpy as np
import pandas as pd
from scipy.stats import binomtest

from rdd_analysis import build_outcome

# Session 10: the within-protein paired analysis from the plan. Match each
# post-cutoff variant to its nearest pre-cutoff variant in the same protein
# (by residue position, "resnum" -- already annotated via AlphaFold/UniProt
# back in session 2), so the comparison is variant-vs-variant within the same
# protein instead of pooled across every protein at once. This controls for
# gene/protein identity directly through matching -- some genes are just
# intrinsically easier to classify than others regardless of when a given
# variant was evaluated, and pooling across genes (like the plain RDD does)
# can't tell that apart from a real temporal effect if gene composition
# shifts over time near a cutoff.
#
# McNemar's test is the standard test for this kind of paired binary
# before/after design. It only looks at pairs where the before and after
# outcomes disagree, and tests whether disagreements are symmetric (equally
# likely to flip from incorrect-to-correct as the reverse) rather than
# needing any linear model to be correctly specified. It shares no machinery
# with the RDD or the permutation test, so this is a genuinely independent
# third check on whether AlphaMissense and ESM-1b's null results hold up.

IN_PATH = "data/clinvar_complete.csv"
BANDWIDTH_DAYS = 730  # same as the RDD's main plotting bandwidth (sessions 6-9)

WITHIN_PROTEIN_PREDICTORS = [
    {"name": "AlphaMissense", "cutoff": pd.Timestamp("2023-09-19")},
    {"name": "ESM-1b", "cutoff": pd.Timestamp("2021-08-01")},
]


def match_within_protein(valid, cutoff, bandwidth_days):
    valid = valid.dropna(subset=["uniprot", "resnum"]).copy()
    valid["days_from_cutoff"] = (valid["LastEvaluated"] - cutoff).dt.days
    window = valid[valid["days_from_cutoff"].abs() <= bandwidth_days]

    before = window[window["days_from_cutoff"] < 0]
    after = window[window["days_from_cutoff"] >= 0]

    # Pre-group "before" variants by protein once, instead of re-filtering
    # the whole before-set for every protein in the loop below.
    before_by_protein = {uid: (g["resnum"].values, g["correct"].values)
                          for uid, g in before.groupby("uniprot")}

    pairs = []
    proteins_matched = set()
    for uniprot_id, after_group in after.groupby("uniprot"):
        if uniprot_id not in before_by_protein:
            continue
        before_resnums, before_correct = before_by_protein[uniprot_id]
        proteins_matched.add(uniprot_id)
        # Nearest-neighbor by residue position, with replacement -- a given
        # "before" variant can be reused as the match for more than one
        # "after" variant. A cleaner without-replacement version is possible
        # but adds real complexity for a first pass at this analysis.
        for resnum, correct in zip(after_group["resnum"].values, after_group["correct"].values):
            idx = np.argmin(np.abs(before_resnums - resnum))
            pairs.append((before_correct[idx], correct))

    return pairs, len(proteins_matched)


def mcnemar(pairs):
    pairs = np.array(pairs)
    before_correct = pairs[:, 0]
    after_correct = pairs[:, 1]

    b = int(np.sum((before_correct == 1) & (after_correct == 0)))  # correct -> incorrect
    c = int(np.sum((before_correct == 0) & (after_correct == 1)))  # incorrect -> correct
    concordant = len(pairs) - b - c

    if b + c == 0:
        return {"n_pairs": len(pairs), "b": b, "c": c, "concordant": concordant,
                "paired_jump": 0.0, "p": 1.0}

    result = binomtest(c, n=b + c, p=0.5)
    paired_jump = (c - b) / len(pairs)
    return {"n_pairs": len(pairs), "b": b, "c": c, "concordant": concordant,
            "paired_jump": paired_jump, "p": result.pvalue}


def main():
    print("Loading dataset...")
    df = pd.read_csv(IN_PATH, low_memory=False)
    df["LastEvaluated"] = pd.to_datetime(df["LastEvaluated"], errors="coerce")

    for pred in WITHIN_PROTEIN_PREDICTORS:
        name, cutoff = pred["name"], pred["cutoff"]
        print(f"\n{'=' * 60}\n{name}\n{'=' * 60}")
        valid = build_outcome(df, name)
        pairs, n_proteins = match_within_protein(valid, cutoff, BANDWIDTH_DAYS)

        if len(pairs) < 20:
            print(f"  Not enough within-protein pairs to test (n={len(pairs)}).")
            continue

        result = mcnemar(pairs)
        print(f"  {n_proteins} proteins matched, {result['n_pairs']} paired variants")
        print(f"  concordant pairs (agreed)={result['concordant']}, "
              f"discordant: correct->incorrect={result['b']}, incorrect->correct={result['c']}")
        print(f"  paired jump estimate={result['paired_jump']:+.4f}, McNemar exact p={result['p']:.4g}")


if __name__ == "__main__":
    main()
