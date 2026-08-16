import numpy as np
import pandas as pd

from rdd_analysis import BANDWIDTHS_DAYS, build_outcome, run_rdd, ols_robust

# Session 10: adds review status, submitter count, and allele frequency as
# covariates directly in the RDD regression -- the standard way to add
# controls to a regression discontinuity design, rather than literal
# matching, since these are continuous/categorical variables that are
# straightforward to include as regression terms. Gene/protein identity is
# controlled separately in within_protein_paired.py via matching instead --
# gene has thousands of categories, far too many for a clean dummy-variable
# approach here, and matching holds protein identity fixed more directly
# anyway.
#
# Logic: rerun the exact same RDD as sessions 6-9, but with these covariates
# added to the design matrix. If AlphaMissense and ESM-1b's jump estimates
# barely move once review status/submitters/allele frequency are accounted
# for, that's real evidence the null results weren't hiding a confound. If
# the jump estimates DO move a lot, that itself is worth knowing -- it would
# mean the raw RDD was missing something these sessions hadn't controlled for.

IN_PATH = "data/clinvar_complete.csv"
AF_COL = "gnomAD_exomes_AF"

CONFOUNDER_PREDICTORS = [
    {"name": "AlphaMissense", "cutoff": pd.Timestamp("2023-09-19")},
    {"name": "ESM-1b", "cutoff": pd.Timestamp("2021-08-01")},
]


def build_covariates(window):
    """Review status as dummies (whatever categories are actually present in
    this window), submitter count standardized, and log10 allele frequency
    with a missing-indicator -- gnomAD doesn't cover every variant, and being
    entirely absent from gnomAD is itself informative, not just a gap to
    drop."""
    review_dummies = pd.get_dummies(window["ReviewStatus"], prefix="review",
                                     drop_first=True, dtype=float)

    submitters = window["NumberSubmitters"].astype(float)
    std = submitters.std()
    submitters_z = (submitters - submitters.mean()) / std if std > 0 else submitters * 0

    af = window[AF_COL].astype(float)
    af_missing = af.isna().astype(float)
    log_af = np.log10(af.fillna(0) + 1e-6)

    covariates = review_dummies.copy()
    covariates["submitters_z"] = submitters_z.values
    covariates["log_af"] = log_af.values
    covariates["af_missing"] = af_missing.values
    return covariates.values.astype(float)


def run_confounder_rdd(valid, cutoff, bandwidth_days):
    valid = valid.copy()
    valid["days_from_cutoff"] = (valid["LastEvaluated"] - cutoff).dt.days
    window = valid[valid["days_from_cutoff"].abs() <= bandwidth_days].copy()
    window = window.dropna(subset=["ReviewStatus", "NumberSubmitters"])

    n = len(window)
    if n < 100:
        return None

    D = (window["days_from_cutoff"] >= 0).astype(float).values
    years = window["days_from_cutoff"].values / 365.0
    y = window["correct"].values.astype(float)
    covariates = build_covariates(window)

    Xmat = np.column_stack([np.ones(n), D, years, D * years, covariates])

    try:
        beta, se, p = ols_robust(Xmat, y)
    except np.linalg.LinAlgError:
        return None

    return {"n": n, "jump": beta[1], "jump_se": se[1], "jump_p": p[1],
            "n_covariates": covariates.shape[1]}


def main():
    print("Loading dataset with gnomAD AF joined (run annotate_gnomad_af.py first if this fails)...")
    df = pd.read_csv(IN_PATH, low_memory=False)
    df["LastEvaluated"] = pd.to_datetime(df["LastEvaluated"], errors="coerce")

    for pred in CONFOUNDER_PREDICTORS:
        name, cutoff = pred["name"], pred["cutoff"]
        print(f"\n{'=' * 60}\n{name}\n{'=' * 60}")
        valid = build_outcome(df, name)

        for bw in BANDWIDTHS_DAYS:
            raw = run_rdd(valid, cutoff, bw)
            adjusted = run_confounder_rdd(valid, cutoff, bw)
            if raw is None or adjusted is None:
                print(f"  Bandwidth +/-{bw / 365:.0f}yr: not enough data")
                continue
            print(f"  Bandwidth +/-{bw / 365:.0f}yr: "
                  f"raw jump={raw['jump']:+.4f} (p={raw['jump_p']:.4g})  ->  "
                  f"confounder-adjusted jump={adjusted['jump']:+.4f} (p={adjusted['jump_p']:.4g}), "
                  f"n={adjusted['n']}, {adjusted['n_covariates']} covariate columns")


if __name__ == "__main__":
    main()
