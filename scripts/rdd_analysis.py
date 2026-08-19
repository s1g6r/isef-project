import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

# Session 6: the actual regression discontinuity design (RDD) the naive pre/post
# splits in sessions 2 and 5 were only ever a warm-up for. Naive splits compare
# two big averages and can't tell a real discontinuity apart from a slow trend
# (sessions 2 and 5 both ran into exactly that). RDD instead only looks at data
# close to the cutoff, fits a separate local linear trend on each side, and
# tests whether the trend jumps right at the cutoff.
#
# Per-observation outcome: AUC itself isn't something you can regress a single
# variant on, so the outcome here is "did the predictor's own published call
# match the ClinVar label" (1/0) for that variant. Regressing this 0/1 outcome
# on time, separately on each side of the release date, and testing for a jump
# at the boundary is the standard way to run an RDD when the underlying metric
# (AUC) is itself an aggregate over many observations.
#
# Implemented as a single OLS regression on the windowed data:
#   correct = b0 + b1*D + b2*years_from_cutoff + b3*(D*years_from_cutoff)
# where D = 1 if after the release date. b1 IS the RDD estimate -- the jump in
# P(correct) right at the cutoff. Standard errors are heteroskedasticity-robust
# (HC1), computed by hand with numpy so this doesn't need statsmodels installed.
#
# This gives a point estimate and a rough significance check. It is NOT the
# permutation test that will confirm whether the jump could be due to chance --
# that's a separate, more careful step for later, once all four predictors are
# in (the OLS p-value here assumes the linear-near-cutoff model is correctly
# specified, which a permutation test doesn't have to assume).

IN_PATH = "data/clinvar_phase1_and_2_complete.csv"
OUT_PLOT = "outputs/rdd_plots.png"
BANDWIDTHS_DAYS = [365, 730, 1095]  # 1yr, 2yr, 3yr -- robustness check, not cherry-picked
PLOT_BANDWIDTH_DAYS = 730

PREDICTORS = [
    {"name": "AlphaMissense", "cutoff": pd.Timestamp("2023-09-19")},
    {"name": "ESM-1b", "cutoff": pd.Timestamp("2021-08-01")},
    {"name": "REVEL", "cutoff": pd.Timestamp("2016-10-01")},
    {"name": "PolyPhen-2", "cutoff": pd.Timestamp("2010-02-01")},
]


def ols_robust(X, y):
    """OLS with HC1 heteroskedasticity-robust standard errors, no statsmodels needed."""
    n, k = X.shape
    # np.errstate here suppresses a known false-positive RuntimeWarning from
    # numpy + Apple's Accelerate BLAS backend on Apple Silicon (arm64): certain
    # matrix-vector shapes trigger a spurious "divide by zero"/"overflow"
    # warning even though the result is correct. Verified by hand against a
    # manual (non-BLAS) computation on this exact data -- values matched
    # exactly, no NaN/Inf. Not something to silently trust elsewhere without
    # checking; confirmed once here specifically.
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        XtX_inv = np.linalg.inv(X.T @ X)
        beta = XtX_inv @ X.T @ y
        resid = y - X @ beta
        meat = X.T @ (X * (resid ** 2)[:, None])
        V = XtX_inv @ meat @ XtX_inv * (n / (n - k))
    se = np.sqrt(np.diag(V))
    p = 2 * (1 - norm.cdf(np.abs(beta / se)))
    return beta, se, p


def build_outcome(df, name):
    """Per-variant 1/0 correctness outcome, using each predictor's own published call.
    REVEL and Polyphen2-HVAR don't have a three-way class the way AlphaMissense
    does, so both use a 0.5 score threshold -- REVEL's score is explicitly
    designed as a pathogenicity probability (Ioannidis et al. 2016), and
    Polyphen2's score is a naive-Bayes posterior probability of "damaging"
    (Adzhubei et al. 2010), so 0.5 matches each model's own intended meaning
    rather than being an arbitrary cutoff."""
    if name == "AlphaMissense":
        valid = df.dropna(subset=["am_class", "label", "LastEvaluated"]).copy()
        valid = valid[valid["am_class"] != "ambiguous"]
        valid["predicted_pathogenic"] = (valid["am_class"] == "likely_pathogenic").astype(int)
    elif name == "ESM-1b":  # no published three-way class, use the -7.5 LLR cutoff directly
        valid = df.dropna(subset=["esm1b_llr", "label", "LastEvaluated"]).copy()
        valid["predicted_pathogenic"] = (valid["esm1b_llr"] <= -7.5).astype(int)
    elif name == "REVEL":
        valid = df.dropna(subset=["REVEL_score", "label", "LastEvaluated"]).copy()
        valid["predicted_pathogenic"] = (valid["REVEL_score"] >= 0.5).astype(int)
    else:  # PolyPhen-2 (HVAR)
        valid = df.dropna(subset=["Polyphen2_HVAR_score", "label", "LastEvaluated"]).copy()
        valid["predicted_pathogenic"] = (valid["Polyphen2_HVAR_score"] >= 0.5).astype(int)
    valid["correct"] = (valid["predicted_pathogenic"] == valid["label"]).astype(int)
    return valid


def run_rdd(valid, cutoff, bandwidth_days):
    valid = valid.copy()
    valid["days_from_cutoff"] = (valid["LastEvaluated"] - cutoff).dt.days
    window = valid[valid["days_from_cutoff"].abs() <= bandwidth_days]

    n = len(window)
    if n < 100:
        return None

    D = (window["days_from_cutoff"] >= 0).astype(float).values
    years = window["days_from_cutoff"].values / 365.0
    y = window["correct"].values.astype(float)

    Xmat = np.column_stack([np.ones(n), D, years, D * years])
    beta, se, p = ols_robust(Xmat, y)

    return {"n": n, "jump": beta[1], "jump_se": se[1], "jump_p": p[1],
            "beta": beta, "window": window}


def main():
    print("Loading Phase 1+2 complete dataset...")
    df = pd.read_csv(IN_PATH, low_memory=False)
    df["LastEvaluated"] = pd.to_datetime(df["LastEvaluated"], errors="coerce")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for ax, pred in zip(axes, PREDICTORS):
        name, cutoff = pred["name"], pred["cutoff"]
        print(f"\n{'=' * 60}\n{name}\n{'=' * 60}")

        valid = build_outcome(df, name)
        print(f"Rows with a valid predicted call and label: {len(valid)}")

        main_result = None
        for bw in BANDWIDTHS_DAYS:
            result = run_rdd(valid, cutoff, bw)
            if result is None:
                print(f"  Bandwidth +/-{bw / 365:.0f}yr: not enough data")
                continue
            print(f"  Bandwidth +/-{bw / 365:.0f}yr: n={result['n']}, "
                  f"jump in P(correct)={result['jump']:+.4f} "
                  f"(robust SE={result['jump_se']:.4f}, p={result['jump_p']:.4g})")
            if bw == PLOT_BANDWIDTH_DAYS:
                main_result = result

        window = main_result["window"]
        beta = main_result["beta"]
        bw_years = PLOT_BANDWIDTH_DAYS / 365

        window = window.copy()
        window["month_bin"] = (window["days_from_cutoff"] // 30) * 30
        binned = window.groupby("month_bin")["correct"].agg(["mean", "count"])
        binned = binned[binned["count"] >= 10]

        ax.scatter(binned.index / 365, binned["mean"],
                   s=(binned["count"] / binned["count"].max() * 150).clip(lower=10),
                   alpha=0.6, color="steelblue", label="monthly bin (size ~ n)")

        before_x = np.linspace(-bw_years, 0, 50)
        after_x = np.linspace(0, bw_years, 50)
        ax.plot(before_x, beta[0] + beta[2] * before_x, color="coral", linewidth=2, label="local linear fit")
        ax.plot(after_x, (beta[0] + beta[1]) + (beta[2] + beta[3]) * after_x, color="coral", linewidth=2)

        ax.axvline(x=0, color="black", linestyle="--", linewidth=1.5, label=f"{name} release")
        ax.set_xlabel("Years relative to release date")
        ax.set_ylabel("P(predictor's own call matches ClinVar label)")
        ax.set_title(f"{name} RDD (bandwidth = {bw_years:.0f}yr)")
        ax.set_ylim(0.5, 1.0)
        ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(OUT_PLOT, dpi=300)
    print(f"\nPlot saved to {OUT_PLOT}")

    print("\n" + "=" * 60)
    print("Reminder: the p-values above assume the local-linear model is right")
    print("near the cutoff. A permutation test (randomizing where the cutoff")
    print("falls and re-running this many times) is the more rigorous check")
    print("for whether these jumps could be due to chance -- that's the next")
    print("real step now that all four predictors are tested here.")
    print("\nAlso remember the session 7 lesson: before trusting any single")
    print("jump, check whether an unrelated predictor that didn't exist yet")
    print("shows the same pattern at that calendar date (like the shared")
    print("Aug 2021 ClinVar batch that explained ESM-1b's naive jump).")


if __name__ == "__main__":
    main()
