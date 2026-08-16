import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from rdd_analysis import IN_PATH, BANDWIDTHS_DAYS, build_outcome, run_rdd

# Session 9: the permutation test session 6-8 kept deferring. The OLS p-values
# in rdd_analysis.py assume the before/after linear model is correctly
# specified right around the cutoff -- if the real relationship is slightly
# curved, a "significant" jump can show up even with no true discontinuity,
# just because a straight line is a bad fit to a curve. A permutation test
# sidesteps that assumption entirely: rerun the exact same RDD procedure at
# many fake ("placebo") cutoff dates drawn from the predictor's own data, and
# see how often a jump that big happens by chance alone. If the real jump
# isn't unusual compared to the placebo jumps, it's not real -- regardless of
# what the OLS p-value said.
#
# REVEL is included too -- its RDD result flips sign across bandwidths in
# session 8, which is the same "looks significant, might just be a bad linear
# fit" situation ESM-1b turned out to be. Rather than write REVEL off from
# bandwidth instability alone, the permutation test is the more rigorous way
# to check whether any of its three bandwidth results actually survive.
# PolyPhen-2 is left out -- it can't be tested at any bandwidth (predates
# ClinVar's 2013 launch, session 8), and the placebo-cutoff floor below would
# rule out almost all of its usable date range anyway.

N_PERMUTATIONS = 2000
EXCLUDE_DAYS = 365  # placebo cutoffs within 1yr of the real cutoff are excluded,
                     # so the placebo distribution isn't contaminated by whatever
                     # is actually happening near the true release date
EARLIEST_CUTOFF = pd.Timestamp("2013-01-01")  # ClinVar's own launch date -- before
                     # this the data is a thin, sparse tail (curation-date
                     # artifacts, not real submission activity), same reason
                     # PolyPhen-2 is untestable in rdd_analysis.py
LATEST_CUTOFF_MARGIN_DAYS = 180  # keep placebo cutoffs far enough from the most
                     # recent data that there's still a real "after" window
RNG_SEED = 42

PERMUTATION_PREDICTORS = [
    {"name": "AlphaMissense", "cutoff": pd.Timestamp("2023-09-19")},
    {"name": "ESM-1b", "cutoff": pd.Timestamp("2021-08-01")},
    {"name": "REVEL", "cutoff": pd.Timestamp("2016-10-01")},
]


def get_candidate_cutoffs(valid, true_cutoff):
    """Placebo cutoffs are drawn from the predictor's own observed LastEvaluated
    dates (not a uniform date range) so the null distribution reflects where
    data actually is -- resampling real dates instead of synthetic ones."""
    dates = valid["LastEvaluated"].dropna()
    latest = dates.max() - pd.Timedelta(days=LATEST_CUTOFF_MARGIN_DAYS)
    in_range = dates[(dates >= EARLIEST_CUTOFF) & (dates <= latest)]
    far_from_real = in_range[(in_range - true_cutoff).abs() > pd.Timedelta(days=EXCLUDE_DAYS)]
    return far_from_real.values  # numpy datetime64[ns] array, duplicates kept (real weighting)


def permutation_test(valid, true_cutoff, bandwidth_days, n_perm, rng):
    real = run_rdd(valid, true_cutoff, bandwidth_days)
    if real is None:
        return None

    candidates = get_candidate_cutoffs(valid, true_cutoff)

    placebo_jumps = []
    draws = rng.choice(candidates, size=n_perm, replace=True)
    for d in draws:
        result = run_rdd(valid, pd.Timestamp(d), bandwidth_days)
        if result is not None:
            placebo_jumps.append(result["jump"])
    placebo_jumps = np.array(placebo_jumps)

    p_perm = np.mean(np.abs(placebo_jumps) >= abs(real["jump"]))
    return {
        "real_jump": real["jump"], "real_n": real["n"], "real_p_ols": real["jump_p"],
        "placebo_jumps": placebo_jumps, "p_perm": p_perm,
    }


def main():
    print("Loading Phase 1+2 complete dataset...")
    df = pd.read_csv(IN_PATH, low_memory=False)
    df["LastEvaluated"] = pd.to_datetime(df["LastEvaluated"], errors="coerce")

    rng = np.random.default_rng(RNG_SEED)

    fig, axes = plt.subplots(len(PERMUTATION_PREDICTORS), len(BANDWIDTHS_DAYS), figsize=(15, 12))

    for row, pred in enumerate(PERMUTATION_PREDICTORS):
        name, cutoff = pred["name"], pred["cutoff"]
        print(f"\n{'=' * 60}\n{name}\n{'=' * 60}")
        valid = build_outcome(df, name)

        for col, bw in enumerate(BANDWIDTHS_DAYS):
            result = permutation_test(valid, cutoff, bw, N_PERMUTATIONS, rng)
            ax = axes[row, col]
            if result is None:
                print(f"  Bandwidth +/-{bw / 365:.0f}yr: not enough data")
                ax.axis("off")
                continue

            print(f"  Bandwidth +/-{bw / 365:.0f}yr: n={result['real_n']}, "
                  f"real jump={result['real_jump']:+.4f}, "
                  f"OLS p={result['real_p_ols']:.4g}, "
                  f"permutation p={result['p_perm']:.4g} "
                  f"({N_PERMUTATIONS} placebo cutoffs)")

            ax.hist(result["placebo_jumps"], bins=40, color="steelblue", alpha=0.7,
                    label="placebo jumps")
            ax.axvline(result["real_jump"], color="coral", linewidth=2,
                       label=f"real jump ({result['real_jump']:+.3f})")
            ax.set_title(f"{name}, +/-{bw / 365:.0f}yr\npermutation p={result['p_perm']:.3g}",
                        fontsize=10)
            ax.set_xlabel("Jump in P(correct)")
            if col == 0:
                ax.set_ylabel("Count (placebo cutoffs)")
            ax.legend(fontsize=7)

    plt.tight_layout()
    plt.savefig("outputs/permutation_test_plots.png", dpi=150)
    print("\nPlot saved to outputs/permutation_test_plots.png")


if __name__ == "__main__":
    main()
