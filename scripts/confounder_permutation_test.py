import numpy as np
import pandas as pd

from rdd_analysis import BANDWIDTHS_DAYS, build_outcome
from confounder_rdd import IN_PATH, CONFOUNDER_PREDICTORS, run_confounder_rdd
from permutation_test import get_candidate_cutoffs

# Session 10 continued: confounder_rdd.py produced two OLS p-values that look
# newly significant compared to the raw RDD -- AlphaMissense at 3yr (raw
# p=0.208 -> adjusted p=0.039) and ESM-1b at 1yr (stays significant, jump
# shrinks from +0.0237 to +0.0179). Session 9 already proved an OLS p-value
# in this exact design can be wrong: ESM-1b's original 1yr raw jump looked
# significant (p=0.0024) and did not survive a permutation test (p=0.32).
# The same skepticism has to apply here too, or the confounder analysis
# would be held to a lower standard than the rest of the project.
#
# Reuses the same placebo-cutoff logic as permutation_test.py (draw fake
# cutoffs from the predictor's own real dates, 2013+ only, excluding dates
# within a year of the true cutoff), just calling the confounder-adjusted
# RDD instead of the raw one at each placebo cutoff. N_PERM is lower than
# permutation_test.py's 2000 (1000 here) since each fit is slower with the
# extra covariates -- still enough resolution for a p-value down to ~0.001.

N_PERM = 1000
RNG_SEED = 42


def permutation_test_confounder(valid, true_cutoff, bandwidth_days, n_perm, rng):
    real = run_confounder_rdd(valid, true_cutoff, bandwidth_days)
    if real is None:
        return None

    candidates = get_candidate_cutoffs(valid, true_cutoff)
    placebo_jumps = []
    draws = rng.choice(candidates, size=n_perm, replace=True)
    for d in draws:
        result = run_confounder_rdd(valid, pd.Timestamp(d), bandwidth_days)
        if result is not None:
            placebo_jumps.append(result["jump"])
    placebo_jumps = np.array(placebo_jumps)

    p_perm = np.mean(np.abs(placebo_jumps) >= abs(real["jump"]))
    return {"real_jump": real["jump"], "real_n": real["n"], "real_p_ols": real["jump_p"], "p_perm": p_perm}


def main():
    print("Loading dataset with gnomAD AF joined...")
    df = pd.read_csv(IN_PATH, low_memory=False)
    df["LastEvaluated"] = pd.to_datetime(df["LastEvaluated"], errors="coerce")

    rng = np.random.default_rng(RNG_SEED)

    for pred in CONFOUNDER_PREDICTORS:
        name, cutoff = pred["name"], pred["cutoff"]
        print(f"\n{'=' * 60}\n{name}\n{'=' * 60}")
        valid = build_outcome(df, name)

        for bw in BANDWIDTHS_DAYS:
            result = permutation_test_confounder(valid, cutoff, bw, N_PERM, rng)
            if result is None:
                print(f"  Bandwidth +/-{bw / 365:.0f}yr: not enough data")
                continue
            print(f"  Bandwidth +/-{bw / 365:.0f}yr: n={result['real_n']}, "
                  f"confounder-adjusted jump={result['real_jump']:+.4f}, "
                  f"OLS p={result['real_p_ols']:.4g}, permutation p={result['p_perm']:.4g} "
                  f"({N_PERM} placebo cutoffs)")


if __name__ == "__main__":
    main()
