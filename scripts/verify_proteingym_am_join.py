import random

import pandas as pd

# Session 11 continued: before trusting AlphaMissense's leak-free AUC
# (built on a join that only matched 28% of ProteinGym's variants), confirm
# the matched rows are actually correct and not silently wrong. The risk:
# DMS papers sometimes number residues against a lab construct instead of
# the full UniProt canonical sequence (e.g. a truncated protein numbered
# starting at 1), which would shift every position and could, in principle,
# produce coincidental wrong matches rather than just failed ones.
#
# Check: for a random sample of matched rows, parse the "mutant" string
# (e.g. "A119D" -> ref=A, position=119, alt=D) and confirm the reference
# amino acid actually appears at that position in the assay's own
# target_seq (from the reference file) -- if the join's position numbering
# were wrong, this would fail for most or all sampled rows.

REFERENCE_PATH = "data/proteingym/DMS_substitutions_reference.csv"
JOINED_PATH = "data/proteingym_leak_free_joined.csv"
N_SAMPLE = 30
SEED = 42


def parse_mutant(mutant):
    ref_aa = mutant[0]
    alt_aa = mutant[-1]
    position = int(mutant[1:-1])
    return ref_aa, position, alt_aa


def main():
    ref = pd.read_csv(REFERENCE_PATH)
    target_seq_by_filename = dict(zip(ref["DMS_filename"], ref["target_seq"]))

    joined = pd.read_csv(JOINED_PATH)
    matched = joined.dropna(subset=["am_pathogenicity"])
    print(f"{len(matched):,} matched rows to sample from.")

    random.seed(SEED)
    sample = matched.sample(n=min(N_SAMPLE, len(matched)), random_state=SEED)

    n_ok, n_bad = 0, 0
    for _, row in sample.iterrows():
        target_seq = target_seq_by_filename.get(row["dms_filename"])
        if target_seq is None:
            print(f"  NO target_seq FOUND for {row['dms_filename']} -- skipping")
            continue

        ref_aa, position, alt_aa = parse_mutant(row["mutant"])
        if position < 1 or position > len(target_seq):
            print(f"  {row['dms_filename']} {row['mutant']}: position {position} "
                  f"out of range for a {len(target_seq)}-residue sequence -- BAD")
            n_bad += 1
            continue

        actual_aa = target_seq[position - 1]
        if actual_aa == ref_aa:
            n_ok += 1
        else:
            n_bad += 1
            print(f"  {row['dms_filename']} {row['mutant']}: expected '{ref_aa}' at position "
                  f"{position}, target_seq has '{actual_aa}' -- BAD")

    print(f"\n{n_ok}/{n_ok + n_bad} sampled matches have the correct reference amino acid "
          f"at the claimed position.")
    if n_bad == 0:
        print("Join looks correct -- the 28% match rate is real coverage, not a numbering bug.")
    else:
        print("Something is wrong with the join's position alignment -- do not trust the "
              "AlphaMissense leak-free AUC until this is resolved.")


if __name__ == "__main__":
    main()
