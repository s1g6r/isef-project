import os
import urllib.request
import zipfile

import pandas as pd

# Session 11: downloads ProteinGym's leak-free DMS (deep mutational scanning)
# benchmark for H2 -- experimentally measured variant effects from lab
# assays, not from ClinVar or any other clinical database, so nothing here
# could have leaked into a predictor's training set the way ClinVar might
# have. Scope: AlphaMissense and ESM-1b only. AlphaMissense's own genome-wide
# release (already downloaded in session 1, data/AlphaMissense_hg38.tsv.gz)
# already covers ProteinGym's human variants directly via its uniprot_id +
# protein_variant columns, so no separate AlphaMissense download is needed
# here. REVEL and PolyPhen-2 aren't in ProteinGym's baselines and would
# require running their actual model pipelines, which is out of scope.
#
# Three downloads:
#   1. reference_files/DMS_substitutions.csv (small, straight from GitHub) --
#      per-assay metadata including UniProt_ID and taxon. Used to filter to
#      human-taxon assays only, since ClinVar and AlphaMissense are both
#      human-only, and to know which per-assay files are worth keeping.
#   2. DMS_ProteinGym_substitutions.zip (~1GB) -- the actual experimental DMS
#      scores, one CSV per assay. Only human-taxon assay files get extracted.
#   3. zero_shot_substitutions_scores.zip (~4.4GB) -- every baseline model's
#      predictions on every assay, 43 models total. Only files with "esm1b"
#      in the name get extracted; the internal folder layout isn't
#      documented anywhere reliable, so this greps the real zip contents by
#      name rather than assuming a structure. Everything else -- including
#      both big zips once extraction is done -- gets deleted immediately,
#      since keeping several GB of scores for 42 models this project will
#      never use has no reason to sit on disk.

BASE_URL = "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3"
REFERENCE_URL = "https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/reference_files/DMS_substitutions.csv"
OUT_DIR = "data/proteingym"


def download(url, dest):
    if os.path.exists(dest):
        print(f"  {dest} already exists, skipping download.")
        return
    print(f"  Downloading {url}")
    print(f"    -> {dest}")
    tmp = dest + ".part"
    urllib.request.urlretrieve(url, tmp)
    os.rename(tmp, dest)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    ref_path = os.path.join(OUT_DIR, "DMS_substitutions_reference.csv")
    download(REFERENCE_URL, ref_path)

    dms_zip = os.path.join(OUT_DIR, "DMS_ProteinGym_substitutions.zip")
    download(f"{BASE_URL}/DMS_ProteinGym_substitutions.zip", dms_zip)

    scores_zip = os.path.join(OUT_DIR, "zero_shot_substitutions_scores.zip")
    download(f"{BASE_URL}/zero_shot_substitutions_scores.zip", scores_zip)

    print("\nReading reference file to find human-taxon assays...")
    ref = pd.read_csv(ref_path)
    human_ref = ref[ref["taxon"].str.lower() == "human"]
    human_filenames = set(human_ref["DMS_filename"])
    print(f"  {len(human_filenames)} of {len(ref)} DMS assays are human-taxon -- "
          f"only these are relevant to a human-ClinVar-based comparison.")

    dms_dir = os.path.join(OUT_DIR, "DMS_assays")
    if os.path.isdir(dms_dir) and os.listdir(dms_dir):
        print(f"\n{dms_dir} already has files, skipping extraction.")
    else:
        print("\nExtracting human-taxon DMS assay files...")
        os.makedirs(dms_dir, exist_ok=True)
        with zipfile.ZipFile(dms_zip) as zf:
            names = zf.namelist()
            matched = [n for n in names if os.path.basename(n) in human_filenames]
            print(f"  Found {len(matched)} matching files in the zip (expected {len(human_filenames)}).")
            for n in matched:
                zf.extract(n, dms_dir)
        print(f"  Extracted to {dms_dir}")

    scores_dir = os.path.join(OUT_DIR, "zero_shot_scores")
    if os.path.isdir(scores_dir) and os.listdir(scores_dir):
        print(f"\n{scores_dir} already has files, skipping extraction.")
    else:
        # First attempt assumed the scores zip had one file per model (like
        # "esm1b/DMS_id.csv") and turned up zero matches -- but the zip has
        # exactly 217 entries, the same as the number of DMS assays, which
        # means it's actually one file PER ASSAY with all 43 models as
        # columns inside (wide format), not one file per model. So this
        # extracts by the same DMS_filename matching the DMS zip already
        # used correctly, and the ESM-1b column gets picked out once a file
        # is actually opened and its real columns are visible.
        print("\nExtracting human-taxon files from the scores zip (same filename matching as the DMS zip)...")
        os.makedirs(scores_dir, exist_ok=True)
        with zipfile.ZipFile(scores_zip) as zf:
            names = zf.namelist()
            matched = [n for n in names if os.path.basename(n) in human_filenames]
            print(f"  Found {len(matched)} matching files out of {len(names)} total in the zip.")
            if matched:
                for n in matched:
                    zf.extract(n, scores_dir)
                print(f"  Extracted to {scores_dir}")
            else:
                print("  0 matches -- NOT deleting the zip. Inspect its real contents before doing anything else:")
                print("  python3 -c \"import zipfile; z = zipfile.ZipFile('data/proteingym/"
                      "zero_shot_substitutions_scores.zip'); print(z.namelist()[:20])\"")

    print("\nDeleting the DMS zip now that its extraction succeeded "
          "(the scores zip is only deleted once its extraction is confirmed working)...")
    if os.path.exists(dms_zip):
        os.remove(dms_zip)
        print(f"  Removed {dms_zip}")
    if os.path.isdir(scores_dir) and os.listdir(scores_dir) and os.path.exists(scores_zip):
        os.remove(scores_zip)
        print(f"  Removed {scores_zip}")

    print("\nDone.")


if __name__ == "__main__":
    main()
