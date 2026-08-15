import os
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

# Downloads dbNSFP v4.1a (academic-use branch, includes REVEL and Polyphen2)
# from its Zenodo mirror -- no institutional-email registration needed, unlike
# the official dbnsfp.org site. Verified this mirror is the full academic
# version (not the commercial-safe v4.1c branch, which strips out Polyphen2,
# REVEL, VEST, ClinPred, CADD, LINSIGHT, and GenoCanyon precisely because
# those are the licensed predictors -- confirmed the record page explicitly
# says so before trusting it).
#
# One download solves both remaining predictors (#3 REVEL, #4 PolyPhen-2) at
# once, in genomic chrom/pos/ref/alt coordinates -- the same join key already
# used for every other predictor in this project, so no coordinate-mapping
# step is needed (unlike PolyPhen-2's own WHRESS download, which uses RefSeq
# protein-relative positions and would need its own translation layer).
#
# ~25.2GB total across per-chromosome files. Resumable, concurrent, stdlib
# only -- same pattern as get_archival_clinvar.py and get_esm1b_scores.py.

RECORD_ID = "4323592"
BASE_URL = f"https://zenodo.org/records/{RECORD_ID}/files"
OUT_DIR = "data/dbnsfp"
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]  # skip chrM, not relevant to nuclear missense variants
MAX_WORKERS = 4


def download_one(chrom):
    fname = f"dbNSFP4.1a_variant.chr{chrom}.gz"
    url = f"{BASE_URL}/{fname}?download=1"
    out_path = os.path.join(OUT_DIR, fname)

    if os.path.exists(out_path):
        return f"chr{chrom}: already downloaded, skipping"

    tmp_path = out_path + ".part"
    try:
        urllib.request.urlretrieve(url, tmp_path)
        os.rename(tmp_path, out_path)
        size_gb = os.path.getsize(out_path) / (1024 ** 3)
        return f"chr{chrom}: done ({size_gb:.2f} GB)"
    except Exception as e:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        return f"chr{chrom}: FAILED ({e})"


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"Downloading dbNSFP v4.1a variant files for {len(CHROMS)} chromosomes "
          f"({MAX_WORKERS} at a time -- this is ~25GB total, will take a while)...")

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(download_one, c): c for c in CHROMS}
        for future in as_completed(futures):
            print(future.result())

    print("\nDone.")


if __name__ == "__main__":
    main()
