import os
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

# Week 3 goal (weekly_plan.md): archival ClinVar releases at 6-month intervals,
# 2018-2026, to later build a first-classification-date table (needed for H3
# and for confounder control -- knowing when a variant was FIRST classified,
# not just when it was LAST evaluated).
#
# NCBI's archive splits by era:
#   2018-2024 -> https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/{year}/variant_summary_{year}-{month}.txt.gz
#   2025+     -> https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_{year}-{month}.txt.gz
# Already have 2019-01 from session 1 (data/variant_summary_2019-01.txt.gz) --
# that one is left where it is; everything else downloads into data/archive/.

BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive"
OUT_DIR = "data/archive"
FLAT_YEARS = {2025, 2026}  # years served directly under archive/, not archive/{year}/

# (year, month) at ~6-month spacing, 2018 through the most recent full release
# available as of this session (2026-07). The "current" file (data/variant_summary.txt.gz,
# no archive suffix) already covers ~Aug 2026 and is used as-is downstream.
RELEASES = [
    (2018, 1), (2018, 7),
    (2019, 1), (2019, 7),
    (2020, 1), (2020, 7),
    (2021, 1), (2021, 7),
    (2022, 1), (2022, 7),
    (2023, 1), (2023, 7),
    (2024, 1), (2024, 7),
    (2025, 1), (2025, 7),
    (2026, 1), (2026, 7),
]


def url_for(year, month):
    fname = f"variant_summary_{year}-{month:02d}.txt.gz"
    if year in FLAT_YEARS:
        return f"{BASE}/{fname}", fname
    return f"{BASE}/{year}/{fname}", fname


def download_one(year, month):
    url, fname = url_for(year, month)
    dest = os.path.join(OUT_DIR, fname)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return fname, "skipped (already have)"
    try:
        urllib.request.urlretrieve(url, dest)
        size_mb = os.path.getsize(dest) / (1024 * 1024)
        return fname, f"downloaded ({size_mb:.1f} MB)"
    except Exception as e:
        if os.path.exists(dest):
            os.remove(dest)  # don't leave a partial/corrupt file behind
        return fname, f"FAILED: {e}"


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"Downloading {len(RELEASES)} archival ClinVar releases into {OUT_DIR}/ ...")
    print("(This is several GB total -- older releases are small, 2025-2026 ones are ~250-450MB each.)\n")

    results = []
    with ThreadPoolExecutor(max_workers=4) as ex:
        futures = {ex.submit(download_one, y, m): (y, m) for y, m in RELEASES}
        for fut in as_completed(futures):
            fname, status = fut.result()
            print(f"  {fname}: {status}")
            results.append((fname, status))

    failed = [r for r in results if r[1].startswith("FAILED")]
    print(f"\nDone. {len(results) - len(failed)}/{len(results)} succeeded.")
    if failed:
        print("Failed (safe to just rerun this script to retry):")
        for fname, status in failed:
            print(f"  {fname}: {status}")


if __name__ == "__main__":
    main()
