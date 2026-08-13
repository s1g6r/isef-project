import os
import urllib.request

# Session 4 (weekly_plan.md Week 4): add ESM-1b scores, predictor #2 of the
# 4-predictor RDD design (AlphaMissense 2023, ESM-1b 2021, REVEL 2016, PolyPhen-2 2010).
#
# ESM-1b (Brandes et al. 2023, Nature Genetics) has no single bulk genome-wide
# download like AlphaMissense's flat tsv.gz -- the official distribution is an
# interactive web portal. BUT the portal's own Hugging Face Space repo hosts the
# complete precomputed catalog as a public file, no login/registration needed:
#   ALL_hum_isoforms_ESM1b_LLR.zip (1.34 GB) -- one CSV per protein isoform,
#   named {uniprot_accession}_LLR.csv, each a position x amino-acid LLR score matrix.
# This is the actual data backing the portal (confirmed by reading the portal's
# own app.py source), not a scraped copy.

FILES = {
    "ALL_hum_isoforms_ESM1b_LLR.zip":
        "https://huggingface.co/spaces/ntranoslab/esm_variants/resolve/main/ALL_hum_isoforms_ESM1b_LLR.zip",
}
OUT_DIR = "data/esm1b"


def download(fname, url):
    dest = os.path.join(OUT_DIR, fname)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"{fname}: already downloaded ({os.path.getsize(dest) / (1024**3):.2f} GB), skipping")
        return
    print(f"{fname}: downloading from {url} ...")
    print("(1.34 GB -- this will take a while depending on connection speed)")

    def progress(block_num, block_size, total_size):
        downloaded = block_num * block_size
        if total_size > 0 and block_num % 200 == 0:
            pct = min(100, downloaded * 100 / total_size)
            print(f"  {pct:.0f}% ({downloaded / (1024**3):.2f} / {total_size / (1024**3):.2f} GB)")

    try:
        urllib.request.urlretrieve(url, dest, reporthook=progress)
        print(f"{fname}: done ({os.path.getsize(dest) / (1024**3):.2f} GB)")
    except Exception as e:
        if os.path.exists(dest):
            os.remove(dest)
        print(f"{fname}: FAILED: {e} -- safe to rerun this script to retry")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for fname, url in FILES.items():
        download(fname, url)


if __name__ == "__main__":
    main()
