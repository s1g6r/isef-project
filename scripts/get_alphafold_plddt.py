import os
import urllib.request
import urllib.error
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# Session 2: real disorder annotation via AlphaFold pLDDT.
# MobiDB (get_disorder.py) is dead - API returns 200 with headers but zero rows.
# AlphaFold DB per-protein CIF files store pLDDT in the B-factor column (B_iso_or_equiv),
# one value per residue. We only need the proteins that actually appear in our
# ClinVar+AlphaMissense joined dataset, not the whole human proteome.

OUT_DIR = "data/alphafold_cif"
URL_TMPL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.cif"
MAX_WORKERS = 8

os.makedirs(OUT_DIR, exist_ok=True)

print("Loading joined dataset to get unique UniProt accessions...")
df = pd.read_csv("data/clinvar_am_joined.csv", usecols=["uniprot"])
accessions = sorted(df["uniprot"].dropna().unique())
print(f"Unique accessions needed: {len(accessions)}")

# Skip anything already downloaded (makes this resumable if interrupted)
todo = [acc for acc in accessions if not os.path.exists(f"{OUT_DIR}/{acc}.cif")]
print(f"Already downloaded: {len(accessions) - len(todo)}")
print(f"Left to fetch: {len(todo)}")


def fetch(acc):
    url = URL_TMPL.format(acc=acc)
    dest = f"{OUT_DIR}/{acc}.cif"
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            data = resp.read()
        with open(dest, "wb") as f:
            f.write(data)
        return acc, "ok"
    except urllib.error.HTTPError as e:
        return acc, f"http_{e.code}"
    except Exception as e:
        return acc, f"error_{type(e).__name__}"


ok_count = 0
missing = []  # accessions with no AlphaFold model (usually HTTP 404)
errors = []   # transient failures worth retrying

with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
    futures = {pool.submit(fetch, acc): acc for acc in todo}
    done = 0
    for future in as_completed(futures):
        acc, status = future.result()
        done += 1
        if status == "ok":
            ok_count += 1
        elif status.startswith("http_404"):
            missing.append(acc)
        else:
            errors.append((acc, status))

        if done % 500 == 0 or done == len(todo):
            print(f"  {done}/{len(todo)} processed | ok={ok_count} missing={len(missing)} errors={len(errors)}")

print("\nDone.")
print(f"Downloaded this run: {ok_count}")
print(f"No AlphaFold model (404): {len(missing)}")
print(f"Other errors (network/timeout - safe to rerun script to retry): {len(errors)}")

with open("data/alphafold_missing.txt", "w") as f:
    f.write("\n".join(missing))
print("Missing accessions saved to data/alphafold_missing.txt")

if errors:
    print("\nSample errors:", errors[:10])
    print("Just rerun this script - it skips files already downloaded and will retry the rest.")
