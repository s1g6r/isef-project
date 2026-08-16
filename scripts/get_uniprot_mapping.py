import time
import urllib.request

import pandas as pd

# Session 11 continued: ProteinGym's "UniProt_ID" reference column turns out
# to be UniProt's mnemonic entry name (e.g. "PAI1_HUMAN"), not the accession
# number (e.g. "P05121") that AlphaMissense's genome-wide file uses as its
# own join key -- confirmed by checking one known example directly against
# UniProt's own REST API rather than assuming the column meant what its name
# implied. This maps the 96 human ProteinGym entry names to their real
# accessions the same way, one REST lookup per ID -- only 96 IDs, so this is
# a quick targeted pass, not a bulk download.

REFERENCE_PATH = "data/proteingym/DMS_substitutions_reference.csv"
OUT_PATH = "data/proteingym/uniprot_entry_name_to_accession.csv"


def lookup_accession(entry_name):
    url = f"https://rest.uniprot.org/uniprotkb/{entry_name}.tsv?fields=accession"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            lines = resp.read().decode().strip().split("\n")
            if len(lines) < 2:
                return None
            # UniProt's .tsv response returns more columns than the
            # requested "fields=accession" alone (entry name, description,
            # etc., all tab-separated on one line) -- the accession is
            # always the first field, so split rather than trust the whole
            # line is just the one value asked for.
            return lines[1].split("\t")[0].strip()
    except Exception as e:
        print(f"  FAILED for {entry_name}: {e}")
        return None


def main():
    ref = pd.read_csv(REFERENCE_PATH)
    human_ref = ref[ref["taxon"].str.lower() == "human"]
    entry_names = sorted(human_ref["UniProt_ID"].unique())
    print(f"Looking up {len(entry_names)} human UniProt entry names...")

    rows = []
    for i, name in enumerate(entry_names, 1):
        accession = lookup_accession(name)
        rows.append({"entry_name": name, "accession": accession})
        if i % 10 == 0:
            print(f"  ...{i}/{len(entry_names)} looked up")
        time.sleep(0.1)

    mapping = pd.DataFrame(rows)
    n_failed = mapping["accession"].isna().sum()
    print(f"\n{len(mapping) - n_failed}/{len(mapping)} resolved successfully.")
    if n_failed:
        print("Failed entries (won't be usable for the AlphaMissense join):")
        print(mapping[mapping["accession"].isna()])

    mapping.to_csv(OUT_PATH, index=False)
    print(f"Saved {OUT_PATH}")


if __name__ == "__main__":
    main()
