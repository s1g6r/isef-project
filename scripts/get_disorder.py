import pandas as pd
import urllib.request
import gzip

# Download MobiDB disorder annotations for human proteome
# These are precomputed per-residue disorder calls - no AlphaFold needed
url = "https://mobidb.org/api/download?proteome=UP000005640&projection=mobidb_lite_disorder_merged&format=tsv"

print("Downloading MobiDB disorder annotations for human proteome...")
print("This may take a few minutes...")

try:
    urllib.request.urlretrieve(url, 'data/mobidb_disorder.tsv')
    print("Downloaded successfully")
    
    df = pd.read_csv('data/mobidb_disorder.tsv', sep='\t')
    print(f"Rows: {len(df)}")
    print("Columns:", df.columns.tolist())
    print(df.head())
    
except Exception as e:
    print(f"Download failed: {e}")
    print("Try the manual fallback below:")
    print("Go to: https://mobidb.org/proteome/UP000005640")
    print("Download the disorder annotations manually")
