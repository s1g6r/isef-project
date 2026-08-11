import pandas as pd

print("Loading AlphaMissense sample...")

# Try standard path — adjust filename if yours differs slightly
try:
    am = pd.read_csv('data/AlphaMissense_hg38.tsv.gz', sep='\t', nrows=1000, comment='#')
except FileNotFoundError:
    # Try without .gz
    am = pd.read_csv('data/AlphaMissense_hg38.tsv', sep='\t', nrows=1000, comment='#')

print("Columns:", am.columns.tolist())
print("\nFirst 5 rows:")
print(am.head())
print("\nChromosome format examples:", am.iloc[:, 0].unique()[:5])  # first col is usually chrom

