import pandas as pd

# Load full file this time (no nrows limit)
print("Loading ClinVar...")
cv = pd.read_csv('data/variant_summary.txt.gz', sep='\t', low_memory=False)
print(f"Total rows: {len(cv)}")

# Keep only GRCh38
cv38 = cv[cv['Assembly'] == 'GRCh38'].copy()
print(f"GRCh38 rows: {len(cv38)}")

# Keep only germline missense-relevant classifications
keep = ['Pathogenic', 'Likely pathogenic', 'Benign', 'Likely benign',
        'Pathogenic/Likely pathogenic', 'Benign/Likely benign']
cv38 = cv38[cv38['ClinicalSignificance'].isin(keep)].copy()
print(f"After classification filter: {len(cv38)}")

# Simplify labels to binary: 1 = pathogenic, 0 = benign
path_map = {
    'Pathogenic': 1,
    'Likely pathogenic': 1,
    'Pathogenic/Likely pathogenic': 1,
    'Benign': 0,
    'Likely benign': 0,
    'Benign/Likely benign': 0
}
cv38['label'] = cv38['ClinicalSignificance'].map(path_map)

# Parse date
cv38['LastEvaluated'] = pd.to_datetime(cv38['LastEvaluated'], errors='coerce')
print(f"Rows with valid date: {cv38['LastEvaluated'].notna().sum()}")

# Keep columns we actually need
cols = ['#AlleleID', 'GeneSymbol', 'Chromosome', 'PositionVCF',
        'ReferenceAlleleVCF', 'AlternateAlleleVCF',
        'label', 'ClinicalSignificance', 'LastEvaluated',
        'ReviewStatus', 'NumberSubmitters']
cv38 = cv38[cols]

# Save clean version
cv38.to_csv('data/clinvar_grch38_clean.csv', index=False)
print("Saved to data/clinvar_grch38_clean.csv")
print(cv38.head())
