import pandas as pd
from sklearn.metrics import roc_auc_score

print("Loading AlphaMissense...")
am = pd.read_csv(
    'data/AlphaMissense_hg38.tsv.gz',
    sep='\t',
    comment='#',
    header=None,
    names=['chrom','pos','ref','alt','genome','uniprot','transcript','aa_change','am_score','am_class']
)
print(f"AlphaMissense rows: {len(am)}")
print(am.head(3))

# Fix chromosome format: 'chr1' -> '1'
am['chrom'] = am['chrom'].str.replace('chr', '', regex=False)

# Make score numeric
am['am_score'] = pd.to_numeric(am['am_score'], errors='coerce')

print("\nLoading ClinVar clean...")
cv = pd.read_csv('data/clinvar_grch38_clean.csv', dtype={'Chromosome': str})

# Build join keys
cv['chrom']  = cv['Chromosome'].astype(str)
cv['pos']    = cv['PositionVCF'].astype(str)
cv['ref']    = cv['ReferenceAlleleVCF'].astype(str)
cv['alt']    = cv['AlternateAlleleVCF'].astype(str)

am['pos'] = am['pos'].astype(str)
am['ref'] = am['ref'].astype(str)
am['alt'] = am['alt'].astype(str)

print("\nJoining...")
merged = cv.merge(am[['chrom','pos','ref','alt','am_score','am_class']],
                  on=['chrom','pos','ref','alt'],
                  how='inner')
print(f"Variants after join: {len(merged)}")
print(f"Label distribution:\n{merged['label'].value_counts()}")

# Sanity check: reproduce known AUC
valid = merged.dropna(subset=['am_score','label'])
auc = roc_auc_score(valid['label'], valid['am_score'])
print(f"\nOverall AUC (AlphaMissense vs ClinVar labels): {auc:.4f}")
print("Expected: ~0.90-0.93 based on published numbers")
print("\nSample of merged data:")
print(merged[['GeneSymbol','chrom','pos','ref','alt','label','am_score','am_class','LastEvaluated']].head(10))

merged.to_csv('data/clinvar_am_joined.csv', index=False)
print("\nSaved to data/clinvar_am_joined.csv")
