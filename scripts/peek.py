import pandas as pd
cv = pd.read_csv('data/variant_summary.txt.gz', sep='\t', nrows=1000)
print(cv.columns.tolist())
print(cv[['GeneSymbol','ClinicalSignificance','Assembly','Chromosome','Start','ReferenceAlleleVCF','AlternateAlleleVCF']].head(10))
