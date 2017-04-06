import pandas as pd

mutations = pd.read_table('data/4dataset_nonsilent.txt',
                          names=['gene', 'patient', 'effect', 'categ'])
mat = pd.crosstab(mutations['patient'], mutations['gene'])
mat = mat.apply(lambda x: x > 0) * 1
mat.to_csv('results/gene_sample_crosstab.tsv', sep='\t')

# test

patients = mutations['patient'].values
genes = mutations['gene'].values
for patient, gene in zip(patients, genes):
    assert mat.loc[patient, gene] != 0
