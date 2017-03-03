import pandas as pd

# read mutation dataset
mutations = pd.read_table('data/4dataset_nonsilent.txt',
                          names=['gene', 'patient', 'effect', 'categ'])
# read set of significant genes
genes_df = pd.read_table('data/SigGenes_005.txt', usecols=['gene'])
genes = set(genes_df.gene.values)

# subset mutation dataset for significant genes
sig_gene_mutations = mutations[mutations['gene'].isin(genes)]

# write RB1 mutations
rb1_sig_gene_mutations = sig_gene_mutations[sig_gene_mutations['gene'] == 'RB1']
rb1_sig_gene_mutations.to_csv('results/RB1_4dataset_nonsilent.tsv',
                              sep='\t', index=False)

# write non-RB1 significant gene mutations
not_rb1_sig_gene_mutations = sig_gene_mutations[sig_gene_mutations['gene'] != 'RB1']
not_rb1_sig_gene_mutations.to_csv('results/SigGenes_005_notRB1_4dataset_nonsilent.tsv',
                                  sep='\t', index=False)
