import pandas as pd

# read mutation dataset
mutations = pd.read_table('data/4dataset_nonsilent.txt',
                          names=['gene', 'patient', 'effect', 'categ'])
# read set of significant genes
genes_df = pd.read_table('data/SigGenes_005.txt', usecols=['gene'])
genes = set(genes_df.gene.values)

# subset mutation dataset for significant genes
sig_gene_mutations = mutations[mutations['gene'].isin(genes)]

# get rb1 samples
rb1_sig_gene_mutations = sig_gene_mutations[sig_gene_mutations['gene'] == 'RB1']
rb1_samples = set(rb1_sig_gene_mutations['patient'].values)
all_samples = set(sig_gene_mutations['patient'].values)
print len(rb1_samples), 'out of', len(all_samples), 'samples have RB1 mutation'

rb1_samples_sig_gene_mutations = sig_gene_mutations[sig_gene_mutations['patient'].isin(rb1_samples)]


# write significant gene mutations from RB1 mutated samples
rb1_samples_sig_gene_mutations.to_csv('results/SigGenes_005_RB1_samples_4dataset_nonsilent.tsv',
                                      sep='\t', index=False)

# write significant gene mutations from non-RB1 mutated samples
no_rb1_samples_sig_gene_mutations = sig_gene_mutations[~sig_gene_mutations['patient'].isin(rb1_samples)]
no_rb1_samples_sig_gene_mutations.to_csv('results/SigGenes_005_not_RB1_samples_4dataset_nonsilent.tsv',
                                         sep='\t', index=False)
