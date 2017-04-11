# Comutation Plot

More information: `python genes_001.py --help`

## Data

* `--mutsig_genes_file`
    * MutSig output file (TSV) containing information on significant genes
    * Used columns: `gene`, `p`
* `--mutation_tsv_file`
    * TSV file containing gene mutation information
    * Columns: `gene`, `patient`, `effect`, `categ`

## Mutation type categories

Definitions for 4th column in `--mutation_tsv_file`

* `3`: C:G transitions
* `4`: C:G transversions
* `5`: A:T transitions
* `6`: A:T transversions
* `7`: null+indel mutations

## Filtering

1. `--p_value` and `--gene_list_file` are disjointly two primary filters
and sort methods
2. `--num_genes` is applied as an additional cutoff
