# gene-pairs

## Quickstart

    $ ./run_all.sh

## Pipeline

`run_pipeline.py` runs all sub-scripts collectively.

* `-i` mutation TSV file with the following columns: (same as one of the input files used for MutSig)
    * `gene`
    * `patient`
    * `effect`
    * `categ`
* `--mutsig_genes_file [GENE LIST]` **or** `-p [PERCENT THRESHOLD]`
* `-f` filter out pairs with 0 common samples

### 1. `gene_pairs.py`

#### Input data
* mutation TSV file
* either a percent threshold or MutSig genes output file defining the significant genes

#### Output
* GenePairsNum
* GenePairs

### 2. `gene_pairs_fisher.R`

#### Input Data
* GenePairsNum file (from `gene_pairs.py`)

#### Output
* Fisher test results (FisherGenePairs)

### 3. `co-occurrence_freq_plot.R`

#### Input Data
* GenePairs file (from `gene_pairs.py`)

#### Output
* co-occurrence frequency plot

## Output File Specifics

### GenePairs

* `GenePairs_{sig_gene_filter}.txt`

Columns:

* `Gene1`
* `Gene1Freq`
* `Gene2`
* `Gene2Freq`
* `PercofSamples`
* `Co_Occurrence` (sorted descending)

### GenePairsNum

* `GenePairsNum_{sig_gene_filter}.txt`
* flat file from `gene_pairs.py` to `gene_pairs_fisher.R`

Columns:

* `Gene1`
* `Gene1Samples`
* `Gene2`
* `Gene2Samples`
* `Common`

### FisherGenePairs

* `FisherGenePairs_{sig_gene_filter}.txt`

Columns:

* `Gene1`
* `Gene1Samples`
* `Gene2`
* `Gene2Samples`
* `Common`
* `P_value` (sort ascending)
* `Adjusted_p`

#### Process

For each gene pair, 4 values are computed:

* `N1 = 228 - (N2 + N3 + N4)`
* `N2 = Gene2Samples - CommonSamples`
* `N3 = Gene1Samples - CommonSamples`
* `N4 = CommonSamples`

A two-sided Fisher's exact test is performed on the matrix:

    N1 N2
    N3 N4

If the `-f` flag is applied, gene pairs are ignored from calculation and output if they share no common samples.
