# gene-pairs

## Quickstart

### HPC

`cd /ufrc/zhou/share/projects/bioinformatics/gene-pairs`

    $ . ve/bin/activate
    $ python -u gene_pairs.py --list [GENE LIST] --alpha [ALPHA LVL]

### Local

    $ ./bootstrap.sh
    $ . ve/bin/activate
    $ python -u gene_pairs.py --list [GENE LIST] --alpha [ALPHA LVL]

## Requirements

* `scipy`
* `statsmodels` (0.8.0rc1)

## Data Source

Stored in `data/` directory

* `genenonsilent200.txt`
* `samples.txt`
* `GeneNumDict.dic`
* `PercDict.dic`
* `Sig43List.txt`/`Sig200List.txt`/`MoreThan2.txt`

## User Input

* gene list
	* `Sig43List`/`Sig200List`/`MoreThan2`
* alpha level
	* `0` to show results

## Output files

Stored in `results/` directory.

### Gene Pairs

* `GenePairs_{list_name}.txt`

Columns:

* `Gene1`
* `Gene1Freq`
* `Gene2`
* `Gene2Freq`
* `PercofSamples`
* `Co_Occurrence` (sort descending)

### Gene Pairs with Fisher Test

* `FisherGenePairs_{list_name}_{alpha}.txt` or
`FisherGenePairs_{list_name}.txt`

Columns:

* `Gene1`
* `Gene1Samples`
* `Gene2`
* `Gene2Samples`
* `Common`
* `P_value` (sort ascending)
* `Adjusted_p`

## Process

For each gene pair, 4 values are computed:

* `N1 = 228 - (N2 + N3 + N4)`
* `N2 = Gene2Samples - CommonSamples`
* `N3 = Gene1Samples - CommonSamples`
* `N4 = CommonSamples`

A two-sided Fisher's exact test is performed on the matrix:

    N1 N2
    N3 N4

Gene pairs are ignored from calculation and output if they share no common samples.

All gene pairs with common samples are written to the first file, sorted by `Co-Occurrence`.

With the user-defined alpha level, significant rows are written to the Fisher file (based on `P_value`). If no alpha level is specified, all gene pairs with common samples will be written.
