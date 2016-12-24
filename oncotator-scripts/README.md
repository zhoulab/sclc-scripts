# oncotator files

Full pipeline description: https://github.com/zhoulab/sclc-scripts/wiki/SCLC-mutation-analysis-pipeline

## `vcf_transform.py`

### `VCFFile`

An extension of `pyvcf.Reader` that provides data access through `pandas.DataFrame`.

## `oncotator_maf_cross_check.py`

Script to validate consistency between two oncotator MAF output
files.

## `pandastsv.py`

### `get_tsv_df()`

helper function used to load files in `oncotator_maf_cross_check.py`.