"""
Script to check if chr/start/end matches gene name
for different MAF outputs of oncotator.

*** MAY REQUIRE HIGH MEMORY USAGE ***

FILE1
    filepath of first file
FILE2
    filepath of second file
KEY_COLS
    list of column names to use as unique key
COL_TO_MATCH
    column name to cross check
"""

from pandastsv import get_tsv_df


FILE1 = 'transformed_SCLC_onco.tsv'
FILE2 = 'oncotator_muTect_122116VS.maf'
KEY_COLS = ['Chromosome', 'Start_position', 'End_position']
COL_TO_MATCH = 'Entrez_Gene_Id'

# get DataFrames
df1 = get_tsv_df(FILE1, low_memory=False)
df2 = get_tsv_df(FILE2, low_memory=False)

# assert same number of rows
assert len(df1) == len(df2), 'Number of rows are different.'

# group by key columns
df1_groups = df1.groupby(KEY_COLS)
df2_groups = df2.groupby(KEY_COLS)

for key, df1_group in df1_groups:
    df2_group = df2_groups.get_group(key)
    val1 = df1_group[COL_TO_MATCH].values[0]
    val2 = df2_group[COL_TO_MATCH].values[0]
    # assert on the specified column
    assert (val1 == val2), ', '.join(['{}={}'.format(x, eval(x))
                                      for x in ['key', 'val1', 'val2']])
print 'These two files have identical information.'
