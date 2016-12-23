#!/usr/bin/env python

import pandas as pd


def get_tsv_df(filename, **kwargs):
    """Return pandas.DataFrame of file.

    * skip lines starting with '##'
    * remove '#' from header line
    * contents should be tab-separated
    """
    with open(filename) as f:
        # skip initial comment lines
        skip_count = 0
        line = f.readline()
        while '\t' not in line:
            skip_count += 1
            line = f.readline()
        header = line.rstrip().split('\t')
        # VCF header rows start with '#'
        if line[0] == '#':
            header[0] = header[0][1:]
        skip_count += 1
    return pd.read_table(filename, header=None, names=header,
                         skiprows=skip_count, **kwargs)
