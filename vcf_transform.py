import pandas as pd
import numpy as np

from collections import OrderedDict
from argparse import ArgumentParser


def get_vcf_df(filename):
    """Return pandas.DataFrame of VCF file.

    * skip lines starting with '##'
    * remove '#' from header line
    """
    with open(filename) as f:
        # skip initial comment lines
        skip_count = 0
        line = f.readline()
        while line[:2] == '##':
            skip_count += 1
            line = f.readline()
        header = line.rstrip().split('\t')
        if line[0] == '#':
            header[0] = header[0][1:]
        skip_count += 1
    return pd.read_table(filename, header=None,
                         names=header, skiprows=skip_count)


def transform_vcf_df(vcf_df):
    """Return pandas.DataFrame of transformed VCF file.

    * replace './.:.:.:.:.:.' with NaN values for mutation

    Variables
    ---------
    transform_columns : OD
        column definitions in the form
        { transformed column : intiial column }
    """
    transform_columns = OrderedDict([('chr', 'CHROM'),
                                     ('start', 'POS'),
                                     ('end', 'POS'),
                                     ('ref_allele', 'REF'),
                                     ('alt_allele', 'ALT'),
                                     ('patient', 'PATIENT')])
    id_cols = map(str, vcf_df.columns[:9])
    transform_df = pd.melt(vcf_df, id_vars=id_cols,
                           var_name='PATIENT',
                           value_name='MUTATION')
    transform_df['MUTATION'].replace('./.:.:.:.:.:.', np.nan,
                                     inplace=True)
    transform_df = transform_df.loc[transform_df['MUTATION'].notnull(),
                                    transform_columns.values()]
    transform_df.columns = transform_columns.keys()
    return transform_df


if __name__ == '__main__':
    parser = ArgumentParser(description='Transform a VCF file.')
    parser.add_argument('filename', help='input filename')
    parser.add_argument('-o', '--output',
                        help='output filename (default: input filename)')
    args = parser.parse_args()
    if not args.output:
        args.output = args.filename
    vcf_df = get_vcf_df(args.filename)
    transform_df = transform_vcf_df(vcf_df)
    transform_df.to_csv(args.output, sep='\t', index=False)
