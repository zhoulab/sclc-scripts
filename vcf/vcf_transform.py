#!/usr/bin/env python

import pandas as pd
import numpy as np

from collections import OrderedDict
from argparse import ArgumentParser

from vcf import Reader


class VCFFile(Reader):
    """Extension of PyVCF's Reader with pandas.DataFrame.


    Additional attributes
    ---------------------
    see PyVCF documentation for Reader attributes.
    df : pandas.DataFrame


    Also handles no-call format strings within the DataFrame.

    """

    def __init__(self, filename, **kwargs):
        """Must take filename for pandas.DataFrame loading."""
        Reader.__init__(self, filename=filename, **kwargs)
        self._init_df(filename)

    def _init_df(self, filename):
        """Initialize pandas.DataFrame of VCF file."""
        self.df = pd.read_table(filename, header=None,
                                names=(self._column_headers + self.samples),
                                skiprows=(len(self._header_lines) + 1))

    def get_maflite_df(self):
        """Return melted pandas.DataFrame of VCF with MAFLITE columns.

        * replace no-call format strings with NaN values.

        Variables
        ---------
        maflite_columns : OrderedDict
            column definitions in the form
            { transformed column : intiial column }
        """
        maflite_columns = OrderedDict([('chr', 'CHROM'),
                                       ('start', 'POS'),
                                       ('end', 'POS'),
                                       ('ref_allele', 'REF'),
                                       ('alt_allele', 'ALT'),
                                       ('patient', 'SAMPLE')])
        maflite_df = pd.melt(self.df, id_vars=self._column_headers,
                             var_name='SAMPLE',
                             value_name='VARIANT')
        for no_call_string in self._get_no_call_strings():
            maflite_df['VARIANT'].replace(no_call_string, np.nan,
                                          inplace=True)
        maflite_df = maflite_df.loc[maflite_df['VARIANT'].notnull(),
                                    maflite_columns.values()]
        maflite_df.columns = maflite_columns.keys()
        return maflite_df

    def _get_no_call_strings(self):
        """Return no-call string for each format string in the format cache."""
        # generate the cache
        for x in self.next():
            pass
        for format_string in self._format_cache:
            yield ':'.join([self._get_no_call_value(self.formats[x])
                            for x in format_string.split(':')])

    def _get_no_call_value(self, format_namedtuple):
        """Return no-call string for a vcf.parser._Format namedtuple.

        default: .

        Specific formats
        ----------------
        GT:  ./.

        """
        if format_namedtuple.id == 'GT':
            return './.'
        return '.'


if __name__ == '__main__':
    parser = ArgumentParser(description='Generate MAFLITE file from VCF.')
    parser.add_argument('filename', help='input VCF file')
    parser.add_argument('output', help='output MAFLITE file')
    args = parser.parse_args()
    vcf_file = VCFFile(args.filename)
    maflite_df = vcf_file.get_maflite_df()
    maflite_df.to_csv(args.output, sep='\t', index=False)
