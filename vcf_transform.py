from pandas import read_table


def get_vcf_df(fname):
    """Return pandas.DataFrame of VCF file."""
    with open(fname) as f:
        # skip initial comment lines
        skip_ct = 0
        line = f.readline()
        while line[:2] == '##':
            skip_ct += 1
            line = f.readline()
        header = line.split('\t')
        if line[0] == '#':
            header[0] = header[0][1:]
        skip_ct += 1
        return read_table(fname, header=None, names=header, skiprows=skip_ct)


if __name__ == '__main__':
    vcf_df = get_vcf_df('SCLC.muTect.vcf')
    # TODO: tranformation
