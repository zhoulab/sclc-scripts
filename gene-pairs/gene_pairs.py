import argparse
import logging
import os
import csv
import sys
import warnings
import itertools

from logging import StreamHandler, Formatter


class GenePair(object):
    def __init__(self, gene1, gene2):
        self.Gene1 = gene1
        self.Gene2 = gene2

    def __getitem__(self, attr):
        return getattr(self, attr)

    def co_occurrence(self):
        pass


class GenePairList(object):
    def __init__(self, pairs):
        self.pairs = list(pairs) if type(pairs) is not list else pairs

    def sort(self, attr, reverse=False):
        self.pairs = sorted(self.pairs, key=lambda (p): p[attr],
                            reverse=reverse)

    def write(self, fpath, header, sort_attr=None, reverse=False, alpha=None):
        sorted_pairs = (sorted(self.pairs, key=lambda (p): p[sort_attr],
                               reverse=reverse)
                        if sort_attr else self.pairs)
        with open(fpath, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(header)
            for pair in sorted_pairs:
                if not alpha or pair.P_value < alpha:
                    writer.writerow(["{:.2E}".format(pair[col])
                                     if col in ['P_value', 'Adjusted_p']
                                     else pair[col]
                                     for col in header])


def get_mutsig_genes(mutsig_genes_file_file):
    """Return iterable of genes from MutSig genes output file"""
    with open(mutsig_genes_file_file) as f:
        reader = csv.reader(f, delimiter='\t')
        col = reader.next().index('gene')
        for row in reader:
            yield row[col]


def get_unique_samples(mutation_tsv_file):
    """Return set of patient IDs from a TSV file.

    Assume column 2 to be the patient column.
    """
    with open(mutation_tsv_file) as f:
        reader = csv.reader(f, delimiter='\t')
        patients = set([row[1] for row in reader])
    return patients


def get_common_samples(count_data, samples, gene1, gene2):
    """Return set of samples which contain both `gene1` and `gene2`"""
    for sample in samples:
        if count_data[gene1][sample] and count_data[gene2][sample]:
            yield sample


def get_count_data(mutation_tsv_file, outer_value):
    """Return count data for a mutation TSV file

    Input file must have:
    * no header
    * gene ID under column 1
    * patient ID under column 2

    Parameters
    ----------
    mutation_tsv_file : filepath
    outer_value : 'gene'/'patient'

    Return
    ------
    genes : set of gene IDs
    patients : set of patient IDs
    count_data :
    example, when outer_value='gene':
    {
        gene : {
                    patient : count
               }
    }
    """
    if outer_value not in ('gene', 'patient'):
        raise ValueError('outer_value must be "gene" or "patient"')
    with open(mutation_tsv_file) as f:
        reader = list(csv.reader(f, delimiter='\t'))
    genes = set([row[0] for row in reader])
    patients = set([row[1] for row in reader])
    if outer_value == 'gene':
        count_data = {gene: {patient: 0 for patient in patients}
                      for gene in genes}
        for row in reader:
            count_data[row[0]][row[1]] += 1
    else:  # outer_value='patient'
        count_data = {patient: {gene: 0 for gene in genes}
                      for patient in patients}
        for row in reader:
            count_data[row[1]][row[0]] += 1
    return genes, patients, count_data


def get_pairs(mutation_tsv_file, sig_genes=None, percent_threshold=None,
              filter_common=False, log=None):
    """
    Return GenePair objects with base attriutes:
    * Gene1Freq
    * Gene2Freq
    * Gene1Samples
    * Gene2Samples
    * Common
    * PercofSamples
    * Co_Occurrence

    Parameters
    ----------
    filter_common : filter out common samples
    log : Logger object
    """

    # get count data
    genes, patients, count_data = get_count_data(mutation_tsv_file,
                                                 outer_value='gene')
    # create dictionaries
    gene_counts = {gene: sum([int(bool(count))
                              for patient, count in count_data[gene].items()])
                   for gene in genes}
    gene_percentages = {gene: (100.0 * count / len(patients))
                        for gene, count in gene_counts.items()}

    if sig_genes:
        sig_genes = set(sig_genes)
        for gene in sig_genes.copy():
            if gene not in genes:
                warnings.warn('{} not in {}.'.format(gene, mutation_tsv_file))
                warnings.warn('Continuing without the significant gene.')
                sig_genes.remove(gene)
    if not sig_genes:
        sig_genes = [gene for gene in genes
                     if gene_percentages[gene] > percent_threshold]

    iter_pairs = list(itertools.combinations(sorted(sig_genes), 2))
    num_pairs = len(list(iter_pairs))
    log.info('analyzing %i pairs', num_pairs)
    log.info('\n{: <100}|'.format('PROGRESS:'))
    for i, (gene1, gene2) in enumerate(iter_pairs):
        if i % (num_pairs / 100) == 0:
            sys.stdout.write('#')
        pair = GenePair(gene1, gene2)
        pair.Gene1Freq = gene_percentages[gene1]
        pair.Gene2Freq = gene_percentages[gene2]
        pair.Gene1Samples = gene_counts[gene1]
        pair.Gene2Samples = gene_counts[gene2]
        pair.Common = 0
        for patient in patients:
            if count_data[gene1][patient] and count_data[gene2][patient]:
                pair.Common += 1
        if filter_common and pair.Common is 0:
            continue
        pair.PercofSamples = (100.0 * pair.Common / len(patients))
        pair.Co_Occurrence = (100.0 * pair.PercofSamples /
                              (pair.Gene1Freq * pair.Gene2Freq))
        yield pair


def write_patient_counts(mutation_tsv_file, filepath, ignore_duplicates=False):
    """Write patient count file to `filepath`

    Parameters
    ----------
    mutation_tsv_file : filepath
        see get_count_data()
    ignore_duplicates : bool
        if True, do not count duplicate instances of (patient, gene)
    """
    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        genes, patients, count_data = get_count_data(mutation_tsv_file,
                                                     outer_value='patient')
        writer.writerow(['patient', 'count'])
        for patient in patients:
            count = sum([int(bool(count)) if ignore_duplicates else count
                         for gene, count in count_data[patient].items()])
            writer.writerow([patient, count])


def main(args):
    log = logging.getLogger()
    handler = StreamHandler(stream=sys.stdout)
    formatter = Formatter(fmt='%(asctime)-15s %(levelname)-6s %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.DEBUG)

    if not args:
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--mutation_tsv_file', required=True,
                            help='Mutation TSV file with columns gene/patient/effect/categ')
        parser.add_argument('--mutsig_genes_file',
                            help='MutSig genes output file')
        parser.add_argument('-p', '--percent_threshold', type=int,
                            help='threshold for patients count / total count')
        parser.add_argument('--calc_out', required=True,
                            help='Output filename')
        parser.add_argument('--num_out', required=True,
                            help='Output filename')
        parser.add_argument('--patient_count_out', required=True,
                            help='Patient count output filepath')
        parser.add_argument('--patient_count_dup_out', required=True,
                            help='Patient count output with duplicate genes')
        parser.add_argument('--filter_common', action='store_true',
                            help='Filter common')
        args = parser.parse_args()

    if args.mutsig_genes_file and args.percent_threshold:
        raise Exception('Cannot accept both a MutSig output file '
                        'and percent threshold.')
    if not any([args.mutsig_genes_file, args.percent_threshold]):
        raise Exception('Required: either a MutSig output file '
                        'or percent threshold.')
    elif args.mutsig_genes_file:
        log.info('Filtering for sig genes from %s...', args.mutsig_genes_file)
        sig_genes = get_mutsig_genes(args.mutsig_genes_file)
        pairs = get_pairs(args.mutation_tsv_file,
                          sig_genes=sig_genes, log=log,
                          filter_common=args.filter_common)
    elif args.percent_threshold:
        log.info('Filtering for sig genes using threshold=%i...',
                 args.percent_threshold)
        pairs = get_pairs(args.mutation_tsv_file,
                          percent_threshold=args.percent_threshold, log=log,
                          filter_common=args.filter_common)
    log.debug('creating GenePairList')
    pairs_list = GenePairList(pairs)
    print

    if args.calc_out:
        log.debug('writing file %s', args.calc_out)
        pairs_list.write(args.calc_out,
                         header=['Gene1', 'Gene1Freq', 'Gene2', 'Gene2Freq',
                                 'PercofSamples', 'Co_Occurrence'],
                         sort_attr='Co_Occurrence', reverse=True)
    if args.num_out:
        log.debug('writing file %s', args.num_out)
        pairs_list.write(args.num_out,
                         header=['Gene1', 'Gene1Samples',
                                 'Gene2', 'Gene2Samples', 'Common'])
    write_patient_counts(args.mutation_tsv_file, args.patient_count_out,
                         ignore_duplicates=True)
    write_patient_counts(args.mutation_tsv_file, args.patient_count_dup_out,
                         ignore_duplicates=False)

if __name__ == '__main__':
    main()
