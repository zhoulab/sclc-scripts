import argparse
import logging
import os
import csv
import sys
import pickle
import itertools

from logging import StreamHandler, Formatter


SAMPLES_FILE = 'data/samples.txt'
GENE_NONSILENT = 'data/genenonsilent200.txt'
PERC_DICT_FILE = 'data/PercDict.dic'
NUM_DICT_FILE = 'data/GeneNumDict.dic'

LISTS_DIR = 'GenePairsSig'


class GenePair(object):
    def __init__(self, gene1, gene2):
        self.Gene1 = gene1
        self.Gene2 = gene2

    def __getitem__(self, attr):
        return getattr(self, attr)

    def co_occurrence(self):
        pass

    def fisher(self, n):
        import scipy.stats as stats

        arr = ([n - self.Gene1Samples - self.Gene2Samples + self.Common,
                self.Gene2Samples - self.Common],
               [self.Gene1Samples - self.Common,
                self.Common])
        self.P_value = stats.fisher_exact(arr)[1]


class GenePairList(object):
    def __init__(self, pairs):
        self.pairs = list(pairs) if type(pairs) is not list else pairs

    def sort(self, attr, reverse=False):
        self.pairs = sorted(self.pairs, key=lambda (p): p[attr],
                            reverse=reverse)

    def run_fisher(self, n, log=None, adjust=False):
        if log:
            log.info('{: <100}|'.format('Progress:'))
        for i, pair in enumerate(self.pairs):
            if log and i % (len(self.pairs) / 100) is 0:
                sys.stdout.write('#')
            pair.fisher(n=n)
        if log:
            print
        if adjust:
            self.adjust_p()

    def adjust_p(self):
        from statsmodels.sandbox.stats.multicomp import multipletests

        p_vals = [pair.P_value for pair in self.pairs]
        adj_p_vals = multipletests(p_vals, method='fdr_bh')[1]
        for pair, adj_p in zip(self.pairs, adj_p_vals):
            pair.Adjusted_p = float(adj_p)

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


class Sample(object):
    def __init__(self, name):
        self.name = name
        self.genes = set()

    def has_gene(self, gene):
        return gene in self.genes


    def __getitem__(self, name):
        return self.samples[name]


def get_lines_from_file(fpath):
    """Return list of lines in a file without newline char"""
    with open(fpath) as file:
        return [line.rstrip() for line in file]


def generate_samples(fpath):
    """
    Return dict of Sample objects from `fpath`
    """
    with open(fpath) as file:
        return {line.strip(): Sample(line.strip()) for line in file}


def populate_sample_genes(fpath, sample_objs, sig_genes=None):
    """Add genes to sample based on `fpath`
    Assumed file columns: [gene, sample name, ...]

    Parameters
    ---------
    fpath : filepath
        line[0] : gene name
        line[1] : sample name
    sample_objs : dict of Sample objects
        {name: Sample(name)}
    sig_genes : list of str
        genes to filter for
    """
    with open(fpath) as file:
        reader = csv.reader(file, delimiter='\t')
        for line in reader:
            if sig_genes:
                if line[0] in sig_genes:
                    sample_objs[line[1]].genes.add(line[0])
            else:
                sample_objs[line[1]].genes.add(line[0])


def get_samples_with_gene(sample_objs, gene):
    """Yield iterable of sample names with `gene`

    Parameters
    ----------
    sample_objs : dict of Sample objects
    gene : str
    """
    for sample in sample_objs.values():
        if sample.has_gene(gene):
            yield sample.name


def get_common_samples(sample_objs, gene1, gene2):
    """Return set of samples which contain both `gene1` and `gene2`

    Parameters
    ----------
    sample_objs : dict of Sample objects
    gene1 : str
    gene2 : str
    """
    gene1_samples = set(get_samples_with_gene(sample_objs, gene1))
    gene2_samples = set(get_samples_with_gene(sample_objs, gene2))
    return gene1_samples & gene2_samples


def get_perc_dict(fpath):
    """Return dict from `fpath` using pickle.load()

    Convert percentages to floats,
    ignore the one case of {'Gene': 'Perc'}
    """
    with open(fpath) as f:
        return {k: float(v) for k, v in pickle.load(f).items()
                if k != 'Gene'}


def get_num_dict(fpath):
    """Return dict from `fpath` using pickle.load()

    Convert percentages to ints,
    ignore the one case of {'Gene': 'Perc'}
    """
    with open(fpath) as f:
        return {k: int(v) for k, v in pickle.load(f).items()
                if k != 'Gene'}


def get_pairs(sample_objs, perc_dict_file, num_dict_file,
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
    sample_objs : dict of Sample objects
    perc_dict_file : filepath
    num_dict_file : filepath
    filter_common : filter out common samples
    log : Logger object
    """
    log = log or logging.getLogger('dummy')
    perc_dict = get_perc_dict(perc_dict_file)
    num_dict = get_num_dict(num_dict_file)
    genes = sorted(set.union(*[sample.genes
                               for sample in sample_objs.values()]))
    num_pairs = len(list(itertools.combinations(genes, 2)))
    log.info('analyzing %i pairs', num_pairs)
    log.info('\n{: <100}|'.format('PROGRESS:'))
    for i, (gene1, gene2) in enumerate(itertools.combinations(genes, 2)):
        if log and i % (num_pairs / 100) is 0:
            sys.stdout.write('#')
        pair = GenePair(gene1, gene2)
        pair.Gene1Freq = perc_dict[gene1]
        pair.Gene2Freq = perc_dict[gene2]
        pair.Gene1Samples = num_dict[gene1]
        pair.Gene2Samples = num_dict[gene2]
        pair.Common = len(get_common_samples(sample_objs, gene1, gene2))
        if filter_common and pair.Common is 0:
            continue
        pair.PercofSamples = (100.0 * pair.Common / len(sample_objs))
        pair.Co_Occurrence = (100.0 * pair.PercofSamples /
                              (pair.Gene1Freq * pair.Gene2Freq))
        yield pair


def main():
    log = logging.getLogger()
    handler = StreamHandler(stream=sys.stdout)
    formatter = Formatter(fmt='%(asctime)-15s %(levelname)-6s %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--list',
                        help='gene list (Sig43List/Sig200List/MoreThan2)')
    # parser.add_argument('-a', '--alpha',
    #                     help='alpha level (0 to include all results)')
    parser.add_argument('--calc_out',
                        help='Output filename')
    parser.add_argument('--num_out',
                        help='Output filename')
    parser.add_argument('--filter_common', action='store_true',
                        help='Filter common')
    args = parser.parse_args()
    genes_file = args.list
    # alpha = float(args.alpha) if args.alpha else 0

    log.debug('generating sample objects')
    sample_objs = generate_samples(SAMPLES_FILE)

    log.debug('transforming GENE_NONSILENT, filter for sig genes')
    genes = get_lines_from_file(os.path.join(LISTS_DIR, genes_file + '.txt'))
    populate_sample_genes(GENE_NONSILENT, sample_objs, sig_genes=genes)
    log.debug('creating GenePairList')
    pairs_list = GenePairList(get_pairs(sample_objs, PERC_DICT_FILE,
                                        NUM_DICT_FILE, log=log,
                                        filter_common=args.filter_common))
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

if __name__ == '__main__':
    main()
