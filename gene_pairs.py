import argparse
import logging
import os
import csv
import sys
import pickle
import itertools
import scipy.stats as stats

from statsmodels.sandbox.stats.multicomp import multipletests

SAMPLES_FILE = 'data/samples.txt'
GENE_NONSILENT = 'data/genenonsilent200.txt'
PERC_DICT_FILE = 'data/PercDict.dic'
NUM_DICT_FILE = 'data/GeneNumDict.dic'

LISTS_DIR = 'GenePairsSig'
CALC_OUT = 'results/GenePairs_{}.txt'
FISHER_OUT = 'results/FisherGenePairs_{}.txt'


class GenePair(object):
    def __init__(self, gene1, gene2):
        self.Gene1 = gene1
        self.Gene2 = gene2

    def __getitem__(self, attr):
        return getattr(self, attr)

    def co_occurrence(self):
        pass

    def fisher(self, n):
        arr = ([n - self.Gene1Samples - self.Gene2Samples + self.Common,
                self.Gene2Samples - self.Common],
               [self.Gene1Samples - self.Common,
                self.Common])
        self.P_value = stats.fisher_exact(arr)[1]


class GenePairList(object):
    def __init__(self, pairs):
        self.pairs = pairs

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


def generate_sample_ids(fpath):
    """
    Return dict of {sample_name: id} values
    Use line number for id (autoincrement starting at 1)
    """
    with open(fpath) as file:
        return {line.strip(): (i + 1) for i, line in enumerate(file)}


def get_sample_ids(fpath):
    """
    Return dict of {sample_name: id} values from a file
    File should be in a 2-col format (id, sample_id).
    """
    with open(fpath) as file:
        reader = csv.reader(file, delimiter='\t')
        return {line[1]: line[0] for line in reader}


def write_sample_ids(sample_ids, out_file):
    """Write sample name, id to out_file sorted by sample name"""
    with open(out_file, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        for s, s_id in sorted(sample_ids.items()):
            writer.writerow([s, s_id])


def get_transform_nonsilent(fpath, sample_ids, sig_genes=None):
    """Yield (gene, sample id) from file
    Assumed file columns: [gene, sample name, ...]
    """
    with open(fpath) as file:
        reader = csv.reader(file, delimiter='\t')
        for line in reader:
            if sig_genes:
                if line[0] in sig_genes:
                    yield (line[0], sample_ids[line[1]])
            else:
                yield (line[0], sample_ids[line[1]])


def write_transform_nonsilent(fpath, data):
    """Write data in two-column format (gene, sample id)
    Parameters
    ----------
    data: list of iterables
    """
    with open(fpath, 'w') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        for row in data:
            writer.writerow(row)


def get_lines_from_file(fpath):
    """Return list of lines in a file without newline char"""
    with open(fpath) as file:
        return [line.rstrip() for line in file]


def get_samples_with_gene(data, gene):
    """Yield iterable of samples with `gene`

    Parameters
    ----------
    data : list of 2-tuples (gene, sample id)
    """
    for g, s_id in data:
        if g == gene:
            yield s_id


def get_common_samples_between(data, gene1, gene2):
    """Return set of samples that contain both `gene1` and `gene2`

    Parameters
    ----------
    data : list of 2-tuples (gene, sample id)
    gene1 : str
    gene2 : str
    """
    gene1_samples = set(get_samples_with_gene(data, gene1))
    gene2_samples = set(get_samples_with_gene(data, gene2))
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


def get_pairs(data, perc_dict_file, num_dict_file,
              filter_common=True, log=None):
    """
    Return gene pairs with base attributes

    Parameters
    ----------
    data : list of 2-tuples (gene, sample id)
    """
    perc_dict = get_perc_dict(perc_dict_file)
    num_dict = get_num_dict(num_dict_file)
    genes = sorted(set([g for g, s_id in data]))
    ids = sorted(set([s_id for g, s_id in data]))

    if log:
        log.info('{: <100}|'.format('Progress:'))
    for i, (gene1, gene2) in enumerate(itertools.combinations(genes, 2)):
        if log and i % (len(list(itertools.combinations(genes, 2))) / 100) is 0:
            sys.stdout.write('#')
        pair = GenePair(gene1, gene2)
        pair.Gene1Freq = perc_dict[gene1]
        pair.Gene2Freq = perc_dict[gene2]
        pair.Gene1Samples = num_dict[gene1]
        pair.Gene2Samples = num_dict[gene2]
        pair.Common = len(get_common_samples_between(data, gene1, gene2))
        if filter_common and pair.Common is 0:
            continue
        pair.PercofSamples = (100.0 * pair.Common / len(ids))
        pair.Co_Occurrence = (100.0 * pair.PercofSamples /
                              (pair.Gene1Freq * pair.Gene2Freq))
        yield pair


def write_files(pairs_list, gene_pairs_file, fisher_file, alpha=None):
    """
    Parameters
    ----------
    pairs_list : GenePairList
    """

    pairs_list.write(gene_pairs_file,
                     header=['Gene1', 'Gene1Freq', 'Gene2', 'Gene2Freq',
                             'PercofSamples', 'Co_Occurrence'],
                     sort_attr='Co_Occurrence', reverse=True)

    pairs_list.write(fisher_file,
                     header=['Gene1', 'Gene1Samples', 'Gene2', 'Gene2Samples',
                             'Common', 'P_value', 'Adjusted_p'],
                     sort_attr='P_value', alpha=alpha)


if __name__ == '__main__':
    log = logging.getLogger()
    handler = logging.StreamHandler(stream=sys.stdout)
    log.addHandler(handler)
    log.setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--list',
                        help='gene list (Sig43List/Sig200List/MoreThan2)')
    parser.add_argument('-a', '--alpha',
                        help='alpha level (0 to include all results)')
    args = parser.parse_args()
    genes_file = args.list
    alpha = float(args.alpha) if args.alpha else 0

    calc_fpath = CALC_OUT.format(genes_file)
    fname_append = ('{}_{}'.format(genes_file, str(alpha))
                    if alpha else genes_file)
    fisher_fpath = FISHER_OUT.format(fname_append)

    log.debug('generating sample IDs')
    sample_ids = generate_sample_ids(SAMPLES_FILE)

    log.debug('transforming GENE_NONSILENT using sample IDs')
    transform_data = get_transform_nonsilent(GENE_NONSILENT, sample_ids)

    log.debug('transforming GENE_NONSILENT using sample IDs, filter for sig genes')
    genes = get_lines_from_file(os.path.join(LISTS_DIR, genes_file + '.txt'))
    genes_samples = get_transform_nonsilent(GENE_NONSILENT, sample_ids,
                                            sig_genes=genes)
    gs_set = set(genes_samples)
    log.debug('creating GenePairList')
    pairs_list = GenePairList(list(get_pairs(gs_set, PERC_DICT_FILE, NUM_DICT_FILE,
                                             log=log)))
    print

    log.debug('running Fisher tests')
    pairs_list.run_fisher(n=228, log=log, adjust=True)

    log.debug('writing files')
    write_files(pairs_list, calc_fpath, fisher_fpath, alpha=alpha)
    log.debug('Done. Result files:\n{}\n{}'.format(calc_fpath, fisher_fpath))
