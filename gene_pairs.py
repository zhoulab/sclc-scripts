import argparse
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
OUT_DIR = 'results'
CALC_OUT = os.path.join(OUT_DIR, 'GenePairs_{}.txt')
FISHER_OUT = os.path.join(OUT_DIR, 'FisherGenePairs_{}.txt')


class GenePair(object):
    def __init__(self, gene1, gene2):
        self.Gene1 = gene1
        self.Gene2 = gene2

    def __getitem__(self, attr):
        return getattr(self, attr)


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


def write_files(data, perc_dict_file, num_dict_file,
                gene_pairs_file, fisher_file, alpha=None):
    """
    Write final gene pairs file with header:
    Gene1, Gene1Freq, Gene2, Gene2Freq, PercofSamples, Co-Occurrence

    Parameters
    ----------
    data : list of 2-tuples (gene, sample id)
    """
    perc_dict = get_perc_dict(perc_dict_file)
    num_dict = get_num_dict(num_dict_file)
    genes = sorted(set([g for g, s_id in data]))
    ids = sorted(set([s_id for g, s_id in data]))

    pairs = set()
    print '{: <100}|'.format('Progress:')
    num_combinations = len(list(itertools.combinations(genes, 2)))
    for i, (gene1, gene2) in enumerate(itertools.combinations(genes, 2)):
        if i % (num_combinations / 100) is 0:
            sys.stdout.write('#')
        pair = GenePair(gene1, gene2)
        pair.Gene1Freq = perc_dict[gene1]
        pair.Gene2Freq = perc_dict[gene2]
        pair.Gene1Samples = num_dict[gene1]
        pair.Gene2Samples = num_dict[gene2]
        pair.Common = len(get_common_samples_between(data, gene1, gene2))
        if pair.Common is 0:
            continue
        pair.PercofSamples = 100.0 * pair.Common / len(ids)
        pair.Co_Occurrence = 100.0 * pair.PercofSamples / (pair.Gene1Freq * pair.Gene2Freq)

        arr = ([228 - pair.Gene1Samples - pair.Gene2Samples + pair.Common,
                pair.Gene2Samples - pair.Common],
               [pair.Gene1Samples - pair.Common,
                pair.Common])
        __, p = stats.fisher_exact(arr)
        pair.P_value = float(p)

        pairs.add(pair)
    print

    with open(gene_pairs_file, 'w') as f:
        header = ['Gene1', 'Gene1Freq', 'Gene2', 'Gene2Freq',
                  'PercofSamples', 'Co_Occurrence']
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
        sorted_pairs = sorted(pairs, key=lambda (p): p.Co_Occurrence, reverse=True)
        for pair in sorted_pairs:
            writer.writerow([pair[col] for col in header])

    with open(fisher_file, 'w') as f:
        header = ['Gene1', 'Gene1Samples', 'Gene2', 'Gene2Samples',
                  'Common', 'P_value', 'Adjusted_p']
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
        sorted_pairs = sorted(pairs, key=lambda (p): p.P_value)
        p_vals = [pair.P_value for pair in sorted_pairs]
        adj_p_vals = multipletests(p_vals, method='fdr_bh')[1]
        for pair, adj_p in zip(sorted_pairs, adj_p_vals):
            pair.Adjusted_p = float(adj_p)
            if not alpha or pair.P_value < alpha:
                writer.writerow(["{:.2E}".format(pair[col])
                                 if type(pair[col]) is float else pair[col]
                                 for col in header])


if __name__ == '__main__':
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

    # generate sample IDs
    sample_ids = generate_sample_ids(SAMPLES_FILE)

    # transform GENE_NONSILENT using sample IDs
    transform_data = get_transform_nonsilent(GENE_NONSILENT, sample_ids)

    # transform GENE_NONSILENT using sample IDs, filter for sig genes
    sig_genes = get_lines_from_file(os.path.join(LISTS_DIR, genes_file + '.txt'))
    genes = get_transform_nonsilent(GENE_NONSILENT, sample_ids,
                                    sig_genes=sig_genes)
    two_perc_genes_set = sorted(set(genes), key=lambda (k, v): (k, v))

    write_files(two_perc_genes_set, PERC_DICT_FILE, NUM_DICT_FILE,
                calc_fpath, fisher_fpath, alpha=alpha)
    print('Done. Result files:\n{}\n{}'.format(calc_fpath, fisher_fpath))
