import argparse
import csv
import itertools
import warnings
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from collections import OrderedDict


def create_proxy(color):
    """Generate color proxy for use in plt.legend"""
    return matplotlib.lines.Line2D([0], [0], linestyle='none',
                                   mfc=color, mec='none', marker='s')


def get_fixed_colors(n):
    """Return visually contrasting RBG hex codes, `n` < 26"""
    # these colors seem to be visually contrasting
    colors = get_spaced_colors(25)
    if not 0 <= n <= 25:
        raise ValueError('unsupported number of colors')
    return colors[:n]


def get_spaced_colors(n):
    """Return list of RGB hex codes that are evenly spaced numerically"""
    max_value = 16581375  # 255**3
    interval = int(max_value / n)
    return ['#' + hex(i)[2:].zfill(6)
            for i in range(0, max_value, interval)][:-1]


def get_mutsig_gene_pvals(mutsig_genes_file, p_value=None):
    """Return dictionary for gene p-values from `mutsig_genes_file`

    Return
    ------
    dict : { gene (str) : p-value (float) }
    """
    with open(mutsig_genes_file, 'rU') as f:
        dialect = csv.Sniffer().sniff(f.read(1024))
        f.seek(0)
        reader = csv.reader(f, delimiter='\t', dialect=dialect)
        col = reader.next().index('p')  # find column index and ignore header
        if p_value:
            return {row[0]: float(row[col])
                    for row in reader
                    if float(row[col]) <= p_value}
        return {row[0]: float(row[col])
                for row in reader}


def get_gene_info(mutation_tsv_file, mutsig_genes_file,
                  p_value=None, gene_list_file=None,
                  num_genes=None):
    """
    Return
    ------
    genes_info (dict of dicts)
    {
        gene : {
                        'p' : p-value (float)
                    'score' : score (float)
                     sample : mutation (str)
               }
    }
    """

    if gene_list_file:
        with open(gene_list_file) as f:
            genes_and_scores = OrderedDict([(line.split()[0], float(line.split()[1])
                                            if len(line.split()) == 2 else None)
                                            for line in f])
    else:
        all_genes_info = get_mutsig_gene_pvals(mutsig_genes_file, p_value)

    genes_info = OrderedDict()
    with open(mutation_tsv_file) as f:
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if not row:
                warnings.warn('Samples file: line {} is empty'.format(i + 1),
                              stacklevel=2)
            elif ((gene_list_file and (row[0] in genes_and_scores)) or
                  ((not gene_list_file) and (row[0] in all_genes_info))):
                # new gene
                if row[0] not in genes_info:
                    genes_info[row[0]] = {}
                    if gene_list_file:
                        genes_info[row[0]]['score'] = genes_and_scores[row[0]]
                    else:
                        genes_info[row[0]]['p'] = all_genes_info[row[0]]
                if row[1] in genes_info[row[0]]:
                    warnings.warn('Gene {} already has a mutation entry ({}) for sample {}.\n'
                                  'Overwriting with value {}.'.format(row[0], genes_info[row[0]][row[1]], row[1], row[3]),
                                  stacklevel=2)
                genes_info[row[0]][row[1]] = row[3]
    if gene_list_file:
        sorted_genes = sorted(genes_and_scores, key=lambda g: genes_and_scores[g],
                              reverse=True)
    else:
        sorted_genes = sorted(genes_info, key=lambda g: genes_info[g]['p'])
    genes_info = OrderedDict((gene, genes_info[gene]) for gene in sorted_genes)
    if num_genes:
        genes_info = genes_info[:num_genes]
    return genes_info


def generate_mut_plot(mutsig_genes_file, mutation_tsv_file, output,
                      sample_id_output, show_p_values, show_percent_graph,
                      mut_type_map, num_genes=None,
                      p_value=None, gene_list_file=None):
    """Generate PDF of mutation type plot.

    Parameters
    ----------
    mutsig_genes_file : filepath
        filepath for MutSig genes file
    mutation_tsv_file : filepath
        filepath for Mutation TSV file with columns gene/patient/effect/categ
    output : filepath
        plot PDF destination
    sample_id_output : filepath
        optional filepath to save sample IDs in the format (id, name)
    show_p_values : boolean
        optional generate p-value graph [-log(p)]
    show_percent_graph : boolean
        optional generate stacked bar graph for mutation distributions
    mut_type_map : OrderedDict
        mutation type names for legend
    num_genes : int
        limit on number of genes to display
    p_value : float
        p_value upper limit filter
    gene_list_file : filepath
        filepath for gene list with {gene name}\t{score}
    """
    genes_info = get_gene_info(mutation_tsv_file, mutsig_genes_file,
                               p_value, gene_list_file)

    genes_list = genes_info.keys()
    samples = list(set(itertools.chain(*[genes_info[key].keys()
                                         for key in genes_info])
                       ) - {'p', 'score'})
    for gene in list(reversed(genes_list)):
        samples = sorted(samples, key=lambda mut: mut not in genes_info[gene])
    sample_nums = [(i + 1, s) for i, s in enumerate(samples)]
    if sample_id_output:
        with open(sample_id_output, 'w') as id_file:
            fwriter = csv.writer(id_file, delimiter='\t')
            for sample_info in sample_nums:
                fwriter.writerow(sample_info)

    ncols = len(genes_list)
    nrows = len(samples)

    image = np.zeros(nrows * ncols)
    image = image.reshape((ncols, nrows))

    for r, gene in enumerate(genes_list):
        for c, sample in enumerate(samples):
            if sample in genes_info[gene]:
                image[r][c] = genes_info[gene][sample]
            else:
                image[r][c] = None

    plt.subplot(132)
    plt.rcParams["figure.figsize"] = [20.0, 5.0]

    fig, ax1 = plt.subplots()
    ax1.tick_params(length=0)

    colors = ['#ad2a1a', '#da621e', '#d3b53d', '#829356', '#0d3d56']
    color_map = ListedColormap(colors)
    proxies = [create_proxy(item) for item in colors]
    plt.legend(proxies, mut_type_map.values(), numpoints=1, loc='lower center',
               markerscale=2, fontsize=12, ncol=5,
               bbox_to_anchor=[0.5, 1.02])

    comut_grid = ax1.imshow(image, aspect='auto', interpolation='nearest')
    comut_grid.set_cmap(color_map)

    ax1.xaxis.tick_top()
    ax1.xaxis.set_ticks(np.arange(-0.5, nrows - 0.5))
    ax1.yaxis.tick_right()
    ax1.yaxis.set_ticks(np.arange(-0.5, ncols - 0.5))
    ax1.yaxis.set_ticks(np.arange(0, ncols), minor=True)
    ax1.set_xticklabels([i for i, __ in sample_nums],
                        rotation='vertical', fontsize=4, ha='left')
    ax1.set_yticklabels([])
    ax1.set_yticklabels(genes_list, fontsize=13, va='center', minor=True)
    ax1.tick_params(axis='both', which='both', length=0)
    ax1.grid(linestyle='solid', color=((.3,) * 3))
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # -log(p)
    if show_p_values:
        ax2 = fig.add_subplot(133)
        y_pos = range(len(genes_list))
        neg_log_p = [math.log(genes_info[gene]['p']) * -1
                     for gene in genes_info]
        ax2.barh(y_pos, sorted(neg_log_p), 1, align='edge',
                 color='#737373', edgecolor='black', linewidth=1, log=True)
        ax2.set_yticklabels([])
        ax2.set_xlabel('-log(p-value)')
        box = ax2.get_position()
        ax2.set_position([box.x0 + box.width * 0.71, box.y0,
                          box.width * 0.2, box.height])
        ax2.set_ylim([0, len(genes_list)])

    # stacked bars
    if show_percent_graph:
        ax3 = fig.add_subplot(131)
        mut_summary = {}
        for gene in genes_list:
            for s, m in genes_info[gene].items():
                if s != 'p':
                    if gene not in mut_summary:
                        mut_summary[gene] = {}
                        for mut in mut_type_map:
                            mut_summary[gene][mut] = 0
                    mut_summary[gene][m] += 1
        totals = [sum([mut_summary[g][m]
                       for m in mut_type_map])
                  for g in genes_list]
        bar_data = {}
        for mut in mut_type_map:
            bar_data[mut] = [(float(n) / total * 100) for n, total in
                             zip([mut_summary[g][mut] for g in genes_list], totals)]
            bar_l = list(reversed([i for i in range(len(genes_list))]))
            bottom = [0] * len(genes_list)
            for m in mut_type_map.keys()[:mut_type_map.keys().index(mut)]:
                bottom = [base + last
                          for base, last in zip(bottom, bar_data[m])]
            ax3.barh(bar_l, bar_data[mut], 1, left=bottom,
                     color=colors[mut_type_map.keys().index(mut)], linewidth=0)
        ax3.set_xlim([0, 100])
        ax3.set_ylim([0, len(genes_list)])
        ax3.tick_params(left='off', labelleft='off')
        ax3.set_xlabel('%')
        ax3.set_xticklabels([0, 20, 40, 60, 80, 100], rotation='vertical')
        ax3.get_xaxis().set_tick_params(direction='out')

    plt.savefig(output, bbox_inches='tight', format='pdf')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate comutation plots')
    # TODO: mutsig_genes_file is not relevant when gene_list_file is specified
    parser.add_argument('--mutsig_genes_file', required=True,
                        help='MutSig genes output file')
    parser.add_argument('--mutation_tsv_file', required=True,
                        help='Mutation TSV file with columns gene/patient/effect/categ')
    parser.add_argument('-o', '--output', required=True,
                        help='output PDF filepath')
    parser.add_argument('--sample_id_output', required=False,
                        help='output sample ID filepath')
    parser.add_argument('--show_p_values', action='store_true',
                        help='option to show p-values')
    parser.add_argument('--show_percent_graph', action='store_true',
                        help='option to show percentage distribution graph')
    parser.add_argument('-n', '--num_genes', type=int, required=False,
                        help='Number of genes to display')
    # specify one of the following:
    parser.add_argument('-p', '--p_value', type=float,
                        help='Sample gene mutation p-value significance cutoff')
    parser.add_argument('-g', '--gene_list_file', required=False,
                        help='Filepath for list of genes')
    args = parser.parse_args()
    mut_type_map = OrderedDict([('3', 'C:G transitions'),
                                ('4', 'C:G transversions'),
                                ('5', 'A:T transitions'),
                                ('6', 'A:T transversions'),
                                ('7', 'null+indel mutations')])
    generate_mut_plot(mut_type_map=mut_type_map, **args.__dict__)
