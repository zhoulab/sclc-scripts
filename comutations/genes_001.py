import argparse
import csv
import itertools
import math
import matplotlib
matplotlib.use('PDF')
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


def get_sig_gene_pvals(sig_genes_file):
    """Return dictionary for gene p-values from `sig_genes_file`

    Return
    ------
    dict : { gene (str) : p-value (float) }
    """
    with open(sig_genes_file, 'rU') as f:
        dialect = csv.Sniffer().sniff(f.read(1024))
        f.seek(0)
        reader = csv.reader(f, delimiter='\t', dialect=dialect)
        col = reader.next().index('p')  # find column index and ignore header
        return {row[0]: float(row[col])
                for row in reader}


def get_gene_muts(gene_mut_file, sig_genes_file):
    """
    Return
    ------
    genes_info (dict of dicts)
    {
        gene : {
                      'p' : p-value (float)
                   sample : mutation (str)
               }
    }
    """
    all_genes_info = get_sig_gene_pvals(sig_genes_file)
    genes_info = {}
    with open(gene_mut_file) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0] in all_genes_info:
                if row[0] not in genes_info:
                    genes_info[row[0]] = {}
                    genes_info[row[0]]['p'] = all_genes_info[row[0]]
                genes_info[row[0]][row[1]] = row[3]
    return genes_info


def generate_mut_plot(sig_genes_file, gene_mut_file, output,
                      sample_id_output, show_p_values, show_percent_graph):
    """Generate PDF of mutation type plot.

    Parameters
    ----------
    sig_genes_file : filepath
        filepath for TSV file of significant genes
    gene_mut_file : filepath
        filepath for TSV file of gene mutation information
    output : filepath
        plot PDF destination
    sample_id_output : filepath
        optional filepath to save sample IDs in the format (id, name)
    show_p_values : boolean
        optional generate p-value graph [-log(p)]
    show_percent_graph : boolean
        optional generate stacked bar graph for mutation distributions
    """
    genes_info = get_gene_muts(gene_mut_file, sig_genes_file)

    genes_list = sorted(genes_info.keys(), key=lambda g: genes_info[g]['p'])
    samples = list(set(itertools.chain(*[genes_info[key].keys()
                                         for key in genes_info])
                       ) - {'p'})
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
    mut_types = OrderedDict([('3', 'C:G transitions'),
                             ('4', 'C:G transversions'),
                             ('5', 'A:T transitions'),
                             ('6', 'A:T transversions'),
                             ('7', 'null+indel mutations')])
    color_map = ListedColormap(colors)
    proxies = [create_proxy(item) for item in colors]
    plt.legend(proxies, mut_types.values(), numpoints=1, loc='lower right',
               markerscale=2, bbox_to_anchor=(0, 1), fontsize=12)

    comut_grid = ax1.imshow(image, aspect='auto', interpolation='nearest')
    comut_grid.set_cmap(color_map)

    ax1.xaxis.tick_top()
    ax1.xaxis.set_ticks(np.arange(-0.5, nrows - 0.5))
    ax1.yaxis.tick_right()
    ax1.yaxis.set_ticks(np.arange(-0.5, ncols - 0.5))
    ax1.set_xticklabels([i for i, __ in sample_nums],
                        rotation='vertical', fontsize=4, ha='left')
    ax1.set_yticklabels(genes_list, fontsize=13, va='top')
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
                        for mut in mut_types:
                            mut_summary[gene][mut] = 0
                    mut_summary[gene][m] += 1
        totals = [sum([mut_summary[g][m]
                       for m in mut_types])
                  for g in genes_list]
        bar_data = {}
        for mut in mut_types:
            bar_data[mut] = [(float(n) / total * 100) for n, total in
                             zip([mut_summary[g][mut] for g in genes_list], totals)]
            bar_l = list(reversed([i for i in range(len(genes_list))]))
            bottom = [0] * len(genes_list)
            for m in mut_types.keys()[:mut_types.keys().index(mut)]:
                bottom = [base + last
                          for base, last in zip(bottom, bar_data[m])]
            ax3.barh(bar_l, bar_data[mut], 1, left=bottom,
                     color=colors[mut_types.keys().index(mut)], linewidth=0)
        ax3.set_xlim([0, 100])
        ax3.set_ylim([0, len(genes_list)])
        ax3.tick_params(left='off', labelleft='off')
        ax3.set_xlabel('%')
        ax3.set_xticklabels([0, 20, 40, 60, 80, 100], rotation='vertical')
        ax3.get_xaxis().set_tick_params(direction='out')

    # adjust box positions
    if show_p_values and show_percent_graph:
        box3 = ax3.get_position()
        ax3.set_position([box3.x0, box3.y0 + box3.height * 0.5,
                          box3.width * 0.2, box3.height * 0.5])
        box1 = ax1.get_position()
        ax1.set_position([box1.x0 + box3.width * 0.2, box1.y0 + box1.height * 0.5,
                          box1.width, box1.height * 0.5])
        box2 = ax2.get_position()
        ax2.set_position([box2.x0 * 1.01, box2.y0 + box2.height * 0.5,
                          box2.width, box2.height * 0.5])

        # lower everything
        for ax in (ax1, ax2, ax3):
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 - box.height * 0.5,
                             box.width, box.height])
    else:
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0 + box.height * 0.25,
                          box.width, box.height * 0.5])

    plt.savefig(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate comutation plots')
    parser.add_argument('--sig_genes_file', required=True,
                        help='TSV file of significant genes')
    parser.add_argument('--gene_mut_file', required=True,
                        help='TSV file of gene mutation information')
    parser.add_argument('-o', '--output', required=True,
                        help='output PDF filepath')
    parser.add_argument('--sample_id_output', required=False,
                        help='output sample ID filepath')
    parser.add_argument('--show_p_values', action='store_true',
                        help='option to show p-values')
    parser.add_argument('--show_percent_graph', action='store_true',
                        help='option to show percentage distribution graph')
    args = parser.parse_args()
    generate_mut_plot(**args.__dict__)
