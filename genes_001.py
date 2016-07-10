import os
import csv
import itertools
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from collections import OrderedDict

import numpy as np

BASE_DIR = os.path.dirname(os.getcwd())

OUTPUT_DIRECTORY = os.path.join(BASE_DIR, 'results')
if not os.path.exists(OUTPUT_DIRECTORY):
    os.makedirs(OUTPUT_DIRECTORY)

FILE = os.path.join(BASE_DIR, 'data/genenonsilent200.txt')
GENES_FILE = os.path.join(BASE_DIR, 'data/SigGenes_001.txt')
MUT_PLOT_DEST = os.path.join(BASE_DIR, 'results/SCLC_comut_plot_001.jpg')
MUT_TYPES_PLOT_DEST = os.path.join(BASE_DIR, 'results/SCLC_comut_type_plot_001.jpg')


def get_genes_info(filepath):
    """
    Returns
    -------
    genes_info : dict
        {gene: p-val}
    """
    with open(filepath, 'rU') as genes:
        reader = csv.reader(genes, delimiter='\t', dialect=csv.excel_tab)
        reader.next()  # ignore header
        genes_info = {}
        for line in reader:
            genes_info[line[0]] = float(line[13])
        return genes_info


def create_proxy(color):
    """Generate color proxy for use in plt.legend"""
    return matplotlib.lines.Line2D([0], [0], linestyle='none',
                                   mfc=color, mec='none', marker='s')


def check_mut_type(s):
    """
    The names of some values vary between different files.
    Put duplicate definitions here.

    Variables
    ---------
    dup_map : dict
        {k: v} *k* is duplicate, *v* is display value
    """
    dup_map = {'Missense_Mutation': 'missense',
               'Nonsense_Mutation': 'nonsense',
               'Silent': 'silent',
               'Splice_Site': 'splice'}
    if s in dup_map:
        return dup_map[s]
    return s


def get_mut_vals(files, col):
    """Return distinct mutation data from all files, under index *col*"""
    if type(files) is str:
        files = [files]
    vals = set()
    for file in files:
        with open(file) as f:
            reader = csv.reader(f, delimiter='\t')
            for line in reader:
                val = check_mut_type(line[col])
                if val not in vals:
                    vals.add(val)
    return vals


def get_gene_muts():
    """
    Return
    ------
    genes_info (dict of dict)
    {
        gene : {
                      'p' : p-value (float)
                   sample : mutation (str)
               }
    }
    """
    all_genes_info = get_genes_info(GENES_FILE)
    genes_info = {}
    with open(FILE) as gene_nonsilent:
        reader = csv.reader(gene_nonsilent, delimiter='\t')
        for line in reader:
            if line[0] in all_genes_info:
                if line[0] not in genes_info:
                    genes_info[line[0]] = {}
                    genes_info[line[0]]['p'] = all_genes_info[line[0]]
                genes_info[line[0]][line[1]] = line[3]

    return genes_info


def get_gene_mut_types(files):
    """
    Return
    ------
    genes_info (dict of dict)
    {
        gene : {
                      'p' : p-value (float)
                   sample : [mut types] (list of str)
               }
    }
    """
    all_genes_info = get_genes_info(GENES_FILE)
    genes_info = {}
    for file in files:
        with open(file) as gene_nonsilent:
            reader = csv.reader(gene_nonsilent, delimiter='\t')
            for line in reader:
                if line[0] in all_genes_info:
                    if line[0] not in genes_info:
                        genes_info[line[0]] = {}
                        genes_info[line[0]]['p'] = all_genes_info[line[0]]
                    if line[1] not in genes_info[line[0]]:
                        genes_info[line[0]][line[1]] = []
                    mut_type = check_mut_type(line[2])
                    genes_info[line[0]][line[1]].append(mut_type)

    return genes_info


def generate_mut_type_plot(save_file, ignore=None):
    """
    Data Structures
    ---------------
    mut_type_counts (dict of dict)
    {
        gene : {
                   type : count
               }
    }

    """
    if not ignore:
        ignore = []
    files = [os.path.join(BASE_DIR, 'data/Peifer2columns.txt'),
             os.path.join(BASE_DIR, 'data/Rudin_59samples.txt')]
    genes_info = get_gene_mut_types(files)

    mut_type_counts = {}
    mut_vals = list(get_mut_vals(files, 2))  # ordered to generate bar
    for gene, d in genes_info.items():
        mut_type_counts[gene] = {}
        for mut_type in mut_vals:
            mut_type_counts[gene][mut_type] = 0
        for k, v in d.items():
            if k != 'p':
                for mut_type in v:
                    mut_type_counts[gene][mut_type] += 1

    genes_list = sorted(genes_info.keys(), key=lambda g: genes_info[g]['p'])

    plt.subplot()
    plt.rcParams["figure.figsize"] = [10.0, 10.0]
    fig, ax = plt.subplots()

    color = OrderedDict()
    proxy = OrderedDict()
    for m, c in zip(mut_vals, get_fixed_colors(len(mut_vals))):
        color[m] = c
        if m not in ignore:
            proxy[m] = create_proxy(c)

    plt.legend(reversed(proxy.values()), reversed(proxy.keys()),
               numpoints=1, loc='upper left', markerscale=2,
               bbox_to_anchor=(1, 1), fontsize=10)

    totals = [sum([mut_type_counts[g][mt]
                   for mt in mut_vals if mt not in ignore])
              for g in genes_list]
    bar_data = {}
    for mut in mut_vals:
        if mut in ignore:
            bar_data[mut] = [0] * len(genes_list)
        else:
            bar_data[mut] = [(float(n) / total * 100) for n, total in
                             zip([mut_type_counts[g][mut]
                                  for g in genes_list], totals)]
        bar_l = range(len(genes_list))
        bottom = [0] * len(genes_list)
        for m in mut_vals[:mut_vals.index(mut)]:
            bottom = [base + last for base, last in zip(bottom, bar_data[m])]
        ax.bar(bar_l, bar_data[mut], 1, bottom=bottom,
               color=color[mut], linewidth=0)

    ax.set_ylim([0, 100])
    ax.get_yaxis().set_tick_params(direction='out')
    ax.set_ylabel('%', rotation=0)
    ax.set_xticklabels(genes_list, ha='left', rotation=45)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0,
                     box.width * 0.8, box.height])

    plt.savefig(save_file)


def generate_mut_plot():
    genes_info = get_gene_muts()

    genes_list = sorted(genes_info.keys(), key=lambda g: genes_info[g]['p'])
    samples = list(set(itertools.chain(*[genes_info[key].keys()
                                         for key in genes_info])
                       ) - {'p'})
    for gene in list(reversed(genes_list)):
        samples = sorted(samples, key=lambda mut: mut not in genes_info[gene])

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
    mut_types = OrderedDict([('3', 'C:G transitions'), ('4', 'C:G transversions'),
                             ('5', 'A:T transitions'), ('6', 'A:T transversions'),
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
    ax1.set_xticklabels(list(samples), rotation='vertical', fontsize=4, ha='left')
    ax1.set_yticklabels(genes_list, fontsize=13, va='top')
    ax1.grid(linestyle='solid', color=((.3,) * 3))
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # -log(p)

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

    totals = [sum([mut_summary[g][m] for m in mut_types]) for g in genes_list]

    bar_data = {}
    for mut in mut_types:
        bar_data[mut] = [(float(n) / total * 100) for n, total in
                         zip([mut_summary[g][mut] for g in genes_list], totals)]
        bar_l = list(reversed([i for i in range(len(genes_list))]))
        bottom = [0] * len(genes_list)
        for m in mut_types.keys()[:mut_types.keys().index(mut)]:
            bottom = [base + last for base, last in zip(bottom, bar_data[m])]
        ax3.barh(bar_l, bar_data[mut], 1, left=bottom,
                 color=colors[mut_types.keys().index(mut)], linewidth=0)

    ax3.set_xlim([0, 100])
    ax3.set_ylim([0, len(genes_list)])
    ax3.tick_params(left='off', labelleft='off')
    ax3.set_xlabel('%')
    ax3.set_xticklabels([0, 20, 40, 60, 80, 100], rotation='vertical')
    ax3.get_xaxis().set_tick_params(direction='out')

    # adjust box positions

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

    plt.savefig(MUT_PLOT_DEST)


def get_fixed_colors(n):
    """Return visually contrasting RBG hex codes n < 26"""
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


if __name__ == "__main__":
    generate_mut_type_plot(MUT_TYPES_PLOT_DEST)
