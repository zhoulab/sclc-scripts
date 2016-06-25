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
PLOT_DEST = os.path.join(BASE_DIR, 'results/SCLC_comut_plot_001.jpg')


def get_genes_info(filepath):
    with open(filepath, 'rU') as genes:
        reader = csv.reader(genes, delimiter='\t', dialect=csv.excel_tab)
        reader.next()
        genes_info = {}
        for line in reader:
            genes_info[line[0]] = float(line[13])
        return genes_info


def create_proxy(color):
    return matplotlib.lines.Line2D([0], [0], linestyle='none',
                                   mfc=color, mec='none', marker='s')


if __name__ == "__main__":
    all_genes_info = get_genes_info(GENES_FILE)
    genes_info = {}

    with open(FILE) as gene_nonsilent:
        reader = csv.reader(gene_nonsilent, delimiter='\t')
        for line in reader:
            if line[0] in all_genes_info.keys():
                if line[0] not in genes_info.keys():
                    genes_info[line[0]] = {}
                    genes_info[line[0]]['p'] = all_genes_info[line[0]]
                genes_info[line[0]][line[1]] = int(line[3])

    genes_list = sorted(genes_info.keys(), key=lambda x: all_genes_info[x])
    samples = sorted(list(set(itertools.chain(*[genes_info[key].keys()
                                                for key in genes_info.keys()])) - {'p'}))
    ncols = len(genes_list)
    nrows = len(samples)

    image = np.zeros(nrows * ncols)
    image = image.reshape((ncols, nrows))

    for r, gene in enumerate(genes_list):
        for c, sample in enumerate(samples):
            if sample in genes_info[gene].keys():
                image[r][c] = genes_info[gene][sample]
            else:
                image[r][c] = None

    plt.subplot(221)
    plt.rcParams["figure.figsize"] = [20.0, 15.0]

    fig, ax1 = plt.subplots()
    ax1.tick_params(length=0)

    colors = ['#ad2a1a', '#da621e', '#d3b53d', '#829356', '#0d3d56']
    mut_types = OrderedDict({3: 'C:G transitions', 4: 'C:G transversions',
                             5: 'A:T transitions', 6: 'A:T transversions',
                             7: 'null+indel mutations'})
    color_map = ListedColormap(colors)
    proxies = [create_proxy(item) for item in colors]
    plt.legend(proxies, mut_types.values(), numpoints=1,
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

    ax2 = fig.add_subplot(2, 2, (2, 4))

    y_pos = range(len(genes_list))
    neg_log_p = [math.log(genes_info[gene]['p']) * -1
                 for gene in genes_info.keys()]

    ax2.barh(y_pos, sorted(neg_log_p), 1, align='edge',
             color='#737373', edgecolor='black', linewidth=1, log=True)
    ax2.set_yticklabels([])
    ax2.set_xlabel('-log(p-value)')

    box = ax2.get_position()
    ax2.set_position([box.x0 + box.width * 0.71, box.y0,
                      box.width * 0.2, box.height])
    ax2.set_ylim([0, len(genes_list)])

    # stacked bars

    ax3 = fig.add_subplot(2, 2, 3)

    sample_info = {}
    for s in samples:
        for gene in genes_list:
            if s in genes_info[gene].keys():
                if s not in sample_info.keys():
                    sample_info[s] = {}
                    for mut in range(3, 8):
                        sample_info[s][mut] = 0
                mut = genes_info[gene][s]
                sample_info[s][mut] += 1

    totals = [sum([sample_info[s][m] for m in range(3, 8)]) for s in samples]

    bar_data = {}
    for mut in range(3, 8):
        bar_data[mut] = [(float(n) / total * 100) for n, total in
                         zip([sample_info[s][mut] for s in samples], totals)]
        bar_l = [i for i in range(len(samples))]
        bottom = [0] * len(samples)
        for x in range(3, mut):
            bottom = [base + last for base, last in zip(bottom, bar_data[x])]
        ax3.bar(bar_l, bar_data[mut], bottom=bottom,
                color=colors[mut - 3], linewidth=0, width=1)

    ax3.set_xlim([0, len(samples)])
    ax3.set_ylim([0, 100])
    ax3.tick_params(bottom='off', labelbottom='off')
    ax3.set_ylabel('%', rotation=0)
    ax3.get_yaxis().set_tick_params(direction='out')

    # adjust box heights

    box1 = ax1.get_position()
    ax1.set_position([box1.x0, box1.y0 + box1.height * 0.2,
                      box1.width, box1.height * 0.2])
    box2 = ax2.get_position()
    ax2.set_position([box2.x0, box2.y0 + box2.height * 0.2,
                      box2.width, box2.height * 0.2])
    box3 = ax3.get_position()
    ax3.set_position([box3.x0, ax1.get_position().y0 - box3.height * 0.3,
                      box1.width, box3.height * 0.3])

    plt.savefig(PLOT_DEST)
