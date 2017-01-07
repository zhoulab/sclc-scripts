#!/usr/bin/env python
import argparse
import os
import subprocess

from gene_pairs import main as generate_gene_pair_files

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--maf_file', required=True,
                    help='MAF file with columns gene/patient/effect/categ')
parser.add_argument('-l', '--sig_genes_list',
                    help='file of significant genes')
parser.add_argument('-p', '--percent_threshold', type=int,
                    help='threshold for patients count / total count')
parser.add_argument('--filter_common', action='store_true',
                    help='Filter common')
parser.add_argument('-a', '--alpha', required=True,
                    help='Fisher test significance level')
args = parser.parse_args()

if bool(args.sig_genes_list) is bool(args.percent_threshold):
    raise Exception('Required: either a gene list or percent threshold.')

# create filepaths
result_dir = os.path.join(os.path.dirname(os.getcwd()), 'results')

if args.sig_genes_list:
    suffix = os.path.basename(args.sig_genes_list)
    suffix = suffix[:suffix.rfind('.')]
else:
    suffix = 'MoreThan{}'.format(args.percent_threshold)

calc_file = 'GenePairs_{}.txt'.format(suffix)
num_file = 'GenePairsNum_{}.txt'.format(suffix)
fisher_file = 'FisherGenePairs_{}.txt'.format(suffix)
fisher_file_filtered = 'FisherGenePairs_{}_{}.txt'.format(suffix, args.alpha)

calc_file = os.path.join(result_dir, calc_file)
num_file = os.path.join(result_dir, num_file)
fisher_file = os.path.join(result_dir, fisher_file)
fisher_file_filtered = os.path.join(result_dir, fisher_file_filtered)

# generate gene pair files
args.calc_out = calc_file
args.num_out = num_file
generate_gene_pair_files(args)

# run Fisher analysis
subprocess.call(['Rscript', 'gene_pairs_fisher.R',
                 '-f', num_file,
                 '-a', args.alpha,
                 '--out_all', fisher_file,
                 '--out_filtered', fisher_file_filtered])
