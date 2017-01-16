#!/usr/bin/env python
import argparse
import os
import subprocess

from gene_pairs import get_unique_samples, main as generate_gene_pair_files

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--maf_file', required=True,
                    help='MAF file with columns gene/patient/effect/categ')
parser.add_argument('--mutsig_genes_file',
                    help='MutSig genes output file')
parser.add_argument('-p', '--percent_threshold', type=int,
                    help='threshold for patients count / total count')
parser.add_argument('--filter_common', action='store_true',
                    help='Filter common')
parser.add_argument('--fisher_alpha', required=True,
                    help='Fisher test significance level')
parser.add_argument('--result_dir', required=True,
                    help='Result directory')
args = parser.parse_args()

if bool(args.mutsig_genes_file) is bool(args.percent_threshold):
    raise Exception('Required: either a MutSig output file or percent threshold.')

# create filepaths
if not os.path.exists(args.result_dir):
    os.makedirs(args.result_dir)
if args.mutsig_genes_file:
    suffix = os.path.basename(args.mutsig_genes_file)
    suffix = suffix[:suffix.rfind('.')]
else:
    suffix = 'MoreThan{}'.format(args.percent_threshold)

filepaths = {
    'calc_file': 'GenePairs_{}.txt'.format(suffix),
    'num_file': 'GenePairsNum_{}.txt'.format(suffix),
    'patient_count_file': 'PatientCounts.txt',
    'patient_count_dup_file': 'PatientCountsDuplicateGenes.txt',
    'fisher_file': 'FisherGenePairs_{}.txt'.format(suffix),
    'fisher_file_filtered' : 'FisherGenePairs_{}_{}.txt'.format(suffix, args.fisher_alpha),
    'co_occurrence_plot': 'FrequencyPlot_{}.pdf'.format(suffix),
    'co_occurrence_table': 'FrequencyTable_{}.txt'.format(suffix)
}

for f, filepath in filepaths.items():
    filepaths[f] = os.path.join(args.result_dir, filepath)

# generate gene pair files
args.calc_out = filepaths['calc_file']
args.num_out = filepaths['num_file']
args.patient_count_out = filepaths['patient_count_file']
args.patient_count_dup_out = filepaths['patient_count_dup_file']
generate_gene_pair_files(args)

# run Fisher analysis
subprocess.call(['Rscript', 'gene_pairs_fisher.R',
                 '-i', filepaths['num_file'],
                 '-a', args.fisher_alpha,
                 '-n', str(len(get_unique_samples(args.maf_file))),
                 '--out_all', filepaths['fisher_file'],
                 '--out_filtered', filepaths['fisher_file_filtered']])
# genereate co-occurrence frequency plot and table
subprocess.call(['Rscript', 'co-occurrence_freq_plot.R',
                 '-i', filepaths['calc_file'],
                 '-o', filepaths['co_occurrence_plot'],
                 '--freq_table_out', filepaths['co_occurrence_table']])
