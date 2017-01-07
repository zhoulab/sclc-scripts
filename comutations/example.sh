#!/usr/bin/env bash

python genes_001.py --show_p_values --show_percent_graph --mutsig_genes_file ../data/SigGenes_001.txt --gene_mut_file ../data/genenonsilent200.txt -o ../results/SCLC_comut_plot_001.pdf --sample_id_out ../results/SampleIDs.txt
