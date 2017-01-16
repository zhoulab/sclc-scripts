#!/usr/bin/env bash

# This script runs the pipeline for each of our desired gene lists.
python -u run_pipeline.py -i ../data/4dataset_nonsilent.txt --mutsig_genes_file ../data/SigGenes_001.txt --result_dir ../results/GenePairs_SigGenes_001/
python -u run_pipeline.py -i ../data/4dataset_nonsilent.txt --mutsig_genes_file ../data/SigGenes_001.txt --result_dir ../results/GenePairs_SigGenes_005/
python -u run_pipeline.py -i ../data/4dataset_nonsilent.txt -p 2 --result_dir ../results/GenePairs_MoreThan2/
