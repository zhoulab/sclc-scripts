#!/bin/sh
cd /ufrc/zhou/share/projects/bioinformatics/gene-pairs
module load python
. ve/bin/activate
python -u gene_pairs.py "$@"
