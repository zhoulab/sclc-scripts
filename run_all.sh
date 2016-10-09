#!/bin/bash

while getopts "l:a:f" opt; do
  case $opt in
    l) gene_list="$OPTARG"
    ;;
    a) alpha="$OPTARG"
    ;;
    f) filter=true
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

calc_file="GenePairs_${gene_list}.txt"
num_file="GenePairsNum_${gene_list}.txt"
fisher_file="FisherGenePairs_${gene_list}.txt"

py_cmd="python -u gene_pairs.py -l $gene_list --calc_out $calc_file --num_out $num_file"
result_dir="results/$gene_list"
if [ "$filter" = true ]; then
    py_cmd="$py_cmd --filter_common"
    result_dir="${result_dir}_remove_common"
fi
mkdir -p "$result_dir"
echo "Results directory: $result_dir"

module load python
. ve/bin/activate
eval $py_cmd

module load R
Rscript gene_pairs_fisher.R -f "$num_file" -a "$alpha" -o "$fisher_file"

declare -a files=("$calc_file" "$num_file" "$fisher_file")
for fname in "${files[@]}"
do
    mv "$fname" -t "$result_dir"
done
