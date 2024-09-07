#!/bin/bash

echo ">>>Test1: Testing $meta_functionality_name"

"$meta_executable" \
  --input "$meta_resources_dir/test_data/transcriptome.fasta" \
  --kallisto_index Kallisto \
  --kmer_size 21 \
  --make_unique 

echo ">>> Checking whether output exists"
[ ! -f "Kallisto" ] && echo "Kallisto index does not exist!" && exit 1
[ ! -s "Kallisto" ] && echo "Kallisto index is empty!" && exit 1

echo ">>>Test2: Testing $meta_functionality_name"

"$meta_executable" \
  --input "$meta_resources_dir/test_data/transcriptome.fasta" \
  --kallisto_index Kallisto \
  --d_list "$meta_resources_dir/test_data/d_list.fasta"

echo "All tests succeeded!"
exit 0
