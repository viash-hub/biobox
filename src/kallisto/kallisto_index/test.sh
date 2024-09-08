#!/bin/bash

echo ">>>Test1: Testing $meta_functionality_name with non-default k-mer size"

"$meta_executable" \
  --input "$meta_resources_dir/test_data/transcriptome.fasta" \
  --kallisto_index Kallisto \
  --kmer_size 21


echo ">>> Checking whether output exists and is correct"
[ ! -f "Kallisto" ] && echo "Kallisto index does not exist!" && exit 1
[ ! -s "Kallisto" ] && echo "Kallisto index is empty!" && exit 1

kallisto inspect Kallisto 2> test.txt
grep "number of k-mers: 2,978" test.txt || { echo "The content of the index seems to be incorrect." && exit 1; }

################################################################################

echo ">>>Test2: Testing $meta_functionality_name with d_list argument"

"$meta_executable" \
  --input "$meta_resources_dir/test_data/transcriptome.fasta" \
  --kallisto_index Kallisto \
  --d_list "$meta_resources_dir/test_data/d_list.fasta"

echo ">>> Checking whether output exists and is correct"
[ ! -f "Kallisto" ] && echo "Kallisto index does not exist!" && exit 1
[ ! -s "Kallisto" ] && echo "Kallisto index is empty!" && exit 1

kallisto inspect Kallisto 2> test.txt
grep "number of k-mers: 3,056" test.txt || { echo "The content of the index seems to be incorrect." && exit 1; }

echo "All tests succeeded!"
exit 0
