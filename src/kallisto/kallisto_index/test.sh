#!/bin/bash

echo ">>>Test 1: Testing $meta_functionality_name with non-default k-mer size"

"$meta_executable" \
  --input "$meta_resources_dir/test_data/transcriptome.fasta" \
  --index Kallisto \
  --kmer_size 21


echo ">>> Checking whether output exists and is correct"
[ ! -f "Kallisto" ] && echo "Kallisto index does not exist!" && exit 1
[ ! -s "Kallisto" ] && echo "Kallisto index is empty!" && exit 1

kallisto inspect Kallisto 2> test.txt
grep "number of k-mers: 989" test.txt || { echo "The content of the index seems to be incorrect." && exit 1; }

################################################################################

echo ">>>Test 2: Testing $meta_functionality_name with d_list argument"

"$meta_executable" \
  --input "$meta_resources_dir/test_data/transcriptome.fasta" \
  --index Kallisto \
  --d_list "$meta_resources_dir/test_data/d_list.fasta"

echo ">>> Checking whether output exists and is correct"
[ ! -f "Kallisto" ] && echo "Kallisto index does not exist!" && exit 1
[ ! -s "Kallisto" ] && echo "Kallisto index is empty!" && exit 1

kallisto inspect Kallisto 2> test.txt
grep "number of k-mers: 959" test.txt || { echo "The content of the index seems to be incorrect." && exit 1; }

echo "All tests succeeded!"
exit 0
