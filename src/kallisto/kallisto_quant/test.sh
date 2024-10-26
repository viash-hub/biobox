#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Test 1: Testing for paired-end reads"
"$meta_executable" \
  --index "$meta_resources_dir/test_data/index/transcriptome.idx" \
  --rf_stranded \
  --output_dir . \
  --input "$meta_resources_dir/test_data/reads/A_R1.fastq;$meta_resources_dir/test_data/reads/A_R2.fastq"

echo ">>> Checking whether output exists"
[ ! -f "run_info.json" ] && echo "run_info.json does not exist!" && exit 1
[ ! -s "run_info.json" ] && echo "run_info.json is empty!" && exit 1
[ ! -f "abundance.tsv" ] && echo "abundance.tsv does not exist!" && exit 1
[ ! -s "abundance.tsv" ] && echo "abundance.tsv is empty!" && exit 1
[ ! -f "abundance.h5" ] && echo "abundance.h5 does not exist!" && exit 1
[ ! -s "abundance.h5" ] && echo "abundance.h5 is empty!" && exit 1

echo ">>> Checking if output is correct"
diff "abundance.tsv" "$meta_resources_dir/test_data/abundance_1.tsv" || { echo "abundance.tsv is not correct"; exit 1; }

rm -rf abundance.tsv abundance.h5 run_info.json

################################################################################

echo ">>> Test 2: Testing for single-end reads"
"$meta_executable" \
  --index "$meta_resources_dir/test_data/index/transcriptome.idx" \
  --rf_stranded \
  --output_dir . \
  --single \
  --input "$meta_resources_dir/test_data/reads/A_R1.fastq" \
  --fragment_length 101 \
  --sd 50

echo ">>> Checking whether output exists"
[ ! -f "run_info.json" ] && echo "run_info.json does not exist!" && exit 1
[ ! -s "run_info.json" ] && echo "run_info.json is empty!" && exit 1
[ ! -f "abundance.tsv" ] && echo "abundance.tsv does not exist!" && exit 1
[ ! -s "abundance.tsv" ] && echo "abundance.tsv is empty!" && exit 1
[ ! -f "abundance.h5" ] && echo "abundance.h5 does not exist!" && exit 1
[ ! -s "abundance.h5" ] && echo "abundance.h5 is empty!" && exit 1

echo ">>> Checking if output is correct"
diff "abundance.tsv" "$meta_resources_dir/test_data/abundance_2.tsv" || { echo "abundance.tsv is not correct"; exit 1; }

rm -rf abundance.tsv abundance.h5 run_info.json

################################################################################

echo "All tests succeeded!"
exit 0
