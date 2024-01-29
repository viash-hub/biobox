#!/bin/bash

set -e

dir_in="${meta_resources_dir%/}/test_data"

echo "> Run PEAR"
"$meta_executable" \
  --forward_fastq "$dir_in/reads.left.fq.gz" \
  --reverse_fastq "$dir_in/reads.right.fq.gz" \
  --output "test" \
  --p_value 0.01 \
  --threads 4 \
  --memory "4G"

echo ">> Checking output"
echo $(ls)
[ ! -f "test.assembled.fastq.gz" ] && echo "Output file test.assembled.fastq.gz does not exist" && exit 1
[ ! -f "test.unassembled.forward.fastq.gz" ] && echo "Output file test.unassembled.forward.fastq.gz does not exist" && exit 1
[ ! -f "test.unassembled.reverse.fastq.gz" ] && echo "Output file test.unassembled.reverse.fastq.gz does not exist" && exit 1
[ ! -f "test.discarded.fastq.gz" ] && echo "Output file ftest.discarded.fastq.gz does not exist" && exit 1

echo "> Test successful"