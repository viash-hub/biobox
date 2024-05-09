#!/bin/bash

test_dir="${meta_resources_dir}/test_data"
out_dir="${meta_resources_dir}/tmp"

############################################################################################

## example 1: samtools fastq -0 /dev/null in_name.bam > all_reads.fq
## example 2: samtools fastq -0 /dev/null -s single.fq -N in_name.bam > paired.fq
## example 3: samtools fastq with fasta output??
## example 4: samtools fastq with compressed input?
## example 5: samtools fastq with no suffix?


echo ">>> Test 1: Sorting a BAM file"

"$meta_executable" \
  --input "$test_dir/a.bam" \
  --output "$test_dir/a.sorted.bam"

echo ">>> Check if output file exists"
[ ] \
    && echo "Output file a.sorted.bam does not exist" && exit 1

echo ">>> Check if output is empty"

echo ">>> Check if output matches expected output"


############################################################################################

############################################################################################


echo "All tests succeeded!"
exit 0