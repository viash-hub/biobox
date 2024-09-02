#!/bin/bash

echo ">>> Testing $meta_executable"

echo ">>> Testing for paired-end reads"
"$meta_executable" \
    --input_1 $meta_resources_dir/test_data/a.3.fastq.gz \
    --input_2 $meta_resources_dir/test_data/a.4.fastq.gz \
    --record_count 3 \
    --seed 1 \
    --output_1  a.1.subsampled.fastq \
    --output_2  a.2.subsampled.fastq 

echo ">> Checking if the correct files are present"
[ ! -f "a.1.subsampled.fastq" ] && echo "Subsampled FASTQ file for read 1 is missing!" && exit 1
[ $(wc -l < a.1.subsampled.fastq) -ne 12 ] && echo "Subsampled FASTQ file for read 1 does not contain the expected number of records" && exit 1
[ ! -f "a.2.subsampled.fastq" ] && echo "Subsampled FASTQ file for read 2 is missing" && exit 1
[ $(wc -l < a.2.subsampled.fastq) -ne 12 ] && echo "Subsampled FASTQ file for read 2 does not contain the expected number of records" && exit 1

rm a.1.subsampled.fastq a.2.subsampled.fastq

echo ">>> Testing for single-end reads"
"$meta_executable" \
    --input_1 $meta_resources_dir/test_data/a.3.fastq.gz \
    --record_count 3 \
    --seed 1 \
    --output_1  a.1.subsampled.fastq 

    
echo ">> Checking if the correct files are present"
[ ! -f "a.1.subsampled.fastq" ] && echo "Subsampled FASTQ file is missing" && exit 1
[ $(wc -l < a.1.subsampled.fastq) -ne 12 ] && echo "Subsampled FASTQ file does not contain the expected number of records" && exit 1

echo ">>> Tests finished successfully"
exit 0

