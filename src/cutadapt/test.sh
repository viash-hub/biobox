#!/bin/bash

set -e

dir_in="$meta_resources_dir/test_data"

echo "> Run cutadapt on single-end data"
"$meta_executable" \
  --report minimal \
  --output output-dir \
  --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --input $dir_in/se/a.fastq \
  --quality_cutoff 20 \
  --json

echo ">> Checking output"
[ ! -f "output-dir/report.txt" ] && echo "report.txt does not exist" && exit 1
[ ! -f "output-dir/report.json" ] && echo "report.json does not exist" && exit 1
[ ! -f "output-dir/1_R1_001.fastq" ] && echo "1_R1_001.fastq does not exist" && exit 1
[ ! -f "output-dir/unknown_R1_001.fastq" ] && echo "1_R1_001.fastq does not exist" && exit 1

echo ">> Check if output is empty"
[ -s "output-dir/1_R1_001.fastq" ] && echo "1_R1_001.fastq should be empty" && exit 1
[ ! -s "output-dir/unknown_R1_001.fastq" ] && echo "unkown_R1_001.fastq is empty" && exit 1

rm -r output-dir

echo "> Run cutadapt on paired-end data"
"$meta_executable" \
  --report minimal \
  --output output-dir \
  --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --input $dir_in/pe/a.1.fastq \
  --input_r2 $dir_in/pe/a.2.fastq \
  --quality_cutoff 20 \
  --json \
  ---cpus 1

echo ">> Checking output"
[ ! -f "output-dir/report.txt" ] && echo "report.txt does not exist" && exit 1
[ ! -f "output-dir/report.json" ] && echo "report.json does not exist" && exit 1
[ ! -f "output-dir/1_R1_001.fastq" ] && echo "1_R1_001.fastq does not exist" && exit 1
[ ! -f "output-dir/1_R2_001.fastq" ] && echo "1_R1_001.fastq does not exist" && exit 1
[ ! -f "output-dir/unknown_R1_001.fastq" ] && echo "1_R1_001.fastq does not exist" && exit 1
[ ! -f "output-dir/unknown_R2_001.fastq" ] && echo "1_R1_001.fastq does not exist" && exit 1

echo ">> Check if output is empty"
[ -s "output-dir/1_R1_001.fastq" ] && echo "1_R1_001.fastq should be empty" && exit 1
[ -s "output-dir/1_R2_001.fastq" ] && echo "1_R2_001.fastq should be empty" && exit 1
[ ! -s "output-dir/unknown_R1_001.fastq" ] && echo "unkown_R1_001.fastq is empty" && exit 1
[ ! -s "output-dir/unknown_R2_001.fastq" ] && echo "unkown_R2_001.fastq is empty" && exit 1

rm -r output-dir

echo "> Test successful"

