#!/bin/bash

set -e

dir_in="${meta_resources_dir%/}/test_data"

echo "> Run lofreq indelqual uniform"
"$meta_executable" \
  --input "$dir_in/a.bam" \
  -u 15 \
  --out "uniform.bam" \

echo ">> Checking output"
[ ! -f "out.bam" ] && echo "Output file out.bam does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "out.bam" ] && echo "Output file out.bam is empty" && exit 1

echo "> run lofreq indelqual dindel"
"$meta_executable" \
  --input "$dir_in/a.bam" \
  --ref "$dir_in/hg38_chr21.fa" \
  --dindel \
  --out "dindel.bam"

echo ">> Checking output"
[ ! -f "dindel.bam" ] && echo "Output file dindel.bam does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "dindel.bam" ] && echo "Output file dindel.bam is empty" && exit 1


echo "> Test successful"
