#!/bin/bash

set -e

dir_in="${meta_resources_dir%/}/test_data"

#############################################
mkdir uniform
cd uniform

echo "> Run lofreq indelqual uniform"
"$meta_executable" \
  --input "$dir_in/test.bam" \
  -u 15 \
  --out "uniform.bam" \

echo ">> Checking output"
[ ! -f "uniform.bam" ] && echo "Output file uniform.bam does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "uniform.bam" ] && echo "Output file uniform.bam is empty" && exit 1

cd ..

#############################################
mkdir dindel
cd dindel

echo "> run lofreq indelqual dindel"
"$meta_executable" \
  --input "$dir_in/test.bam" \
  --ref "$dir_in/test.fa" \
  --dindel \
  --out "dindel.bam"

echo ">> Checking output"
[ ! -f "dindel.bam" ] && echo "Output file dindel.bam does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "dindel.bam" ] && echo "Output file dindel.bam is empty" && exit 1

cd ..

#############################################

echo "> Test successful"
