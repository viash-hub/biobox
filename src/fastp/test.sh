#!/bin/bash

set -e

## VIASH START
meta_executable="target/docker/fastp/fastp"
meta_resources_dir="src/fastp"
## VIASH END

#########################################################################################
mkdir fastp_se
cd fastp_se

echo "> Run fastp on SE"
"$meta_executable" \
  --in1 "$meta_resources_dir/test_data/se/a.fastq" \
  --out1 "trimmed.fastq" \
  --failed_out "failed.fastq" \
  --json "report.json" \
  --html "report.html" \
  --adapter_sequence ACGGCTAGCTA

echo ">> Check if output exists"
[ ! -f "trimmed.fastq" ] && echo ">> trimmed.fastq does not exist" && exit 1
[ ! -f "failed.fastq" ] && echo ">> failed.fastq does not exist" && exit 1
[ ! -f "report.json" ] && echo ">> report.json does not exist" && exit 1
[ ! -f "report.html" ] && echo ">> report.html does not exist" && exit 1

#########################################################################################
cd ..
mkdir fastp_pe_minimal
cd fastp_pe_minimal

echo ">> Run fastp on PE with minimal parameters"
"$meta_executable" \
  --in1 "$meta_resources_dir/test_data/pe/a.1.fastq" \
  --in2 "$meta_resources_dir/test_data/pe/a.2.fastq" \
  --out1 "trimmed_1.fastq" \
  --out2 "trimmed_2.fastq"

echo ">> Check if output exists"
[ ! -f "trimmed_1.fastq" ] && echo ">> trimmed_1.fastq does not exist" && exit 1
[ ! -f "trimmed_2.fastq" ] && echo ">> trimmed_2.fastq does not exist" && exit 1

#########################################################################################
cd ..
mkdir fastp_pe_many
cd fastp_pe_many

echo ">> Run fastp on PE with many parameters"
"$meta_executable" \
  --in1 "$meta_resources_dir/test_data/pe/a.1.fastq" \
  --in2 "$meta_resources_dir/test_data/pe/a.2.fastq" \
  --out1 "trimmed_1.fastq" \
  --out2 "trimmed_2.fastq" \
  --failed_out "failed.fastq" \
  --json "report.json" \
  --html "report.html" \
  --adapter_sequence ACGGCTAGCTA \
  --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --merge \
  --merged_out "merged.fastq"

echo ">> Check if output exists"
[ ! -f "trimmed_1.fastq" ] && echo ">> trimmed_1.fastq does not exist" && exit 1
[ ! -f "trimmed_2.fastq" ] && echo ">> trimmed_2.fastq does not exist" && exit 1
[ ! -f "failed.fastq" ] && echo ">> failed.fastq does not exist" && exit 1
[ ! -f "report.json" ] && echo ">> report.json does not exist" && exit 1
[ ! -f "report.html" ] && echo ">> report.html does not exist" && exit 1
[ ! -f "merged.fastq" ] && echo ">> merged.fastq does not exist" && exit 1

#########################################################################################

echo "> Test successful"