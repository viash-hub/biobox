#!/bin/bash

set -e

## VIASH START
meta_executable="target/executable/seqtk/seqtk_sample"
meta_resources_dir="src/seqtk/seqtk_sample"
## VIASH END

#########################################################################################
mkdir seqtk_sample_se
cd seqtk_sample_se

echo "> Run seqtk_sample on fastq SE"
"$meta_executable" \
  --input "$meta_resources_dir/test_data/reads/a.fastq" \
  --seed 42 \
  --fraction_number 3 \
  --output "sampled.fastq"

echo ">> Check if output exists"
if [ ! -f "sampled.fastq" ]; then
    echo ">> sampled.fastq.gz does not exist"
    exit 1
fi

#########################################################################################
cd ..
mkdir seqtk_sample_pe_number
cd seqtk_sample_pe_number

echo ">> Run seqtk_sample on fastq.gz PE with number of reads"
"$meta_executable" \
  --input "$meta_resources_dir/test_data/reads/a.1.fastq.gz" \
  --seed 42 \
  --fraction_number 3 \
  --output "sampled_1.fastq"

"$meta_executable" \
  --input "$meta_resources_dir/test_data/reads/a.2.fastq.gz" \
  --seed 42 \
  --fraction_number 3 \
  --output "sampled_2.fastq"

echo ">> Check if output exists"
if [ ! -f "sampled_1.fastq" ] || [ ! -f "sampled_2.fastq" ]; then
    echo ">> One or both output files do not exist"
    exit 1
fi

echo ">> Compare reads"
# Extract headers
headers1=$(grep '^@' sampled_1.fastq | sed -e's/ 1$//' | sort)
headers2=$(grep '^@' sampled_2.fastq | sed -e 's/ 2$//' | sort)

# Compare headers
diff <(echo "$headers1") <(echo "$headers2") || { echo "Mismatch detected"; exit 1; }

#########################################################################################
cd ..
mkdir seqtk_sample_pe_fraction
cd seqtk_sample_pe_fraction

echo ">> Run seqtk_sample on fastq.gz PE with fraction of reads"
"$meta_executable" \
  --input "$meta_resources_dir/test_data/reads/a.1.fastq.gz" \
  --seed 42 \
  --fraction_number 0.5 \
  --output "sampled_1.fastq"

"$meta_executable" \
  --input "$meta_resources_dir/test_data/reads/a.2.fastq.gz" \
  --seed 42 \
  --fraction_number 0.5 \
  --output "sampled_2.fastq"

echo ">> Check if output exists"
if [ ! -f "sampled_1.fastq" ] || [ ! -f "sampled_2.fastq" ]; then
    echo ">> One or both output files do not exist"
    exit 1
fi

echo ">> Compare reads"
# Extract headers
headers1=$(grep '^@' sampled_1.fastq | sed -e's/ 1$//' | sort)
headers2=$(grep '^@' sampled_2.fastq | sed -e 's/ 2$//' | sort)

# Compare headers
diff <(echo "$headers1") <(echo "$headers2") || { echo "Mismatch detected"; exit 1; }

