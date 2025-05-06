#!/bin/bash

set -e

## VIASH START
## VIASH END

# create temporary directory
tmp_dir=$(mktemp -d "${meta_temp_dir}/${meta_name}-XXXXXXXX")
function clean_up {
    rm -rf "${tmp_dir}"
}
trap clean_up EXIT

echo "Copy test data from Cell Ranger installation directory"
mkdir -p "${tmp_dir}/test_data/"
cp -r "/opt/cellranger-8.0.0/external/cellranger_tiny_fastq/" "${tmp_dir}/test_data"
cp -r "/opt/cellranger-8.0.0/external/cellranger_tiny_ref/" "${tmp_dir}/test_data"
input_dir="${tmp_dir}/test_data/cellranger_tiny_fastq"
reference_dir="${tmp_dir}/test_data/cellranger_tiny_ref"


## TEST 1: run with folder input
echo "Running ${meta_name} with folder input"
output_dir="${tmp_dir}/test1/"
mkdir -p "${output_dir}"

"${meta_executable}" \
  --fastqs "${input_dir}" \
  --transcriptome "${reference_dir}" \
  --output "${output_dir}" \
  --lanes 1

[[ $? != 0 ]] && echo "Non zero exit code: $?" && exit 1
[[ ! -f "${output_dir}/filtered_feature_bc_matrix.h5" ]] && echo "Output file could not be found!" && exit 1
[[ -f "${output_dir}/possorted_genome_bam.bam" ]] && echo "Output file should not be found!" && exit 1
[[ ! -d "${output_dir}/analysis" ]] && echo "Analysis output directory should exist!" && exit 1

## TEST 2: run with individual file input
echo "Running ${meta_name} with individual file input"
output_dir="${tmp_dir}/test2/"
mkdir -p "${output_dir}"
"${meta_executable}" \
  --fastqs "${input_dir}/tinygex_S1_L001_R1_001.fastq.gz" \
  --fastqs "${input_dir}/tinygex_S1_L001_R2_001.fastq.gz" \
  --transcriptome "${reference_dir}" \
  --output "${output_dir}" \
  --no_secondary \
  --create_bam

[[ $? != 0 ]] && echo "Non zero exit code: $?" && exit 1
[[ ! -f "${output_dir}/filtered_feature_bc_matrix.h5" ]] && echo "Output file could not be found!" && exit 1
[[ ! -f "${output_dir}/possorted_genome_bam.bam" ]] && echo "Output file could not be found!" && exit 1
[[ -d "${output_dir}/analysis" ]] && echo "Analysis output directory should not exist!" && exit 1

echo "All tests succeeded!"