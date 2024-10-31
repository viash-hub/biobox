#!/bin/bash

set -e

## VIASH START
## VIASH END

# Define the output directory and test data paths
output_dir="${meta_resources_dir}/output"

mkdir test_data
cp -r "/opt/cellranger-8.0.0/external/cellranger_tiny_fastq/" test_data/
cp -r "/opt/cellranger-8.0.0/external/cellranger_tiny_ref/" test_data/

input_dir="${meta_resources_dir}/test_data/cellranger_tiny_fastq"
reference_dir="${meta_resources_dir}/test_data/cellranger_tiny_ref"

# Run the tests
echo "> Running tests for ${meta_executable_name}"

# Test with folder input
"$meta_executable" \
  --input "$input_dir" \
  --reference "$reference_dir" \
  --output "$output_dir" \
  --lanes 1

# Check if output file exists
if [ ! -f "$output_dir/filtered_feature_bc_matrix.h5" ]; then
  echo "Test failed: No output was created for folder input."
  exit 1
fi

# Test with fastq files
"$meta_executable" \
  --input "$input_dir/tinygex_S1_L001_R1_001.fastq.gz" \
  --input "$input_dir/tinygex_S1_L001_R2_001.fastq.gz" \
  --reference "$reference_dir" \
  --output "$output_dir"

# Check if output file exists
if [ ! -f "$output_dir/filtered_feature_bc_matrix.h5" ]; then
  echo "Test failed: No output was created for fastq files."
  exit 1
fi

# Additional tests can be added here following the same pattern

echo "All tests succeeded!"