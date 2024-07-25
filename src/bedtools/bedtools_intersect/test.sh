#!/bin/bash

# exit on error
set -e

## VIASH START
meta_executable="target/executable/bedtools/bedtools_intersect"
meta_resources_dir="src/bedtools/bedtools_intersect"
## VIASH END

# Create directories for tests
echo "Creating Test Data..."
mkdir test_data

# Create and populate input_a.bed

# Create and populate input_b.bed

# Create and populate expected.bed

# Run basic test
mkdir test1
cd test1

echo "> Run bedtools_intersect on BED files"
"$meta_executable" \
  --input_a "../test_data/input_a.bed" \
  --input_b "../test_data/input_b.bed" \
  --output "output.bed"



