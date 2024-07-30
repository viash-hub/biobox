#!/bin/bash

# viash ns build -q bedtools_bamtofastq --setup cb

# exit on error
set -e

## VIASH START
meta_executable="target/executable/bedtools/bedtools_sort/bedtools_bamtofastq"
meta_resources_dir="src/bedtools/bedtools_bamtofastq"
## VIASH END

#############################################
# helper functions
assert_file_exists() {
  [ -f "$1" ] || { echo "File '$1' does not exist" && exit 1; }
}
assert_file_not_empty() {
  [ -s "$1" ] || { echo "File '$1' is empty but shouldn't be" && exit 1; }
}
assert_file_contains() {
  grep -q "$2" "$1" || { echo "File '$1' does not contain '$2'" && exit 1; }
}
assert_identical_content() {
  diff -a "$2" "$1" \
    || (echo "Files are not identical!" && exit 1)
}
#############################################

# Create directories for tests
echo "Creating Test Data..."
mkdir -p test_data

# Create and populate example files

# Create expected output files


# Test 1: 
mkdir test1
cd test1

echo "> Run bedtools bamtofastq on BAM file"
"$meta_executable" \
  --input "example.bam" \
  --output_fq "output.fastq"

# checks
assert_file_exists "output.fastq"
assert_file_not_empty "output.fastq"
assert_identical_content "output.fastq" "../test_data/expected.fastq"
echo "- test1 succeeded -"

cd ..

# Test 2:
