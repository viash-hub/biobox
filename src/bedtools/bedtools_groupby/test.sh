#!/bin/bash

# exit on error
set -e

## VIASH START
meta_executable="target/executable/bedtools/bedtools_groupby/bedtools_groupby"
meta_resources_dir="src/bedtools/bedtools_groupby"
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

# Create and populate featuresA.bed

# Create and populate expected output files for different tests

# Test 1: Default intersect
mkdir test1
cd test1

# echo "> Run bedtools_intersect on BED files with default intersect"
# "$meta_executable" \
#   --input_a "../test_data/featuresA.bed" \
#   --input_b "../test_data/featuresB.bed" \
#   --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_default.bed"
echo "- test1 succeeded -"

cd ..



echo "---- All tests succeeded! ----"
exit 0
