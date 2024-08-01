#!/bin/bash

# exit on error
set -e

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

# Test 1: Default sort on BED file
mkdir test1
cd test1

# echo "> Run bedtools_sort on BED file"
# "$meta_executable" \
#   --input "../test_data/featureA.bed" \
#   --output "output.bed"

# # checks
# assert_file_exists "output.bed"
# assert_file_not_empty "output.bed"
# assert_identical_content "output.bed" "../test_data/expected_sorted_A.bed"
# echo "- test1 succeeded -"

cd ..

echo "---- All tests succeeded! ----"
exit 0
