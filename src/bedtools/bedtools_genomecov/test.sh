#!/bin/bash

# exit on error
set -e

## VIASH START
meta_executable="target/executable/bedtools/bedtools_intersect/bedtools_intersect"
meta_resources_dir="src/bedtools/bedtools_intersect"
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

# Create and populate input files
printf "chr1\t248956422\nchr2\t198295559\nchr3\t242193529\n" > "test_data/genome.txt"
printf "chr2\t128\t228\tmy_read/1\t37\t+\nchr2\t428\t528\tmy_read/2\t37\t-\n" > "test_data/example.bed"
printf "chr2\t128\t228\tmy_read/1\t60\t+\t128\t228\t255,0,0\t1\t100\t0\nchr2\t428\t528\tmy_read/2\t60\t-\t428\t528\t255,0,0\t1\t100\t0\n" > "test_data/example.bed12"

# expected output
cat > "test_data/expected_default.bed" <<EOF
chr2	0	198295359	198295559	0.999999
chr2	1	200	198295559	1.0086e-06
chr1	0	248956422	248956422	1
chr3	0	242193529	242193529	1
genome	0	689445310	689445510	1
genome	1	200	689445510	2.90088e-07
EOF

# Test 1: 
mkdir test1
cd test1

echo "> Run bedtools_genomecov on BED file"
"$meta_executable" \
  --input "../test_data/example.bed" \
  --genome "../test_data/genome.txt" \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_default.bed"
echo "- test1 succeeded -"

cd ..

# Test 2: ibam option and pair end option and fragment size option

# Test 3: depth option

# Test 4: strand option

# Test 5: 5' end option

# Test 6: max option

# Test 7: scale option

# Test 8: trackopts option


echo "---- All tests succeeded! ----"
exit 0
