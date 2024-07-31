#!/bin/bash

# exit on error
set -e

# directory of the bam file
test_data="$meta_resources_dir/test_data"

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

mkdir -p data

# Generate expected files for comparison
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t60\t+\nchr2:172936693-172938111\t428\t528\tmy_read/2\t60\t-\n" > data/expected.bed
printf "chr2:172936693-172938111\t128\t228\tchr2:172936693-172938111\t428\t528\tmy_read\t60\t+\t-\n" > data/expected.bedpe
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t60\t+\t128\t228\t255,0,0\t1\t100\t0\nchr2:172936693-172938111\t428\t528\tmy_read/2\t60\t-\t428\t528\t255,0,0\t1\t100\t0\n" > data/expected.bed12
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t0\t+\nchr2:172936693-172938111\t428\t528\tmy_read/2\t0\t-\n" > data/expected_ed.bed
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t60\t+\t128\t228\t255,0,0\t1\t100\t0\nchr2:172936693-172938111\t428\t528\tmy_read/2\t60\t-\t428\t528\t255,0,0\t1\t100\t0\n" > data/expected_color.bed12
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t60\t+\t100M\nchr2:172936693-172938111\t428\t528\tmy_read/2\t60\t-\t100M\n" > data/expected_cigar.bed

ls data/

# Test 1: 
mkdir test1
cd test1

echo "> Run bedtools bamtobed on BAM file"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output.bed" \

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../data/expected.bed"
echo "- test1 succeeded -"

cd ..

# Test 2:
mkdir test2
cd test2

echo "> Run bedtools bamtobed on BAM file with -bedpe"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output.bedpe" \
  --bedpe

# checks
assert_file_exists "output.bedpe"
assert_file_not_empty "output.bedpe"
assert_identical_content "output.bedpe" "../data/expected.bedpe"
echo "- test2 succeeded -"

cd ..

# Test 3:
mkdir test3
cd test3

echo "> Run bedtools bamtobed on BAM file with -bed12"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output.bed12" \
  --bed12

# checks
assert_file_exists "output.bed12"
assert_file_not_empty "output.bed12"
assert_identical_content "output.bed12" "../data/expected.bed12"
echo "- test3 succeeded -"

cd ..

# Test 4:
mkdir test4
cd test4

echo "> Run bedtools bamtobed on BAM file with -ed"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output_ed.bed" \
  --ed

# checks
assert_file_exists "output_ed.bed"
assert_file_not_empty "output_ed.bed"
assert_identical_content "output_ed.bed" "../data/expected_ed.bed"
echo "- test4 succeeded -"

cd ..

# Test 5:
mkdir test5
cd test5

echo "> Run bedtools bamtobed on BAM file with -color"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output_color.bed12" \
  --color

# checks
assert_file_exists "output_color.bed12"
assert_file_not_empty "output_color.bed12"
assert_identical_content "output_color.bed12" "../data/expected_color.bed12"
echo "- test5 succeeded -"

cd ..

# Test 6:
mkdir test6
cd test6

echo "> Run bedtools bamtobed on BAM file with -cigar"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output_cigar.bed" \
  --cigar

# checks
assert_file_exists "output_cigar.bed"
assert_file_not_empty "output_cigar.bed"
assert_identical_content "output_cigar.bed" "../data/expected_cigar.bed"
echo "- test6 succeeded -"

cd ..

echo "---- All tests succeeded! ----"
exit 0
