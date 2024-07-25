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

# Create and populate featuresA.bed
printf "chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n" > "test_data/featuresA.bed"

# Create and populate featuresB.bed
printf "chr1\t180\t280\nchr1\t290\t390\nchr1\t500\t600\n" > "test_data/featuresB.bed"

# Create and populate expected output files for different tests
printf "chr1\t180\t200\nchr1\t180\t250\nchr1\t300\t390\n" > "test_data/expected_default.bed"
printf "chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n" > "test_data/expected_wa.bed"
printf "chr1\t180\t200\tchr1\t180\t280\nchr1\t180\t250\tchr1\t180\t280\nchr1\t300\t390\tchr1\t290\t390\n" > "test_data/expected_wb.bed"
printf "chr1\t100\t200\tchr1\t180\t280\nchr1\t150\t250\tchr1\t180\t280\nchr1\t300\t400\tchr1\t290\t390\n" > "test_data/expected_loj.bed"
printf "chr1\t100\t200\tchr1\t180\t280\t20\nchr1\t150\t250\tchr1\t180\t280\t70\nchr1\t300\t400\tchr1\t290\t390\t90\n" > "test_data/expected_wo.bed"
printf "chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n" > "test_data/expected_u.bed"
printf "chr1\t100\t200\t1\nchr1\t150\t250\t1\nchr1\t300\t400\t1\n" > "test_data/expected_c.bed"
printf "chr1\t180\t250\nchr1\t300\t390\n" > "test_data/expected_f50.bed"
printf "chr1\t180\t250\nchr1\t300\t390\n" > "test_data/expected_f30.bed"
printf "chr1\t180\t200\nchr1\t180\t250\nchr1\t300\t390\n" > "test_data/expected_f10.bed"
printf "chr1\t180\t200\nchr1\t180\t250\nchr1\t300\t390\n" > "test_data/expected_r.bed"

# Test 1: -wa option
mkdir test1
cd test1

echo "> Run bedtools_intersect on BED files with default intersect"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_default.bed"
echo "- test1 succeeded -"

cd ..

# Test 2: Default intersect
mkdir test2
cd test2

echo "> Run bedtools_intersect on BED files with -wa option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --write_a true

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_wa.bed"
echo "- test2 succeeded -"

cd ..

# Test 3: -wb option
mkdir test3
cd test3

echo "> Run bedtools_intersect on BED files with -wb option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --write_b true

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_wb.bed"
echo "- test3 succeeded -"

cd ..

# Test 4: -loj option
mkdir test4
cd test4

echo "> Run bedtools_intersect on BED files with -loj option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --left_outer_join true

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_loj.bed"
echo "- test4 succeeded -"

cd ..

# Test 5: -wo option
mkdir test5
cd test5

echo "> Run bedtools_intersect on BED files with -wo option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --write_A_and_B true

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_wo.bed"
echo "- test5 succeeded -"

cd ..

# Test 6: -u option
mkdir test6
cd test6

echo "> Run bedtools_intersect on BED files with -u option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --report_A_if_no_overlap true

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_u.bed"
echo "- test6 succeeded -"

cd ..

# Test 7: -c option
mkdir test7
cd test7

echo "> Run bedtools_intersect on BED files with -c option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --number_of_overlaps_A true

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_c.bed"
echo "- test7 succeeded -"

cd ..

# Test 8: -f 0.50 option
mkdir test8
cd test8

echo "> Run bedtools_intersect on BED files with -f 0.50 option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --min_overlap_A 0.50

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_f50.bed"
echo "- test8 succeeded -"

cd ..

# Test 9: -f 0.30 option
mkdir test9
cd test9

echo "> Run bedtools_intersect on BED files with -f 0.30 option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --min_overlap_A 0.30

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_f30.bed"
echo "- test9 succeeded -"

cd ..

# Test 10: -f 0.10 option
mkdir test10
cd test10

echo "> Run bedtools_intersect on BED files with -f 0.10 option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --min_overlap_A 0.10

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_f10.bed"
echo "- test10 succeeded -"

cd ..

# Test 11: -r option
mkdir test11
cd test11

echo "> Run bedtools_intersect on BED files with -r option"
"$meta_executable" \
  --input_a "../test_data/featuresA.bed" \
  --input_b "../test_data/featuresB.bed" \
  --output "output.bed" \
  --reciprocal_overlap true

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_r.bed"
echo "- test11 succeeded -"

cd ..
