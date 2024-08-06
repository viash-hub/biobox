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

# Create and populate example.bed
cat << EOF > test_data/example.bed
# Header
chr21	9719758	9729320	variant1	chr21	9719768	9721892	ALR/Alpha	1004	+
chr21	9719758	9729320	variant1	chr21	9721905	9725582	ALR/Alpha	1010	+
chr21	9719758	9729320	variant1	chr21	9725582	9725977	L1PA3	3288	+
chr21	9719758	9729320	variant1	chr21	9726021	9729309	ALR/Alpha	1051	+
chr21	9729310	9757478	variant2	chr21	9729320	9729809	L1PA3	3897	-
chr21	9729310	9757478	variant2	chr21	9729809	9730866	L1P1	8367	+
chr21	9729310	9757478	variant2	chr21	9730866	9734026	ALR/Alpha	1036	-
chr21	9729310	9757478	variant2	chr21	9734037	9757471	ALR/Alpha	1182	-
chr21	9795588	9796685	variant3	chr21	9795589	9795713	(GAATG)n	308	+
chr21	9795588	9796685	variant3	chr21	9795736	9795894	(GAATG)n	683	+
chr21	9795588	9796685	variant3	chr21	9795911	9796007	(GAATG)n	345	+
chr21	9795588	9796685	variant3	chr21	9796028	9796187	(GAATG)n	756	+
chr21	9795588	9796685	variant3	chr21	9796202	9796615	(GAATG)n	891	+
chr21	9795588	9796685	variant3	chr21	9796637	9796824	(GAATG)n	621	+
EOF

# Create and populate expected output files for different tests
cat << EOF > test_data/expected.bed
chr21	9719758	9729320	6353
chr21	9729310	9757478	14482
chr21	9795588	9796685	3604
EOF
cat << EOF > test_data/expected_max.bed
chr21	9719758	9729320	variant1	3288
chr21	9729310	9757478	variant2	8367
chr21	9795588	9796685	variant3	891
EOF

# Test 1: without operation option, default operation is sum
mkdir test1
cd test1

echo "> Run bedtools groupby on BED file"
"$meta_executable" \
  --input "../test_data/example.bed" \
  --groupby "1,2,3" \
  --column "9" \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected.bed"
echo "- test1 succeeded -"

cd ..

# Test 2: with operation option
mkdir test2
cd test2

echo "> Run bedtools groupby on BED file with max operation"
"$meta_executable" \
  --input "../test_data/example.bed" \
  --groupby "1-4" \
  --column "9" \
  --operation "max" \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../test_data/expected_max.bed"
echo "- test2 succeeded -"

cd ..

# Test 3: 






echo "---- All tests succeeded! ----"
exit 0
