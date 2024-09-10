#!/bin/bash

# exit on error
set -eo pipefail

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
TMPDIR=$(mktemp -d "$meta_temp_dir/XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# Create and populate example.bed
cat << EOF > $TMPDIR/example.bed
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
cat << EOF > $TMPDIR/expected.bed
chr21	9719758	9729320	6353
chr21	9729310	9757478	14482
chr21	9795588	9796685	3604
EOF
cat << EOF > $TMPDIR/expected_max.bed
chr21	9719758	9729320	variant1	3288
chr21	9729310	9757478	variant2	8367
chr21	9795588	9796685	variant3	891
EOF
cat << EOF > $TMPDIR/expected_full.bed
chr21	9719758	9729320	variant1	chr21	9719768	9721892	ALR/Alpha	1004	+	6353
chr21	9729310	9757478	variant2	chr21	9729320	9729809	L1PA3	3897	-	14482
chr21	9795588	9796685	variant3	chr21	9795589	9795713	(GAATG)n	308	+	3604
EOF
cat << EOF > $TMPDIR/expected_delimited.bed
chr21	9719758	9729320	variant1	1004;1010;3288;1051
chr21	9729310	9757478	variant2	3897;8367;1036;1182
chr21	9795588	9796685	variant3	308;683;345;756;891;621
EOF
cat << EOF > $TMPDIR/expected_precision.bed
chr21	9719758	9729320	variant1	1.6e+03
chr21	9729310	9757478	variant2	3.6e+03
chr21	9795588	9796685	variant3	6e+02
EOF

# Test 1: without operation option, default operation is sum
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "> Run bedtools groupby on BED file"
"$meta_executable" \
  --input "../example.bed" \
  --groupby "1,2,3" \
  --column "9" \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected.bed"
echo "- test1 succeeded -"

popd > /dev/null

# Test 2: with operation max option
mkdir "$TMPDIR/test2" && pushd "$TMPDIR/test2" > /dev/null

echo "> Run bedtools groupby on BED file with max operation"
"$meta_executable" \
  --input "../example.bed" \
  --groupby "1-4" \
  --column "9" \
  --operation "max" \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected_max.bed"
echo "- test2 succeeded -"

popd > /dev/null

# Test 3: full option
mkdir "$TMPDIR/test3" && pushd "$TMPDIR/test3" > /dev/null

echo "> Run bedtools groupby on BED file with full option"
"$meta_executable" \
  --input "../example.bed" \
  --groupby "1-4" \
  --column "9" \
  --full \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected_full.bed"
echo "- test3 succeeded -"

popd > /dev/null

# Test 4: header option
mkdir "$TMPDIR/test4" && pushd "$TMPDIR/test4" > /dev/null

echo "> Run bedtools groupby on BED file with header option"
"$meta_executable" \
  --input "../example.bed" \
  --groupby "1-4" \
  --column "9" \
  --header \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_file_contains "output.bed" "# Header"
echo "- test4 succeeded -"

popd > /dev/null

# Test 5: Delimiter and collapse
mkdir "$TMPDIR/test5" && pushd "$TMPDIR/test5" > /dev/null

echo "> Run bedtools groupby on BED file with delimiter and collapse options"
"$meta_executable" \
  --input "../example.bed" \
  --groupby "1-4" \
  --column "9" \
  --operation "collapse" \
  --delimiter ";" \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected_delimited.bed"
echo "- test5 succeeded -"

popd > /dev/null

# Test 6: precision option
mkdir "$TMPDIR/test6" && pushd "$TMPDIR/test6" > /dev/null

echo "> Run bedtools groupby on BED file with precision option"
"$meta_executable" \
  --input "../example.bed" \
  --groupby "1-4" \
  --column "9" \
  --operation "mean" \
  --precision 2 \
  --output "output.bed"

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected_precision.bed"
echo "- test6 succeeded -"

popd > /dev/null

echo "---- All tests succeeded! ----"
exit 0
