#!/bin/bash

# exit on error
set -eo pipefail

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

echo "Creating Test Data..."
TMPDIR=$(mktemp -d "$meta_temp_dir/XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# Generate expected files for comparison
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t60\t+\nchr2:172936693-172938111\t428\t528\tmy_read/2\t60\t-\n" > "$TMPDIR/expected.bed"
printf "chr2:172936693-172938111\t128\t228\tchr2:172936693-172938111\t428\t528\tmy_read\t60\t+\t-\n" > "$TMPDIR/expected.bedpe"
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t60\t+\t128\t228\t255,0,0\t1\t100\t0\nchr2:172936693-172938111\t428\t528\tmy_read/2\t60\t-\t428\t528\t255,0,0\t1\t100\t0\n" > "$TMPDIR/expected.bed12"
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t0\t+\nchr2:172936693-172938111\t428\t528\tmy_read/2\t0\t-\n" > "$TMPDIR/expected_ed.bed"
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t60\t+\t128\t228\t250,250,250\t1\t100\t0\nchr2:172936693-172938111\t428\t528\tmy_read/2\t60\t-\t428\t528\t250,250,250\t1\t100\t0\n" > "$TMPDIR/expected_color.bed12"
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t60\t+\t100M\nchr2:172936693-172938111\t428\t528\tmy_read/2\t60\t-\t100M\n" > "$TMPDIR/expected_cigar.bed"
printf "chr2:172936693-172938111\t128\t228\tmy_read/1\t85\t+\nchr2:172936693-172938111\t428\t528\tmy_read/2\t85\t-\n" > "$TMPDIR/expected_tag.bed"


# Test 1: 
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "> Run bedtools bamtobed on BAM file"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output.bed" \

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected.bed"
echo "- test1 succeeded -"

popd > /dev/null

# Test 2:
mkdir "$TMPDIR/test2" && pushd "$TMPDIR/test2" > /dev/null

echo "> Run bedtools bamtobed on BAM file with -bedpe"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output.bedpe" \
  --bedpe

# checks
assert_file_exists "output.bedpe"
assert_file_not_empty "output.bedpe"
assert_identical_content "output.bedpe" "../expected.bedpe"
echo "- test2 succeeded -"

popd > /dev/null

# Test 3:
mkdir "$TMPDIR/test3" && pushd "$TMPDIR/test3" > /dev/null

echo "> Run bedtools bamtobed on BAM file with -bed12"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output.bed12" \
  --bed12

# checks
assert_file_exists "output.bed12"
assert_file_not_empty "output.bed12"
assert_identical_content "output.bed12" "../expected.bed12"
echo "- test3 succeeded -"

popd > /dev/null

# Test 4:
mkdir "$TMPDIR/test4" && pushd "$TMPDIR/test4" > /dev/null

echo "> Run bedtools bamtobed on BAM file with -ed"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output_ed.bed" \
  --edit_distance

# checks
assert_file_exists "output_ed.bed"
assert_file_not_empty "output_ed.bed"
assert_identical_content "output_ed.bed" "../expected_ed.bed"
echo "- test4 succeeded -"

popd > /dev/null

# Test 5:
mkdir "$TMPDIR/test5" && pushd "$TMPDIR/test5" > /dev/null

echo "> Run bedtools bamtobed on BAM file with -color"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output_color.bed12" \
  --bed12 \
  --color "250,250,250" \
  
# checks
assert_file_exists "output_color.bed12"
assert_file_not_empty "output_color.bed12"
assert_identical_content "output_color.bed12" "../expected_color.bed12"
echo "- test5 succeeded -"

popd > /dev/null

# Test 6:
mkdir "$TMPDIR/test6" && pushd "$TMPDIR/test6" > /dev/null

echo "> Run bedtools bamtobed on BAM file with -cigar"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output_cigar.bed" \
  --cigar

# checks
assert_file_exists "output_cigar.bed"
assert_file_not_empty "output_cigar.bed"
assert_identical_content "output_cigar.bed" "../expected_cigar.bed"
echo "- test6 succeeded -"

popd > /dev/null

# Test 7:
mkdir "$TMPDIR/test7" && pushd "$TMPDIR/test7" > /dev/null

echo "> Run bedtools bamtobed on BAM file with -tag"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output_tag.bed" \
  --tag "XT"

# checks
assert_file_exists "output_tag.bed"
assert_file_not_empty "output_tag.bed"
assert_identical_content "output_tag.bed" "../expected_tag.bed"
echo "- test7 succeeded -"

popd > /dev/null

# Test 8: 
mkdir "$TMPDIR/test8" && pushd "$TMPDIR/test8" > /dev/null

echo "> Run bedtools bamtobed on BAM file with other options"
"$meta_executable" \
  --input "$test_data/example.bam" \
  --output "output.bed" \
  --bedpe \
  --mate1 \
  --split \
  --splitD \

# checks
assert_file_exists "output.bed"
assert_file_not_empty "output.bed"
assert_identical_content "output.bed" "../expected.bedpe"
echo "- test8 succeeded -"

popd > /dev/null

echo "---- All tests succeeded! ----"
exit 0
