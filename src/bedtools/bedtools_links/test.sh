#!/bin/bash

# exit on error
set -eo pipefail

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

# Create test data
cat <<EOF > "$TMPDIR/genes.bed"
chr21	9928613	10012791	uc002yip.1	0	-
chr21	9928613	10012791	uc002yiq.1	0	-
chr21	9928613	10012791	uc002yir.1	0	-
chr21	9928613	10012791	uc010gkv.1	0	-
chr21	9928613	10061300	uc002yis.1	0	-
chr21	10042683	10120796	uc002yit.1	0	-
chr21	10042683	10120808	uc002yiu.1	0	-
chr21	10079666	10120808	uc002yiv.1	0	-
chr21	10080031	10081687	uc002yiw.1	0	-
chr21	10081660	10120796	uc002yix.2	0	-
EOF

# Test 1: Default Use
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "> Run bedtools_links on BED file"
"$meta_executable" \
  --input "../genes.bed" \
  --output "genes.html"

# checks
assert_file_exists "genes.html"
assert_file_not_empty "genes.html"
assert_file_contains "genes.html" "uc002yip.1"
echo "- test1 succeeded -"

popd > /dev/null

# Test 2: BED12 file
mkdir "$TMPDIR/test2" && pushd "$TMPDIR/test2" > /dev/null

echo "> Run bedtools_links with base option"
"$meta_executable" \
  --input "../genes.bed" \
  --output "genes.html" \
  --base_url "http://genome.ucsc.edu"

# checks
assert_file_exists "genes.html"
assert_file_not_empty "genes.html"
assert_file_contains "genes.html" "uc002yip.1"
echo "- test2 succeeded -"

popd > /dev/null

# Test 3: Uncompressed BAM file
mkdir "$TMPDIR/test3" && pushd "$TMPDIR/test3" > /dev/null

echo "> Run bedtools_links with organism option and genome database build"
"$meta_executable" \
  --input "../genes.bed" \
  --output "genes.html" \
  --base_url "http://genome.ucsc.edu" \
  --organism "mouse" \
  --database "mm9"

# checks
assert_file_exists "genes.html"
assert_file_not_empty "genes.html"
assert_file_contains "genes.html" "uc002yip.1"
echo "- test3 succeeded -"

popd > /dev/null

echo "---- All tests succeeded! ----"
exit 0
