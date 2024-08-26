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

# Create example BED12 file
cat <<EOF > "$TMPDIR/example.bed12"
chr21	10079666	10120808	uc002yiv.1	0	-	10081686	1	0	1	2	0	6	0	8	0	4	528,91,101,215,	0,1930,39750,40927,
chr21	10080031	10081687	uc002yiw.1	0	-	10080031	1	0	0	8	0	0	3	1	0	2	200,91,	0,1565,
chr21	10081660	10120796	uc002yix.2	0	-	10081660	1	0	0	8	1	6	6	0	0	3	27,101,223,	0,37756,38913,
EOF

# Expected output bed6 file
cat <<EOF > "$TMPDIR/expected.bed6"
chr21	10079666	10120808	uc002yiv.1	0	-
chr21	10080031	10081687	uc002yiw.1	0	-
chr21	10081660	10120796	uc002yix.2	0	-
EOF
# Expected output bed6 file
cat <<EOF > "$TMPDIR/expected_n.bed6"
chr21	10079666	10120808	uc002yiv.1	1	-
chr21	10080031	10081687	uc002yiw.1	1	-
chr21	10081660	10120796	uc002yix.2	1	-
EOF

# Test 1: Default conversion BED12 to BED6
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "> Run bedtools_bed12tobed6 on BED12 file"
"$meta_executable" \
  --input "../example.bed12" \
  --output "output.bed6"

# checks
assert_file_exists "output.bed6"
assert_file_not_empty "output.bed6"
assert_identical_content "output.bed6" "../expected.bed6"
echo "- test1 succeeded -"

popd > /dev/null

# Test 2: Conversion BED12 to BED6 with -n option
mkdir "$TMPDIR/test2" && pushd "$TMPDIR/test2" > /dev/null

echo "> Run bedtools_bed12tobed6 on BED12 file with -n option"
"$meta_executable" \
  --input "../example.bed12" \
  --output "output.bed6" \
  --n_score

# checks
assert_file_exists "output.bed6"
assert_file_not_empty "output.bed6"
assert_identical_content "output.bed6" "../expected_n.bed6"
echo "- test2 succeeded -"

popd > /dev/null

echo "---- All tests succeeded! ----"
exit 0
