#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
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
cat <<EOF > "$TMPDIR/example.vcf"
EOF

# Create expected output
cat <<EOF > "$TMPDIR/expected_output.vcf"
EOF

# Test 1: Default Use
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "> Run bcftools_sort on VCF file"
"$meta_executable" \
  --input "../example.vcf" \
  --output "output.vcf"

# checks
assert_file_exists "output.vcf"
assert_file_not_empty "output.vcf"
assert_identical_content "output.vcf" "../expected_output.vcf"
echo "- test1 succeeded -"

popd > /dev/null

# Test 2: Output type   


# Test 3: Output type


echo "---- All tests succeeded! ----"
exit 0
