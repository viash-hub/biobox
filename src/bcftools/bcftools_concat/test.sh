#!/bin/bash

## VIASH START
## VIASH END

# Exit on error
set -eo pipefail

#test_data="$meta_resources_dir/test_data"

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
##fileformat=VCFv4.1
##contig=<ID=1,length=249250621,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	752567	llama	G	C,A	.	.	.	.	1/2
1	752722	.	G	A,AAA	.	.	.	.	./.
EOF

bgzip -c $TMPDIR/example.vcf > $TMPDIR/example.vcf.gz
tabix -p vcf $TMPDIR/example.vcf.gz

cat <<EOF > "$TMPDIR/example_2.vcf"
##fileformat=VCFv4.1
##contig=<ID=1,length=249250621,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	752569	cat	G	C,A	.	.	.	.	1/2
1	752739	.	G	A,AAA	.	.	.	.	./.
EOF

# Test 1: Default test
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "> Run bcftools_concat default test"
"$meta_executable" \
  --input "../example.vcf" \
  --input "../example_2.vcf" \
  --output "concatenated.vcf" \

#  &> /dev/null

# checks
assert_file_exists "concatenated.vcf"
assert_file_not_empty "concatenated.vcf"
assert_file_contains "concatenated.vcf" "concat -o concatenated.vcf ../example.vcf ../example_2.vcf"
echo "- test1 succeeded -"

popd > /dev/null

echo "---- All tests succeeded! ----"
exit 0

echo
cat concatenated.vcf
echo

