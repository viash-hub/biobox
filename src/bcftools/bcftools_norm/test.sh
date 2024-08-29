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
1	752567	llama	A	C	.	.	.	.	.
1	752722	.	G	A	.	.	.	.	.
EOF

bgzip -c $TMPDIR/example.vcf > $TMPDIR/example.vcf.gz
tabix -p vcf $TMPDIR/example.vcf.gz

# Test 1: Remove ID annotations
mkdir "$TMPDIR/test1" && pushd "$TMPDIR/test1" > /dev/null

echo "> Run bcftools_norm"
"$meta_executable" \
  --input "../example.vcf" \
  --output "normalized.vcf" \

# checks
assert_file_exists "normalized.vcf"
assert_file_not_empty "normalized.vcf"
assert_file_contains "normalized.vcf" "##fileformat=VCFv4.2"
echo "- test1 succeeded -"

popd > /dev/null


echo "---- All tests succeeded! ----"
exit 0


