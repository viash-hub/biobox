#!/bin/bash

set -eo pipefail

test_dir="${meta_resources_dir}/test_data"

# create temporary directory and clean up on exit
TMPDIR=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXX")
function clean_up {
 [[ -d "$TMPDIR" ]] && rm -rf "$TMPDIR"
}
trap clean_up EXIT

cd "$TMPDIR"

mkdir test1
pushd test1

echo "> Run $meta_name with test data and --emblmygff3"
"$meta_executable" \
  --gff "$test_dir/1.gff" \
  --output "output.txt" \

echo ">> Checking output"
[ ! -f "output.txt" ] && echo "Output file output.txt does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "output.txt" ] && echo "Output file output.txt is empty" && exit 1

echo ">> Check if output matches expected output"
diff "output.txt" "$test_dir/agat_sp_statistics_1.txt"
if [ $? -ne 0 ]; then
  echo "Output file output.txt does not match expected output"
  exit 1
fi

echo "> Test successful"


popd
mkdir test2
pushd test2

cat <<EOF > genome.fasta
>sample_sequence
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
EOF

echo "> Run $meta_name with both gs_size and gs_fasta"
error_message=$("$meta_executable" \
  --gff "$test_dir/1.gff" \
  --output "output.txt" \
  --gs_size "1000000" \
  --gs_fasta "genome.fasta" 2>&1 || true)

expected_error="[error] Please provide only one of the following options to set genome size: --gs_size or --gs_fasta"
if [[ "$error_message" != *"$expected_error"* ]]; then
  echo "Output error message: $error_message does not match expected error message: $expected_error"
  exit 1
fi

echo "> Error test successful"

echo "---- All tests succeeded! ----"
exit 0