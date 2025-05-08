#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

test_dir="${meta_resources_dir}/test_data"

# create temporary directory and clean up on exit
TMPDIR=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXX")
function clean_up {
 [[ -d "$TMPDIR" ]] && rm -rf "$TMPDIR"
}
trap clean_up EXIT

echo "> Run $meta_name with test data"
"$meta_executable" \
  --gff "$test_dir/1.gff" \
  --attribute protein_id \
  --output "$TMPDIR/output.txt" 

echo ">> Checking output"
[ ! -f "$TMPDIR/output_protein_id.txt" ] && echo "Output file output_protein_id.txt does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$TMPDIR/output_protein_id.txt" ] && echo "Output file output_protein_id.txt is empty" && exit 1

echo ">> Check if output matches expected output"
diff "$TMPDIR/output_protein_id.txt" "$test_dir/agat_sp_extract_attributes_1.txt"
if [ $? -ne 0 ]; then
  echo "Output file output_protein_id.txt does not match expected output"
  exit 1
fi

echo "> Test successful"