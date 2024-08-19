#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

test_dir="${meta_resources_dir}/test_data"
out_dir="${meta_resources_dir}/out_data"

echo "> Run $meta_name with test data"
"$meta_executable" \
  --gff "$test_dir/1.gff" \
  --attribute protein_id \
  --output "$out_dir/output.txt" 

echo ">> Checking output"
[ ! -f "$out_dir/output_protein_id.txt" ] && echo "Output file output_protein_id.txt does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$out_dir/output_protein_id.txt" ] && echo "Output file output_protein_id.txt is empty" && exit 1

echo ">> Check if output matches expected output"
diff "$out_dir/output_protein_id.txt" "$test_dir/agat_sp_extract_attributes_1.txt"
if [ $? -ne 0 ]; then
  echo "Output file output_protein_id.txt does not match expected output"
  exit 1
fi

echo "> Test successful"