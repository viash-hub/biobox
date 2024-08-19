#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

test_dir="${meta_resources_dir}/test_data"
out_dir="${meta_resources_dir}/out_data"

echo "> Run $meta_name with test data"
"$meta_executable" \
  --gff "$test_dir/1.gff" \
  --output "$out_dir/output.tsv" 

echo ">> Checking output"
[ ! -f "$out_dir/output.tsv" ] && echo "Output file output.tsv does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$out_dir/output.tsv" ] && echo "Output file output.tsv is empty" && exit 1

echo ">> Check if output matches expected output"
diff "$out_dir/output.tsv" "$test_dir/agat_convert_sp_gff2tsv_1.tsv"
if [ $? -ne 0 ]; then
  echo "Output file output.tsv does not match expected output"
  exit 1
fi

echo "> Test successful"