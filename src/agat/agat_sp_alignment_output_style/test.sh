#!/bin/bash

## VIASH START
## VIASH END

test_dir="${meta_resources_dir}/test_data"
out_dir="${meta_resources_dir}/out_data"

echo "> Run $meta_name with test data"
"$meta_executable" \
  --gff "$test_dir/0.gff" \
  --output "$out_dir/output.gff" 

echo ">> Checking output"
[ ! -f "$out_dir/output.gff" ] && echo "Output file output.gff does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$out_dir/output.gff" ] && echo "Output file output.gff is empty" && exit 1

echo ">> Check if output matches expected output"
diff "$out_dir/output.gff" "$test_dir/agat_sp_alignment_output_style_1.gff
.gff"
if [ $? -ne 0 ]; then
  echo "Output file output.gff does not match expected output"
  exit 1
fi

echo "> Test successful"