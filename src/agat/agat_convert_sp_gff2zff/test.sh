#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

test_dir="${meta_resources_dir}/test_data"

echo "> Run $meta_name with test data"
"$meta_executable" \
  --gff "$test_dir/1.gff" \
  --fasta "$test_dir/1.fa" \
  --output "output" \
  --output_dir "output_dir"

echo ">> Checking output"
[ ! -f "output_dir/output.ann" ] && echo "Output file output.ann does not exist" && exit 1
[ ! -f "output_dir/output.dna" ] && echo "Output file output.dna does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "output_dir/output.ann" ] && echo "Output file output.ann is empty" && exit 1
[ ! -s "output_dir/output.dna" ] && echo "Output file output.dna is empty" && exit 1

echo ">> Check if output matches expected output"
diff "output_dir/output.ann" "$test_dir/agat_convert_sp_gff2zff_1.gff"
if [ $? -ne 0 ]; then
  echo "Output file output.ann does not match expected output"
  exit 1
fi

echo "> Test successful"