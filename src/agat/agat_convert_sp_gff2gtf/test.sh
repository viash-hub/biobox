#!/bin/bash

## VIASH START
## VIASH END

test_dir="${meta_resources_dir}/test_data"

echo "> Run $meta_name with test data"
"$meta_executable" \
  --gff "$test_dir/0_test.gff" \
  --output "output.gtf" 

echo ">> Checking output"
[ ! -f "output.gtf" ] && echo "Output file output.gtf does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "output.gtf" ] && echo "Output file output.gtf is empty" && exit 1

rm output.gtf

echo "> Test successful"