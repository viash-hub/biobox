#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

test_dir="${meta_resources_dir}/test_data"

# create temporary directory and clean up on exit
TMPDIR=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXX")
function clean_up {
 [[ -d "$TMPDIR" ]] && rm -rf "$TMPDIR"
}
trap clean_up EXIT


echo "> Run $meta_name with test data"
"$meta_executable" \
  --gff "$test_dir/1.gff" \
  --output "$TMPDIR/output.txt" 
  
echo ">> Checking output"
[ ! -f "$TMPDIR/output.txt" ] && echo "Output file output.txt does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$TMPDIR/output.txt" ] && echo "Output file output.txt is empty" && exit 1

echo ">> Check if output matches expected output"
diff "$TMPDIR/output.txt" "$test_dir/agat_sq_stat_basic_1.gff"
if [ $? -ne 0 ]; then
  echo "Output file output.txt does not match expected output"
  exit 1
fi

echo "> Test successful"