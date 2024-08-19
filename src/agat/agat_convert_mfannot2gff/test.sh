#!/bin/bash

set -eo pipefail

## VIASH START
## VIASH END

test_dir="${meta_resources_dir}/test_data"

# create temporary directory and clean up on exit
TMPDIR=$(mktemp --tmpdir "$meta_temp_dir")
function clean_up {
  rm -rf "$TMPDIR"
}
trap clean_up EXIT

echo "> Run $meta_name with test data"
"$meta_executable" \
  --mfannot "$test_dir/test.mfannot" \
  --gff "$TMPDIR/output.gff" 

echo ">> Checking output"
[ ! -f "$TMPDIR/output.gff" ] && echo "Output file output.gff does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$TMPDIR/output.gff" ] && echo "Output file output.gff is empty" && exit 1

echo ">> Check if output matches expected output"
diff "$TMPDIR/output.gff" "$test_dir/agat_convert_mfannot2gff_1.gff"
if [ $? -ne 0 ]; then
  echo "Output file output.gff does not match expected output"
  exit 1
fi

echo "> Test successful"