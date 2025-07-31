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
  --ref "$test_dir/25_test.gff" \
  --add "$test_dir/9_test.gff" \
  --output "$TMPDIR/output.gff" 

echo ">> Checking output"
[ ! -f "$TMPDIR/output.gff" ] && echo "Output file output.gff does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$TMPDIR/output.gff" ] && echo "Output file output.gff is empty" && exit 1

echo ">> Check if output matches expected output"
diff "$TMPDIR/output.gff" "$test_dir/agat_sp_complement_annotations_1.gff"
if [ $? -ne 0 ]; then
  echo "Output file output.gff does not match expected output"
  exit 1
fi

rm -rf "$TMPDIR/output.gff"

echo "> Run $meta_name with test data"
"$meta_executable" \
  --ref "$test_dir/agat_sp_complement_annotations_ref.gff" \
  --add "$test_dir/agat_sp_complement_annotations_add.gff" \
  --output "$TMPDIR/output.gff" 

echo ">> Checking output"
[ ! -f "$TMPDIR/output.gff" ] && echo "Output file output.gff does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$TMPDIR/output.gff" ] && echo "Output file output.gff is empty" && exit 1

echo ">> Check if output matches expected output"
diff "$TMPDIR/output.gff" "$test_dir/agat_sp_complement_annotations_2.gff"
if [ $? -ne 0 ]; then
  echo "Output file output.gff does not match expected output"
  exit 1
fi
echo "> Test successful"