#!/bin/bash

## VIASH START
## VIASH END

test_dir="${meta_resources_dir}/test_data"

# create temporary directory and clean up on exit
TMPDIR=$(mktemp --tmpdir "$meta_temp_dir")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

echo "> Run $meta_name with test data 1"
"$meta_executable" \
  --gff "$test_dir/file1.gff;$test_dir/file2.gff" \
  --output "$TMPDIR/output.gff"

echo ">> Checking output"
[ ! -f "$TMPDIR/output.gff" ] && echo "Output file output.gff does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$TMPDIR/output.gff" ] && echo "Output file output.gff is empty" && exit 1

echo ">> Check if output matches expected output"
diff "$TMPDIR/output.gff" "$test_dir/agat_sp_merge_annotations_1.gff"
if [ $? -ne 0 ]; then
  echo "Output file output.gff does not match expected output"
  exit 1
fi

echo ">> cleanup"
rm -rf "$TMPDIR/output.gff"

echo "> Run $meta_name with test data 2"
"$meta_executable" \
  --gff "$test_dir/fileA.gff;$test_dir/fileB.gff" \
  --output "$TMPDIR/output.gff"

echo ">> Checking output"
[ ! -f "$TMPDIR/output.gff" ] && echo "Output file output.gff does not exist" && exit 1

echo ">> Check if output is empty"
[ ! -s "$TMPDIR/output.gff" ] && echo "Output file output.gff is empty" && exit 1

echo ">> Check if output matches expected output"
diff "$TMPDIR/output.gff" "$test_dir/agat_sp_merge_annotations_2.gff"
if [ $? -ne 0 ]; then
  echo "Output file output.gff does not match expected output"
  exit 1
fi

echo "> Test successful"